# Library
require(tidyverse)
require(ggtext)
require(vegan)
require(ggpubr)
require(reshape2)
# Determine seed, theme and colors used in the analysis 
set.seed(1996)
# Data
# Rarefied data
otu_rare = readRDS('data/otu_rare.RDS')
otu_rel <- apply(decostand(otu_rare %>% column_to_rownames('Group') %>% t(), method="total", MARGIN=2),1, mean)
# Import taxonomy table 
taxtab = readRDS('data/taxtab.RDS')
# Import metadata 
metadata = readRDS('data/metadata.RDS')

# Colors 
colm = c('#47B66A')
cols = c('#D9A534')

# Figure 1 
# In which 'biota' was an OTU present and what is its taxonomy+relative abundance?
otutabM = otu_rare %>% filter(str_detect(Group, '^M')) %>%
  column_to_rownames('Group') %>% t()
otu_PA <- 1*((otutabM>0)==1)
otu_occM <- rowSums(otu_PA)/ncol(otu_PA) 

otutabS = otu_rare %>% filter(str_detect(Group, '^S')) %>%
  column_to_rownames('Group') %>% t()
otu_PA <- 1*((otutabS>0)==1)
otu_occS <- rowSums(otu_PA)/ncol(otu_PA) 

otu_where = data.frame(otu_occM, otu_occS) %>%
  mutate(where = ifelse(otu_occS == 0 & otu_occM > 0, 'Microbiota', 
                        ifelse(otu_occM == 0 & otu_occS > 0, 'Sporobiota', 
                               ifelse(otu_occM > 0 & otu_occS > 0, 'Both', 'None'))))
   
tax_occ = rownames_to_column(otu_where, 'name') %>%
  left_join(taxtab, by='name') %>%
  mutate(Phylum = str_replace(Phylum, 
                             '(.*)_unclassified', 'Unclassified *\\1*'), 
         Phylum = str_replace(Phylum, 
                             '^(\\S*)$', '*\\1*'))

ggplot(tax_occ, aes(x=where, fill=Phylum)) +
  geom_bar(stat='count') +
  theme_bw(base_size = 18) +
  theme(axis.title.x = element_blank(), 
        legend.title = element_blank(), 
        legend.text = element_markdown(), 
        legend.key.size = unit(18, 'pt'), 
        text = element_text(family = "Calibri")
        ) +
  labs(y='Number of OTUs')
ggsave('plots/where_taxonomy.png', dpi=600)

# Figure 2
tax_occ_rel = data.frame(tax_occ, otu_rel)
ggplot(tax_occ_rel, aes(x=where, y=otu_rel, fill=where)) +
  geom_boxplot() +
  scale_y_log10() +
  scale_fill_manual(values = c('#1dbeda', colm , cols)) +
  theme_bw(base_size = 18) +
  theme(axis.title.x = element_blank(), 
        legend.position = 'none', 
         text = element_text(family = "Calibri")) +
  labs(y='Relative abundance of OTUs')
ggsave('plots/where_relabund.png', dpi=600) 

# Figure 4 
otuPA = otu_rare %>% column_to_rownames('Group')
otuPA[otuPA > 0] = 1
otuPA = rownames_to_column(otuPA, 'Group')
otuPA_meta <- left_join(otuPA, metadata, by=join_by('Group' == 'samples'))

# Graph that takes into account that some OTUs are present than not and than back again - so all that are completely new
new_otus = otuPA_meta %>% 
  filter(time_point < 13) %>%
    # Select only the columns I need and trasform into longer format
  select(person, time_point, biota, starts_with('Otu')) %>%
  pivot_longer(names_to = 'name', values_to = 'PA', cols = starts_with('Otu')) %>%
  # Group the dataframe by person and otu (OTUs)
  group_by(person, name, biota) %>% 
  # Arrange by day
  arrange(time_point, .by_group = TRUE) %>% 
  # Create new column otu_sum is 1 if the OTU is present (PA > 0) on the current day and was not present on any of the previous days
   # If otu_sum is 1 or more than 1, that means that OTU was present on this day and days before
  # If otu_sum is more than 1, it means it was present in the provious days, so turn that into 0 
  mutate(otu_sum = cumsum(PA), 
         new_otu = ifelse(otu_sum == 1 & lag(otu_sum, default = 0) == 0, 1, 0)) %>%
  #filter(date != min(date)) %>%
  ungroup() %>%
  group_by(person, time_point, biota) %>%
  # percentage of new OTUs in 1 day, for each person
  summarise(., new= sum(new_otu)) %>%
  ungroup() %>%
  group_by(time_point, biota) %>%
  mutate(mean =median(new), 
         sd=sd(new)) %>%
  ungroup()

ggplot(new_otus, aes(x=time_point)) +
  geom_point(aes(y=new, color=biota),size=3) +
  geom_line(mapping=aes(y=mean, color=biota), linewidth=1.5) +
  #geom_ribbon(mapping= aes(ymin=mean-sd, ymax=mean +sd),  fill ='grey',  alpha=.2) +
  scale_color_manual(values = c(colm, cols)) +
  scale_x_continuous(breaks = seq(1, 14)) +
  labs(x='Sampling point', y='Number of new OTUs', color='Type of sample') +
  theme_bw(base_size=18) +
  theme( text = element_text(family = "Calibri"))
ggsave('plots/newOTUs.png', dpi=600)

# Is the correlation linear? 
corr_new = new_otus %>% filter(biota == 'Microbiota') %>% left_join(new_otus %>%filter(biota == 'Sporobiota'), by=join_by('person' == 'person', 'time_point' == 'time_point'))

ggscatter(corr_new, x='new.x', y='new.y', 
          add='reg.line', conf.int = TRUE, 
          cor.coef = TRUE, cor.method = 'pearson', 
          xlab='Acquisition of new OTUs in microbiota samples', ylab= 'Acquisition of new OTUs in sporobiota samples')

# Is data normaly distributed? 
shapiro.test(corr_new$new.y)
# Q-Q plot 
ggqqplot(corr_new$new.x, ylab='Acquisition of new OTUs in microbiota samples')
ggqqplot(corr_new$new.y, ylab='Acquisition of new OTUs in sporobiota samples')

# Person's Corelation of aquisition of new OTUs in Microbiota and Sporobiota 
cor.test(corr_new$new.x, corr_new$time_point, method = 'pearson')
cor.test(corr_new$new.y, corr_new$time_point, method = 'pearson')
# There is significaly significant(p-value = 2e-16), strong positive correlation (0.86) between the acquisition of new OTUs between samples in microbiota and sporobiota. 

# Only for me, not for presentation
# What is the taxonomic determination and relative abundance of the OTUs that are new in later time-points?
new_tax = otuPA_meta %>% 
  filter(time_point < 13) %>%
  select(person, time_point, biota, starts_with('Otu')) %>%
  pivot_longer(names_to = 'name', values_to = 'PA', cols = starts_with('Otu')) %>%
  group_by(person, name, biota) %>% 
  arrange(time_point, .by_group = TRUE) %>% 
  mutate(otu_sum = cumsum(PA), 
         new_otu = ifelse(otu_sum == 1 & lag(otu_sum, default = 0) == 0, 1, 0)) %>%
  filter(time_point != 1) %>%
  filter(new_otu == 1) %>%
  left_join(taxtab, by='name') %>%
  left_join(data.frame(otu_rel) %>% rownames_to_column('name'), by='name') %>%
  group_by(Phylum, biota) %>%
  summarise(sum_relabund =sum(otu_rel)*100) %>%
  mutate(Phylum = if_else(sum_relabund < 0.005, 'Other (< 0.5%)', Phylum), 
         Phylum = str_replace(Phylum, 
                             '(.*)_unclassified', 'Unclassified \\1'))

ggplot(new_tax, aes(x=Phylum, y=sum_relabund, fill=Phylum)) +
  geom_col() +
  scale_y_log10() +
  labs(x='', y='Mean relative abundance (%)', fill='') +
  coord_flip() +
  theme_bw() +
    theme(legend.position = 'none', 
          axis.text.x = element_markdown(), 
          legend.key.size = unit(18, 'pt'))
ggsave('plots/new_taxonomy.png', dpi=600)

# What is the total abundance of OTUs that are newly acquired in that sample?
new_rel = otuPA_meta %>% 
  filter(time_point < 13) %>%
  select(person, time_point, biota, starts_with('Otu')) %>%
  pivot_longer(names_to = 'name', values_to = 'PA', cols = starts_with('Otu')) %>%
  group_by(person, name, biota) %>% 
  arrange(time_point, .by_group = TRUE) %>% 
  mutate(otu_sum = cumsum(PA), 
         new_otu = ifelse(otu_sum == 1 & lag(otu_sum, default = 0) == 0, 1, 0)) %>%
  filter(time_point != 1) %>%
  filter(new_otu == 1) %>%
  left_join(data.frame(otu_rel) %>% rownames_to_column('name'), by='name') %>%
  # Sum the total abundance of this OTUs
  group_by(name, otu_rel) %>%
  summarise(sum_otus=sum(new_otu)) %>%
  ungroup()

ggplot(new_rel, aes(x=otu_rel)) +
  geom_histogram() +
  scale_x_log10()+
  #coord_cartesian(xlim = c(0,200)) +
  labs(x='log10(Relative abundance)', y='Number of OTUs') +
  theme_bw()
ggsave('plots/new_relabund_log10.png', dpi=600)


# Figure 5
# Occupancy plot (code by Shade and Stopnisek)
otutab = otu_rare %>% filter(str_detect(Group, '^M')) %>%
  column_to_rownames('Group') %>% t() 

otu_PA <- 1*((otutab>0)==1)                                              # presence-absence data 
otu_occ <- rowSums(otu_PA)/ncol(otu_PA)                                # occupancy calculation
otu_rel <- apply(decostand(otutab, method="total", MARGIN=2),1, mean)     # mean relative abundance
occ_abun <- data.frame(otu_occ=otu_occ, otu_rel=otu_rel) %>%           # combining occupancy and abundance data frame
  rownames_to_column('name') %>%
  left_join(taxtab, by='name')

# Occupancy abundance plot:
plot1 = ggplot(data=occ_abun, aes(x=log10(otu_rel), y=otu_occ)) +
  geom_point(pch=21, fill=colm) +
  labs(x="log10(mean relative abundance)", y="Occupancy") +
  theme_bw(base_size = 18) +
  theme( text = element_text(family = "Calibri"))

# For sporobiota
otutab = otu_rare %>% filter(str_detect(Group, '^S')) %>%
  column_to_rownames('Group') %>% t() 

otu_PA <- 1*((otutab>0)==1)                                             
otu_occ <- rowSums(otu_PA)/ncol(otu_PA)                               
otu_rel <- apply(decostand(otutab, method="total", MARGIN=2),1, mean)     
occ_abun <- data.frame(otu_occ=otu_occ, otu_rel=otu_rel) %>%        
  rownames_to_column('name') %>%
  left_join(taxtab, by='name')

# Occupancy abundance plot:
plot2 = ggplot(data=occ_abun, aes(x=log10(otu_rel), y=otu_occ)) +
  geom_point(pch=21, fill=cols) +
  labs(x="log10(mean relative abundance)", y="Occupancy") +
  theme_bw(base_size = 18) +
  theme( text = element_text(family = "Calibri"))
  #facet_wrap(~Phylum, ncol=3, nrow = 4)

ggarrange(plot1 + rremove("ylab") + rremove("xlab"), 
          plot2 + rremove("ylab"), 
          labels = NULL, 
          ncol=1)
ggsave('plots/occupancy_mean_relabund.png', dpi=600)

# Figure 3
otu_rel <- decostand(otu_rare %>% column_to_rownames('Group') , method="total", MARGIN=1)
  
merged = left_join(metadata, otu_rel %>% rownames_to_column('samples'), by='samples') %>%
  filter(sample_type == 'regular') %>%
  column_to_rownames('samples') %>%
  select(person, biota, starts_with('Otu'))

final <- data.frame()
for (persona in unique(merged$person)) {
  for (bioti in unique(merged$biota)) {
  merged_sub <- merged[merged$biota == bioti & merged$person == persona,]
  merged_sub <- as.data.frame(t(merged_sub[, 3:2996]))
  
  merged_sub2 <- merged_sub
  merged_sub2[merged_sub2>0] <- 1
  merged_sub2$prevalence <- rowSums(merged_sub2)
  merged_sub$person <- persona
  merged_sub$biota <- bioti
  merged_sub$name <- rownames(merged_sub)
  merged_sub2$name <- rownames(merged_sub2)

  melt1 <- melt(merged_sub2, id.vars = c('prevalence', 'name'))
  melt2 <- melt(merged_sub, id.vars = c('person', 'biota', 'name'))
  melt2$prevalence <- melt1$prevalence
    
  final <- rbind(final, melt2)
  }
}

final$count <- 1/12  
# Group_by biota, person and prevalence and summarize count (because each OTU could be present in 12 time points 1/12)
final_agg <- aggregate(count ~ biota + person + prevalence, data = final, FUN = sum)

final_agg_mean = filter(final_agg, prevalence != 0) %>% group_by(biota, prevalence) %>%
  summarise(mean = mean(count), sd=sd(count))

ggplot(final_agg[final_agg$prevalence != 0,], aes(x = prevalence, y = count, color=biota)) +
  geom_point(size=3) +
  geom_line(final_agg_mean, mapping=aes(y=mean, color=biota), linewidth=1.5) +
  scale_color_manual(values = c(colm, cols)) +
  scale_x_continuous(breaks = seq(0,12, by=1)) +
  labs(x='Occupancy by person', y= 'Number of observed new OTUs', color='Type of sample') +
  theme_bw(base_size = 18) +
  theme(text = element_text(family = "Calibri")) +
  coord_flip()

ggsave('plots/occupancy_count.png', dpi=600)
