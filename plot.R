'''
PLOTS FOR THE BIFURCATED INTERACTIONS PROJECT

'''
library(ggplot2)
library(ggridges)
# ---------------------------------
# CHAIN COUNTS
df <- read.csv('/home/user/revathy/bifur/programs/new_code/test_files/chain_count.csv', header = TRUE)
df$count <- as.factor(df$count)
df$chain_num <- as.factor(df$chain_num)

ggplot(df, aes(x = chain_num, y = count)) +
  geom_segment(aes(x = chain_num, xend = chain_num, y = 0, yend = count), color="grey") +
  scale_x_continuous(limits = c(3, 40)) + 
  scale_y_continuous(limits = c(1, 300)) +
  geom_point( color="orange", size=1) +
  theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  xlab("Number of chains") +
  ylab("Count")


ggplot(df, aes(x = chain_num, y = count)) +
  geom_density_ridges()

ggplot(df, aes(x = chain_num)) + geom_histogram(bins = 50) + theme_classic() 

# make the bubble plot - coloured by helix pair
ggplot(df, aes(x = chain_num, y = count)) + 
  geom_point(aes(size = dist, fill = helix), alpha = 0.75, shape = 21) + 
  scale_size_continuous(limits = c(0.0, 5.0), range = c(0, 20), breaks = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5)) + 
  labs(x = "Dimer type", y = "Helix pair", size = "Inter-helix distance", fill = "") + 
  theme(legend.key = element_blank(), 
        axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 90, vjust = 0.3, hjust = 1), 
        axis.text.y = element_text(colour = "black", face = "bold", size = 11), 
        legend.text = element_text(size = 10, face ="bold", colour ="black"), 
        legend.title = element_text(size = 12, face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2), legend.position = "right") +  
  scale_fill_manual(values = colours, guide = FALSE) + 
  scale_y_discrete(limits = rev(levels(df2$helix)))





library(ggplot2)

df <- read.csv('/home/user/revathy/bifur/Programs/new_code/test_files/nofilter_rbi_int_type.csv', header = TRUE)
ggplot(df, aes(x=inttype, y=relval)) +
  geom_bar(stat="identity", position = 'dodge') + 
  scale_x_discrete(limits=c('van der Waals interaction', 'hydrophobic interaction', 'hydrogen bond', 'salt-bridge', 'cation-pi interaction',  'anion-pi interaction', 'pi-pi interaction', 'amino-pi interaction', 'pi-sulphur bond', 'disulphide bond')) +
  scale_fill_brewer(palette="Paired") + 
  labs(title = "RBI interaction types", xlab = "Interaction type", ylab = "Normalised frequency") + 
  coord_flip() + 
  theme_classic()




# DUMBBELL PLOT
# create data sets
df <- read.csv('/home/user/revathy/bifur/Programs/new_code/test_files/hs_foldx_inter.csv', header = TRUE)
ylabel <- c('GLY', 'ALA', 'VAL', 'LEU', 'MET', 'ILE', 'PHE', 'TYR', 'TRP', 'SER', 'THR', 'CYS', 'PRO', 'ASN', 'GLN', 'LYS', 'ARG', 'HIS', 'ASP', 'GLU')

x1 <- 
x2 <- c(4,3,6,7,5,4)

datamain <- data.frame(ylabel,x1,x2)

# import ggplot2, ggalt and tidyverse
library(ggplot2)  
library(ggalt)    
library(tidyverse)

# Draw dumbbell plot
ggplot() +
  geom_dumbbell(data = datamain, aes(y = ylabel,
                                     x = x1, 
                                     xend = x2),
                size = 1.5, color = "blue", size_x = 7,
                size_xend = 7, colour_x = "green",
                colour_xend = "yellow")