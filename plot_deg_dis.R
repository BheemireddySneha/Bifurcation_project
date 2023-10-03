library(ggplot2)

df1 <- read.csv('/home/revathy/lab/bifur/bifur/new_code/interface_files/results/interface_deg.csv')
ggplot(df1, aes(x=res, y=Degree)) +
  geom_violin(width=2.1, size=0.2) +
  theme(
    legend.position="none"
  ) +
  coord_flip() + # This switch X and Y axis and allows to get the horizontal version
  xlab("") +
  ylab("Assigned Probability (%)")

df1 <- read.csv('/home/revathy/lab/bifur/bifur/new_code/interface_files/results/interface_deg.csv')
ggplot(df1, aes(x=trbi, y=Degree)) +
  geom_violin() +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("A boxplot with jitter") +
  xlab("")

transform(df1, Degree = as.numeric(Degree))
df1$Degree = as.numeric(as.character(df1$Degree))
ggplot(df1, aes(x=res, y=Degree, fill=res)) +
  geom_violin(width=1.4) +
  geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  scale_x_discrete(limits=rev(c('GLY', 'ALA', 'VAL', 'LEU', 'MET', 'ILE', 'PHE', 'TYR', 'TRP', 'SER', 'THR', 'CYS', 'PRO', 'ASN', 'GLN', 'LYS', 'ARG', 'HIS', 'ASP', 'GLU'))) +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("A Violin wrapping a boxplot") +
  xlab("") + 
  coord_flip()
