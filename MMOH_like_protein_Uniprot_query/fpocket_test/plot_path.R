library(dplyr)
library(ggplot2)

link <- "C:\\Users\\lisap\\Lab_code\\MMOH_like_protein_Uniprot_query\\fpocket_test\\path_tracing.txt"

table = read.table(link, sep = ' ', header = TRUE)

table %>% 
  ggplot(aes(x = pdb, y = path_length)) + 
  geom_bar(aes(fill = as.factor(round(min_sphere, 2)), color = pdb), position = "dodge", stat="identity") + 
  scale_fill_brewer(palette = "Spectral") +  
  scale_color_manual(values=c("4GAM"="black","1MTY"="black", "6YD0"="black", "6YDI"="black", "6VK4"="black"), guide='none') + 
  theme_bw() + 
  xlab('PDB code') + 
  ylab('Number of pockets \nconnecting active site to bulk solvent') + 
  guides(fill=guide_legend(title="Alpha sphere minimum radii")) +
  theme(axis.title.y = element_text(size = 14, face = 'bold')) + 
  theme(axis.title.x = element_text(size = 14, face = 'bold')) + 
  theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12))
  


  