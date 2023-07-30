library(ggplot2)
library(dplyr)

plt <- ggplot(data) + 
       geom_histogram(aes(x = V1), stat='count')

theme_set(theme_minimal())

data %>% 
  count(V1, sort=TRUE) %>%
  filter(n >= 1000 ) %>%
  ggplot() + 
  geom_bar(aes(x = reorder(V1, n), y = n, fill = log(n)), stat='identity') + 
  xlab('Representative structure with n >= 1000') + 
  ylab('Membership count') + 
  theme(axis.text.x = element_text(hjust=1,angle = 45)) +
  theme(axis.title.y = element_text(vjust = 2)) + 
  theme(axis.title.x = element_text(vjust=2))
  
  
data %>% 
  count(V1) %>%
  ggplot() + 
  geom_bar(aes(x = n)) + 
  ylim(0, 20) + 
  xlim(0, 500) + 
  xlab('Cluster member count')

