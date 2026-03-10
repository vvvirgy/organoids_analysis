library(tidyverse)

?seq

df = tibble(
  x = seq(from = 0, to = 10), 
  y = seq(from = 0, to = 10)
)

df %>% 
  ggplot(aes(
    x, 
    y 
  )) + 
  geom_point(color = 'steelblue4') + 
  geom_line(linetype = 'dashed', colour = 'steelblue4') + 
  theme_bw() + 
  xlim(-0.1,10) + 
  ylim(-0.1,10) + 
  geom_vline(xintercept = 2, color = 'firebrick', linetype = 'dashed') + 
  scale_x_continuous(breaks = c(0,2,4,6,8,10)) + 
  scale_y_continuous(breaks = c(0,2,4,6,8,10)) + 
  labs(
    x = 'ploidy', 
    y = 'expression'
  ) + 
  geom_rect(
    xmin = 0, 
    xmax = 2, 
    ymin = 0, 
    ymax = 2
  )
  

  