# Violin Plot
library(tidyverse)
library(viridis)
combine %>%
  ggplot( aes(x=gr, y=combine, fill=gr)) +
    geom_violin() +
    scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
    theme_ipsum() +
    theme(
      legend.position="none",
      plot.title = element_text(size=11)
    ) +
    ggtitle("Violin chart") +
    xlab("")

data1 <- data %>% 
  gather(key="text", value="value") %>%
  mutate(text = gsub("\\.", " ",text)) %>%
  mutate(value = round(as.numeric(value),0)) %>%
  filter(text %in% c("B1_SPU","B1_LAV","B2_SPU","B2_LAV","B3_SPU", "B3_LAV", "B4_SPU", "B4_LAV", "B5_SPU","B5_LAV"))

# Plot
p <- data1 %>%
  mutate(text = fct_reorder(text, value)) %>% # Reorder data
  ggplot( aes(x=text, y=value, fill=text, color=text)) +
    geom_violin(width=2.1, size=0.2) +
    scale_fill_viridis(discrete=TRUE) +
    scale_color_viridis(discrete=TRUE) +
    theme_ipsum() +
    theme(
      legend.position="none"
    ) +
    coord_flip() + # This switch X and Y axis and allows to get the horizontal version
    xlab("") +
    ylab("")

p

