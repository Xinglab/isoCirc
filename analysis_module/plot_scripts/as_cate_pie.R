library(ggplot2)
library(dplyr)
library(ggpubr)

argv<-commandArgs(trailingOnly = TRUE)

group<-c('BSJ', 'FSJ', 'Both')
data <- data.frame(
  cate=group,
  value=c(3920, 5913, 72769),
  labels=c("3,920", "5,913", "72,769")
)
data$cate <- factor(data$cate, levels =group, labels= group, ordered = TRUE)

# Compute the position of labels
data <- data %>% 
  arrange(desc(group)) %>%
  mutate(prop = value / sum(data$value) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

fig <- argv[1]

blank_theme <- theme_minimal()+
     theme(
         axis.title.x = element_blank(),
         axis.title.y = element_blank(),
         axis.text.x=element_blank(),
         panel.border = element_blank(),
         panel.grid=element_blank(),
         axis.ticks = element_blank(),
         plot.title=element_blank(), 
         legend.text = element_text(size = 36, color="black", family="sans"), 
         legend.title = element_blank(),
         legend.position = "bottom"
     )

ggplot(data, aes(x="", y=value, fill=cate)) +
    geom_bar(stat="identity", width=0.5) +
    blank_theme + 
    coord_polar("y", start=0, direction= -1) + 
    scale_fill_brewer(palette="Set1") + 
    guides(fill = guide_legend(ncol=3))
ggsave(fig, plot=last_plot(), dpi=300, width=7, height=6)

