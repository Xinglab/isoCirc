library(ggplot2)
library(RColorBrewer)

argv<-commandArgs(trailingOnly = TRUE)

data<-read.table(argv[1], header=T)
fig<-argv[2]

# theme
theme_jack <- function (base_size = 20, base_family = "") {
 theme_gray(base_size = base_size, base_family = base_family) %+replace%
   theme(
     panel.background = element_rect(fill=NULL),
     #plot.background = element_rect(fill=NULL,size = 0.5,linetype = “dash”,colour = “black”),# element_line(linetype = “dash”),
     panel.grid.minor.y = element_line(size=0),
     panel.grid.minor.x = element_line(size=0),
     panel.grid.major = element_line(colour = "grey80",size=.2)#,
     #axis.line = element_line(colour = NULL)
   )
}
theme_set(theme_jack())

# parameters
# color
high="black"
mid="red"
low="yellow"
midpoint=0.5

# font size
title_size=20
font_size=9
text_font_size=4
font_color="black"
font_family="sans"
plot_title="Pairwise comparison"
legend_title="Degree of\nsimilarity"

data$Var1<- factor(data$Var1, levels = c('nano_43', 'nano_44', 'nano_45', 'nano_58', 'nano_59', 'nano_60', 'Illumina_1', 'Illumina_2', 'Illumina_3'), 
    labels= c('Rep1_1', 'Rep1_2', 'Rep1_3', 'Rep2_1', 'Rep2_2', 'Rep2_3', 'Illumina_1', 'Illumina_2', 'Illumina_3'), ordered=TRUE)
data$Var2<- factor(data$Var2, levels = c('nano_43', 'nano_44', 'nano_45', 'nano_58', 'nano_59', 'nano_60', 'Illumina_1', 'Illumina_2', 'Illumina_3'), 
    labels= c('Rep1_1', 'Rep1_2', 'Rep1_3', 'Rep2_1', 'Rep2_2', 'Rep2_3', 'Illumina_1', 'Illumina_2', 'Illumina_3'), ordered=TRUE)
# data<-data[which((data$Var1 == 'Illumina_1' | data$Var1 == 'Illumina_2' | data$Var1 == 'Illumina_3') & 
    # (data$Var2 == 'Illumina_1' | data$Var2 == 'Illumina_2' | data$Var2 == 'Illumina_3')),]


data$Cnt <- factor(data$Cnt, levels = c('1', '2', '3'), labels=c('Read count \u2265 1', 'Read count \u2265 2', 'Read count \u2265 3'))

ggplot(data, aes(x=Var1, y=Var2)) +
    geom_tile(aes(fill=Value)) +
    geom_text(aes(label=Value),size = text_font_size) + 
    scale_colour_gradient2(low = low, high = high, space = "Lab", 
    mid=mid, midpoint=midpoint,
    na.value = "grey50", guide = "colourbar", aesthetics = "fill",
    limits=c(0,1)) +
    # scale_fill_continuous(high = "red", low = "red") + 
    # scale_fill_distiller(palette = "Reds") +
    ggtitle(plot_title) +
    xlab(NULL) + 
    ylab(NULL) +
    labs(fill=legend_title) + 
    facet_wrap(. ~ Cnt, scales = "free", nrow=3) +
    theme(legend.title.align=0.5, text=element_text(size=font_size, family=font_family), strip.text=element_text(hjust=0.5, size=title_size), plot.title = element_text(hjust = 0.5, size=title_size),
        axis.text.x = element_text(size=font_size, color=font_color, family=font_family), #, angle=90, hjust=1), 
        axis.text.y=element_text(size=font_size, color=font_color, family=font_family)) # put title in the center
# ggsave(fig, plot=last_plot(), dpi=300, height=23, width=9)
ggsave(fig, plot=last_plot(), dpi=300, height=12, width=5.2)