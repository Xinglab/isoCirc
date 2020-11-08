library(ggplot2)

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
color_palatte<-"Set1" # http://mkweb.bcgsc.ca/brewer/swatches/brewer-palettes-swatches.pdf
# font size
title_size=20
font_size=14
font_color="black"
font_family="sans"
text_font_size=8
plot_title="Error rate of raw read and consensus sequence with different copy numbers" #Fraction of BSJs"
legend_title="" #Category"
xlab=""
ylab="Error rate"
raw_name<-c('Rep1_1', 'Rep1_2', 'Rep1_3', 'Rep2_1', 'Rep2_2', 'Rep2_3')
data$Sample <- factor(data$Sample, levels = raw_name, 	labels = raw_name, ordered=T)
data$Type <- factor(data$Type, levels = c('Raw', '2', '3', '4', '5', '6', '11'), labels=c('Raw read', '2', '3', '4', '5', '6-10', '> 10'), ordered=T)

scal_f<-function(x){x} 
ys = seq(0,0.3,by=0.05)
ylims = c(0, 0.3)

ggplot(data, aes(x=Type, y=Error)) +
    geom_violin() +
    geom_boxplot(width=0.1) +  
    # geom_boxplot(show.legend = FALSE)+
    # coord_flip() +  
    scale_y_continuous(labels = scales::percent_format(accuracy=1), limits=ylims, breaks=ys) + 
    # scale_y_continuous(breaks=ys, limits=ylims) +  
    ggtitle(plot_title) +
    xlab(xlab) +
    ylab(ylab) +
    labs(fill=legend_title) + 
    facet_wrap(. ~ Sample, scales = "free", nrow=3) +
    theme(text=element_text(size=font_size, family=font_family), 
          plot.title = element_text(hjust = 0.5, size=title_size), axis.text=element_text(size=font_size, color="black"),
          axis.text.x = element_text(color=font_color, size=font_size, family=font_family),
          legend.position= "bottom") +
    # + scale_fill_brewer(palette = color_palatte)
ggsave(fig, plot=last_plot(), dpi=300, width=14, height=9)
# ggsave(fig, plot=last_plot(), dpi=300, width=12, height=5)