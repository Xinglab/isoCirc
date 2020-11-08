library(ggplot2)
library(ggrepel)
library(ggpubr)

argv<-commandArgs(trailingOnly = TRUE)

data<-read.table(argv[1], header=T)
fig<-argv[2]

theme_jack <- function (base_size = 12, base_family = "") {
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
font_size=20
font_family="sans"
font_color="black"
text_font_size=6
plot_title=""
legend_title=""
xlab="log10 normalized BSJ number (short-read)"
ylab="log10 normalized BSJ number (isoCirc long read)"

.labs <- rownames(data)

data$Human_tissue <- factor(data$Human_tissue, 
  levels=c('Brain', 'Testis', 'Kidney', 'Prostate', 'Liver', 'Heart', 'Lung', 'Skeletal_muscle'), 
  labels=c('Brain', 'Testis', 'Kidney', 'Prostate', 'Liver', 'Heart', 'Lung', 'Skeletal Muscle'), order=T)

pearson <- round(cor(log10(data$SR_ratio), log10(data$LR_ratio), method="pearson"), digits=2)
spearman <- round(cor(log10(data$SR_ratio), log10(data$LR_ratio), method="spearman"), digits=2)

annotation <- data.frame(
   x = c(-4.3,-4.3),
   y = c(-2.2,-2.25),
   label = c(paste("Pearson     r = ", pearson), paste('Spearman ρ = ', spearman))
)

ggplot(data, aes(x=log10(SR_ratio), y =log10(LR_ratio))) +
 geom_point(size=3, color="red") + 
 scale_y_continuous(breaks=seq(-3.2, -2.2, 0.25), limits=c(-3.2,-2.2)) +
 scale_x_continuous(breaks=seq(-4.4, -3.2, 0.4), limits=c(-4.4, -3.2)) + 
 xlab(xlab) + ylab(ylab) +
 geom_text(data=annotation, aes(x=x, y=y, label=label, hjust = 0), size=text_font_size, family=font_family, color=font_color, fontface="bold")+
 # geom_text(aes(label=Human_tissue), size = text_font_size, family=font_family) +
 geom_text_repel(aes(label = Human_tissue), size = text_font_size, family=font_family)+
 theme(text=element_text(size=font_size, family=font_family), 
          plot.title = element_text(vjust = 1.0), axis.text=element_text(color=font_color)
          # axis.text.x = element_text(size=font_size, color=font_color, family=font_family, angle=0, hjust=0)
          ) 

ggsave(fig, plot=last_plot(), dpi=300, width=8, height=7)