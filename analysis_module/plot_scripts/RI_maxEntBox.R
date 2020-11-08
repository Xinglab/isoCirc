library(ggpubr)
library(RColorBrewer)

argv<-commandArgs(trailingOnly = TRUE)

data<-read.table(argv[1], header=T)
fig <- argv[2]

# font size
font_size=15
font_family="sans"
font_face="bold"
plot_title="Maxent scores of retained introns and non-retained introns"
legend_title=""
scal_f <- function(x){x} #{log2(x)}
# xlab="log10(read length)"
xlab=""
ylab="Maxent score"
color_palatte<-"Set1" # http://mkweb.bcgsc.ca/brewer/swatches/brewer-palettes-swatches.pdf

# theme
theme_jack <- function (base_size = font_size, base_family = "") {
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

data<-subset(data, Anno=='All')
data$Type <- factor(data$Type, levels = c('5P','3P'), labels= c('5 prime', '3 prime'), ordered = TRUE)
data$Cate <- factor(data$Cate, levels = c('NonIR', 'RI'), labels= c('Non-RI in circRNA', 'RI in circRNA'), ordered = TRUE)

p <- ggboxplot(data, x = "Type",  y = "Scores", color = "Category") + 
	 stat_compare_means(aes(group=Category), label = "p.format", method="wilcox.test") +
	 scale_y_continuous(limits=c(-50,16), breaks=c(-50, -20, -10, 0, 10, 16)) + 
     theme(axis.title.x=element_blank()) +
	 ylab(ylab)  +
     theme(text=element_text(size=font_size, family=font_family), 
			plot.title = element_text(hjust = 0.5),
			legend.position='top') + 
     scale_color_brewer(palette = color_palatte)#, name=legend_title)

ggsave(fig, plot=last_plot(), dpi=300, width=6, height=6)