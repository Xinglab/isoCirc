library(ggpubr)
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
# xlab="# of replicate"
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

data$Type <- factor(data$Type, levels = c('5P','3P'), labels= c('5 prime', '3 prime'), ordered = TRUE)
data$Known<- factor(data$Known, levels = c('Known','Novel'),  ordered = TRUE)
data$RepCnt <- factor(data$RepCnt, levels = c('1', '2'), ordered = TRUE)
# data$Type <- factor(data$Type, levels = c('5P','3P'), labels= c('5 prime', '3 prime'), ordered = TRUE)

p <- ggboxplot(data, x = "RepCnt",  y = "Scores", color = "Known", palette = "jco") + 
	 stat_compare_means(aes(group = Known), label = "p.format", method="wilcox.test") +
	 ylim(c(-50,20)) + 
	 xlab(xlab) + 
	 ylab(ylab) + 
	 facet_wrap(. ~ Type, scales = "free", nrow=1) +
     # ggtitle(plot_title) +
     theme(text=element_text(size=font_size, family=font_family), 
			plot.title = element_text(hjust = 0.5),
			legend.position='top') + 
     scale_color_brewer(palette = color_palatte, name=legend_title)

ggsave(fig, plot=last_plot(), dpi=300, width=6, height=6)