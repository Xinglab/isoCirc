library(ggplot2)
argv<-commandArgs(trailingOnly = TRUE)

data<-read.table(argv[1], header=T)
fig <- argv[2]

# font size
font_size=15
font_family="sans"
font_face="bold"
font_color="black"
plot_title="Maxent scores of non-retained introns and retained introns"
line_size=1

legend_title=""
scal_f <- function(x){x} #{log2(x)}
# xlab="log10(read length)"
xlab="Maxent score"
ylab="Cumulative fraction"
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
data$Cate <- factor(data$Cate, levels = c('NonRI', 'RI'), labels= c('Non-RI in circRNA', 'RI in circRNA'), ordered = TRUE)
# data$Category <- factor(data$Category, levels = c('Anno','NonIR', 'RI'), labels= c('Annotated RI in linear transcript', 'Non-RI in circRNA', 'RI in circRNA'), ordered = TRUE)

# data<-data[which(data$RepCnt == '2'),]
# data$Known<- factor(data$Known, levels = c('Known','Novel'),  ordered = TRUE)
# data$RepCnt <- factor(data$RepCnt, levels = c('1', '2'), ordered = TRUE)
# data$Type <- factor(data$Type, levels = c('5P','3P'), labels= c('5 prime', '3 prime'), ordered = TRUE)

# ggplot(data, aes(scal_f(Scores), colour=Known)) +
ggplot(data, aes(scal_f(Scores), colour=Cate)) +
	stat_ecdf(size=line_size) +
	ggtitle(plot_title) +
	xlab(xlab) +
	ylab(ylab) +
	theme(
		plot.title = element_text(hjust = 0.5), legend.title.align=0.5, 
		text=element_text(size=font_size, color=font_color, family=font_family), 
		axis.text.x = element_text(size=font_size, color=font_color, family=font_family), 
		axis.text.y=element_text(size=font_size, color=font_color, family=font_family),
		legend.position= "bottom") + # put title in the center
		scale_x_continuous(limits=c(-50, 16), breaks=c(-50, -20, -10, 0, 10, 16)) + # for including 0s
	facet_wrap(. ~ Type, scales = "free", nrow=1) +
	# scale_fill_brewer(palette = color_palatte, name=legend_title, labels=c("SampA"=samp1, "SampB"=samp2)) +
	scale_color_brewer(palette = color_palatte, name=legend_title)

ggsave(fig, plot=last_plot(), dpi=300, width=12, height=6)