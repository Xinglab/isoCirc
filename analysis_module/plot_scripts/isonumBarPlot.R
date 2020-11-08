# plot Stacked barplot
# created by Yan Gao, 06/14/2017
library(ggplot2)
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
# fill_color="red"
line_color="black"
line_size=0.5
# font size
font_size=20
font_family="sans"
font_color="black"
plot_title=''#Number of full-length isoform of 12 tissues"
legend_title=""
xlab=""
ylab="Number of isoforms"
# ylab="Number of exon-intron-exon isoforms"

data$Sample<- factor(data$Sample, levels = c('Testis', 'Brain', 'Blood', 'Kidney', 'Adrenal', 'Lung', 'SmoothMuscle', 'Adipose', 'Liver', 'Heart', 'SkeletalMuscle', 'Prostate'), labels= c('Testis', 'Brain', 'Blood', 'Kidney', 'Adrenal', 'Lung', 'Smooth Muscle', 'Adipose', 'Liver', 'Heart', 'Skeletal Muscle', 'Prostate'), ordered=TRUE)
# data<-data[which(data$Type=="ISONum1" | data$Type=="ISONum2"),]
# data<-data[which(data$Type=="EIEISONum1" | data$Type=="EIEISONum2"),]

data$Type <- factor(data$Type, levels=c("ISONum1", 'ISONum2', 'ISONum3'), labels=c('Read count = 1', 'Read count = 2', 'Read count \u2265 3'), ordered=TRUE)
# data$Type <- factor(data$Type, levels=c("EIEISONum1", 'EIEISONum2'), labels=c('Read count = 1', 'Read count \u2265 2'), ordered=TRUE)

ys = c(0, 30000, 60000, 90000, 120000)
ylabs= c("0", "30,000", "60,000",  "90,000", "120,000")
ylims = c(0,120000)
# ys = c(0, 600, 1200, 1800)
# ylabs = c("0", "600", "1,200", "1,800")
# ylims = c(0, 1800)

p <- ggplot(data, aes(x=Sample, y=Number, fill=Type)) + 
 	 geom_bar(stat="identity", width=1, size=line_size, color=line_color) +
 	 ggtitle(plot_title) +
 	 xlab(xlab) +
 	 ylab(ylab) +
 	 labs(fill=legend_title) + 
 	 scale_y_continuous(breaks=ys, labels=ylabs, limits=ylims) +
 	 scale_fill_brewer(palette = color_palatte) + 
 	 theme(legend.position=c(0.7,0.7), legend.title = element_blank(), legend.background = element_blank(), text=element_text(size=font_size, family=font_family), 
		plot.title = element_text(hjust = 0.5), axis.text=element_text(color=font_color), axis.text.x = element_text(angle=30, hjust=1), axis.ticks.x=element_blank()) # put title in the center

ggsave(fig, plot=last_plot(), dpi=300, width=8, height=7, limitsize=FALSE)