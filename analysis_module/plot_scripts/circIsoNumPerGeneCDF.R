library(ggplot2)
# library(patchwork) # https://github.com/thomasp85/patchwork

argv<-commandArgs(trailingOnly = TRUE)

data<-read.table(argv[1], header=T)
data<-na.omit(data)
fig<-argv[2]

# theme
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
line_size=1
# font
font_size=15
font_family="sans"
font_color="black"
plot_title="Full-length circRNA isoform number per gene"
legend_title="Number of replicates" #
xlab="Isoform number"
ylab="Cumulative fraction"

# iso_cnt_lim=230
# data <- data[which(data$IsoCnt <= iso_cnt_lim),]
data$RepCnt <- factor(data$RepCnt, levels=c("1", "2", "3", "4", "5", "6"))
data$Type <- factor(data$Type, levels = c('Gene', 'BSJ'), labels= c('Gene', 'BSJ'), ordered=TRUE)

data<-data[which(data$Type=="Gene"),]

# ggplot(data, aes(scal_f(IsoCnt), colour=RepCnt)) +
ggplot(data, aes(IsoCnt)) +
	stat_ecdf(size=line_size) +
	ggtitle(plot_title) +
	xlab(xlab) +
	ylab(ylab) +
	# HEK293
	# scale_x_continuous(limits=c(0,90), breaks=c(0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90)) + #, labels=c("0","0.5", "1")) +
	# tissue
	scale_x_continuous(limits=c(0,200), breaks=c(0, 10, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200)) + 
	theme(text=element_text(size=font_size, family=font_family),
		plot.title = element_text(hjust = 0.5),
		legend.position= "bottom", axis.text=element_text(color=font_color)) + # put title in the center
	facet_wrap(. ~ Type, scales = "free", nrow=1) +
	# scale_fill_brewer(palette = color_palatte, name=legend_title, labels=c("SampA"=samp1, "SampB"=samp2)) +
	scale_color_brewer(palette = color_palatte, name=legend_title)

ggsave(fig, plot=last_plot(), dpi=300, width=7, height=6)