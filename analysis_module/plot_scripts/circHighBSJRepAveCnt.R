library(ggplot2)
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
font_color="black"
font_family="sans"
plot_title="BSJ read count"
legend_title=""
xlab="Read count"
ylab="Cumulative fraction"

rep1<-data[which(data$RepAnno == "1Known" | data$RepAnno =="1Novel"),]
rep2<-data[which(data$RepAnno == "2Known" | data$RepAnno =="2Novel"),]
rep2$Count <- (rep2$Count)/2
data<-rbind(rep1, rep2)

data$RepAnno <- factor(data$RepAnno, levels=c("2Known", "2Novel", "1Known", "1Novel"), labels=c('Known, Detected in 2 replicates', 'Novel, Detected in 2 replicates', 'Known, Detected in 1 replicate', 'Novel, Detected in 1 replicate'), ordered=T)
# data$Anno <- factor(data$Anno, levels=c("Known", "Novel"), ordered=T)

scal_f <- function(x){log2(x)} # log2(x)

ggplot(data, aes(scal_f(Count), colour=RepAnno)) +
	stat_ecdf(size=line_size) +
	ggtitle(plot_title) +
	xlab(xlab) +
	ylab(ylab) +
	scale_x_continuous(breaks=c(0, 3.321928094887362, 5.643856189774724, 6.643856189774724, 9.965784284662087, 10.965784284662087), labels=c("1", "10", "50", "100", "1,000", "2,000"), limits=c(0,10.965784284662087)) +
	# scale_x_continuous(breaks=c(0, 3.321928094887362, 5.643856189774724, 6.643856189774724, 9.965784284662087, 10.965784284662087), labels=c("1", "10", "50", "100", "1,000", "2,000"), limits=c(0,10.965784284662087)) +
	# scale_x_continuous(breaks=c(1.584962500721156, 3.321928094887362, 5.643856189774724, 6.643856189774724, 9.965784284662087, 11.965784284662087), labels=c("3", "10", "50", "100", "1,000", "4,000"), limits=c(1.584962500721156,11.965784284662087)) +
	theme(text=element_text(size=font_size, color=font_color, family=font_family),
		plot.title = element_text(hjust = 0.5), legend.title = element_blank(),
		axis.text=element_text(color=font_color),
		legend.position= c(0.7,0.2)) + # put title in the center
	# facet_wrap(. ~ Sample, scales = "free", nrow=1) +
	scale_color_brewer(palette = color_palatte)
	# scale_color_brewer(palette = color_palatte, name=legend_title)

ggsave(fig, plot=last_plot(), dpi=300, width=7, height=6)