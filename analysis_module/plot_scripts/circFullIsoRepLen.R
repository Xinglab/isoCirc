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
font_size=18
font_family="sans"
font_color="black"
plot_title=""#Cumulative distribution of transcript length for full-length circRNA isoforms"
# legend_title="Read count bins"
legend_title="BSJ-FSJ category"
xlab="Transcript length"
ylab="Cumulative fraction"

cates=c('FSM/NIC-FSM/NIC', 'FSM/NIC-NNC', 'NNC-FSM/NIC', 'NNC-NNC', 'All-All')
cnts=c()
for (cate in cates) {
	sub_data<-data[which(data$Cate==cate),]
	c = format(length(sub_data$Cate), big.mark=",", scientific=F)
	cnts=append(cnts, c)
}
new_cates=paste0(cates, ' (', cnts, ')')

# data$Rep <- factor(data$Rep, levels=c("1", "2", "3"))
# data$Sample<- factor(data$Sample, levels = c('rep1', 'rep2'), labels= c('Rep#1', 'Rep#2'), ordered=TRUE)
# data$Cate <- factor(data$Cate, levels=c('FSM-FSM', 'FSM-NIC/NNC', 'NIC/NNC-FSM', 'NIC/NNC-NIC/NNC', 'All-All'), ordered=TRUE)
data$Cate <- factor(data$Cate, levels=cates, labels=new_cates, ordered=TRUE)
# data$CountBin <- factor(data$CountBin, levels=c('1', '2', '3~5','6~10', '11~50', '>50'), ordered=TRUE)

scal_f <- function(x){log10(x)} # x

# data<-data[which(data$Cate=='FSM-FSM'),]
# ggplot(data, aes(scal_f(Len), colour=Rep)) + # # of replicates that have this isoform
ggplot(data, aes(scal_f(Len), colour=Cate)) + # BSJ-Inter catelogs
# ggplot(data, aes(scal_f(Len), colour=CountBin)) + # read count bins
	stat_ecdf(size=line_size) +
	# ggtitle(plot_title) +
	xlab(xlab) +
	ylab(ylab) +
	scale_x_continuous(breaks=c(1, 2, 2.4771212547196626, 2.6989700043360187, 3, 3.6989700043360187), labels=c("10", "100", "300",  "500", "1,000", "5,000"), limits=c(1, 3.6989700043360187)) +
	theme(text=element_text(size=font_size, family=font_family), 
		  plot.title = element_text(hjust = 0.5), 
		  legend.position= c(0.8,0.2), 
		  axis.text=element_text(color=font_color)) + # put title in the center
	# facet_wrap(. ~ Sample, scales = "free", nrow=1) +
	# scale_fill_brewer(palette = color_palatte, name=legend_title, labels=c("SampA"=samp1, "SampB"=samp2)) +
	scale_color_manual(values=c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#000000"), name=legend_title)
	# scale_color_brewer(palette = color_palatte, name=legend_title)

ggsave(fig, plot=last_plot(), dpi=300, width=9.3, height=6)