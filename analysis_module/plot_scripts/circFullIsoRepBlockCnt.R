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
plot_title="" #Cumulative distribution of exon number for full-length circRNA isoforms"
# legend_title="Read count bins"
legend_title="BSJ-FSJ category"
xlab="Exon number"
ylab="Cumulative fraction"

cates=c('FSM/NIC-FSM/NIC', 'FSM/NIC-NNC', 'NNC-FSM/NIC', 'NNC-NNC', 'All-All')
cnts=c()
for (cate in cates) {
	sub_data<-data[which(data$Cate==cate),]
	c = format(length(sub_data$Cate), big.mark=",", scientific=F)
	cnts=append(cnts, c)
}
new_cates=paste0(cates, ' (', cnts, ')')

data$Cate <- factor(data$Cate, levels=cates, labels=new_cates, ordered=TRUE)

scal_f <- function(x){x} 

xlim=25

ggplot(data, aes(scal_f(ExonNum), colour=Cate)) +
	stat_ecdf(size=line_size) +
	# ggtitle(plot_title) +
	xlab(xlab) +
	ylab(ylab) +
	xlim(c(0,xlim)) + 
	guides(color=guide_legend(ncol=1)) + 
	theme(text=element_text(size=font_size, family=font_family),
		  plot.title = element_text(hjust = 0.5),
		  legend.position=c(0.8,0.2), 
		  axis.text=element_text(color=font_color)) + # put title in the center
	scale_color_manual(values=c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#000000"), name=legend_title)
	# scale_color_brewer(palette = color_palatte, name=legend_title)

ggsave(fig, plot=last_plot(), dpi=300, width=9.3, height=6)