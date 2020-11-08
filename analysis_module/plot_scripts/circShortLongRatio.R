library(ggplot2)
library(scales)
argv<-commandArgs(trailingOnly = TRUE)
data<-read.table(argv[1], header=T)
fig<-argv[2]

# 4.00% -> 4%
# Style from Figure 3A. CircRNA isoform count per tissue for 12 tissues.
# Y axis label: Proportion of circRNA read
# Library -> Library type
# Add X axis
# polyA (illumina)
# RNaseR (illumina)
# RNaseR (isoCirc)
# RNaseR (isoCirc)
# IsoCirc -> isoCirc

color="Set1"
line_size=0.5
font_size=16
font_family="sans"
font_color="black"
theme_jack <- function (base_size = 12, base_family = "") {
  theme_gray(base_size = base_size, base_family = base_family) %+replace%
    theme(
      panel.background = element_rect(fill=NULL),
      #plot.background = element_rect(fill=NULL,size = 0.5,linetype = "dash",colour = "black"),# element_line(linetype = "dash"),
      panel.grid.minor.y = element_line(size=0),
      panel.grid.minor.x = element_line(size=0),
      panel.grid.major = element_line(colour = "grey80",size=.2)#,
      #axis.line = element_line(colour = NULL)
    )
}

old_lib<-c('polyA+Illumina', 'RNaseR+Illumina', 'RNaseR+isoCirc:sample1', 'RNaseR+isoCirc:sample2')
new_lib<-c('polyA (Illumina)', 'RNase R (Illumina)', 'RNase R (isoCirc)', 'RNase R (isoCirc)')

theme_set(theme_jack())
dodge <- position_dodge(width=0.9)
ylab <- "Proportion of circRNA reads"

ggplot(data, aes(y=circRNA_read_ratio, x=Library, fill=Library, group=technical_replicate)) +
	geom_bar(size=line_size, color="black", stat = "identity", position = dodge) + # position_dodge(width = 0.95), width=0.8, ) +
	ylab(ylab) +
	scale_y_continuous(labels = scales::percent_format(accuracy=1)) +
    scale_fill_brewer(palette=color, labels=new_lib) + 
    scale_x_discrete(labels=new_lib)+
	theme(legend.position=c(0.2,0.6), legend.text=element_text(size=font_size), text=element_text(size=font_size, color=font_color, family=font_family), axis.text=element_text(size=font_size,color=font_color), 
		axis.title.x=element_blank(), 
		# axis.text.x=element_blank() ,
		axis.ticks.x=element_blank()) 
# ggplot(R6r, aes(y=circRNA_read_ratio,x=Library,group=technical_replicate)) + geom_bar(aes(color=Library,fill=Library),size=line_size, position = dodge, stat = "identity",alpha=.5)+scale_y_continuous(labels = scales::percent)+ylab('BSJ_read_percentage') +theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + scale_color_brewer(palette=color + xlab(xlab))
# ggplot(R6r, aes(y=circRNA_read_ratio,x=Library,group=technical_replicate)) + geom_bar(aes(color=Library),size=line_size, position = dodge, stat = "identity",fill="white")+scale_y_continuous(labels = scales::percent)+ylab('BSJ_read_percentage') +theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + scale_color_brewer(palette=color) + xlab(xlab)


ggsave(fig, plot=last_plot(), dpi=300, width=9, height=7)
