library(ggplot2)
library(RColorBrewer)

argv<-commandArgs(trailingOnly = TRUE)

data<-read.table(argv[1], header=T)
fig<-argv[2]

# theme
theme_jack <- function (base_size = 20, base_family = "") {
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
high="red"
mid="yellow"
low="white"
midpoint=2200#1500
max=4400#3000

# font size
title_size=25
font_size=20
font_color="black"
label_size=4.5
font_family="sans"
plot_title=""#Pairwise comparison across 12 human tissues"
legend_title='Count' #"Isoform count"
# xlab="Number of replicates"
# ylab="Percentage"

raw_names=c('Adipose', 'Adrenal', 'Blood', 'Brain', 'Heart', 'Kidney', 'Liver', 'Lung', 'Prostate', 'SkeletalMuscle', 'SmoothMuscle', 'Testis')
names=c('Adipose', 'Adrenal', 'Blood', 'Brain', 'Heart', 'Kidney', 'Liver', 'Lung', 'Prostate', 'Skeletal Muscle', 'Smooth Muscle', 'Testis')
# removed Pancreas, Colon, Thyroid
data$Var1<- factor(data$Var1, levels = raw_names, labels=names, ordered=TRUE)
data$Var2<- factor(data$Var2, levels = raw_names, labels=names, ordered=TRUE)

# data$Cnt <- factor(data$Cnt, levels = c('1', '2', '3'), labels=c('Read count >= 1', 'Read count >= 2', 'Read count >= 3'))

ls=c("0", "1,100", "2,200", "3,300", "4,400")
bs=seq(0, 4400, 1100)

data <- data[which(data$Cate== "Isoform"),]
ggplot(data, aes(x=Var1, y=Var2)) +
	geom_tile(aes(fill=Value)) +
	geom_text(aes(label = format(as.numeric(Value), big.mark=",", scientific=F)), size=label_size) +
	scale_colour_gradient2(low = low, high = high, space = "Lab", mid=mid, midpoint=midpoint, na.value = "grey50", guide = "colourbar", aesthetics = "fill", limits=c(0,max), labels=ls, breaks=bs) +
	# scale_fill_continuous(high = "red", low = "red") + 
	# scale_fill_distiller(palette = "Reds") +
    # ggtitle(plot_title) +
    xlab(NULL) + 
    ylab(NULL) +
    labs(fill=legend_title) + 
	facet_wrap(. ~ Cate, scales = "free", nrow=1) +
    theme(text=element_text(size=font_size, family=font_family), strip.text=element_text(hjust=0.5, size=title_size), plot.title = element_text(hjust = 0.5, size=title_size), axis.text.x = element_text(color=font_color, angle=30, hjust=1), axis.text.y = element_text(color=font_color)) # put title in the center
    # axis.text.x = element_text(size=font_size, family=font_family, angle=45, hjust=1)) # put title in the center
ggsave(fig, plot=last_plot(), dpi=300, height=8, width=10)
# ggsave(fig, plot=last_plot(), dpi=300, height=8, width=16)