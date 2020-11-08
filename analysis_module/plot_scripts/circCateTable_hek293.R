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
midpoint=12000
max=24000

# font size
title_size=25
font_size=20
font_color="black"
label_size=8
font_family="sans"
plot_title="BSJ/FSJ Category"
legend_title="Isoform count"


data <- data[which(data$ReadCnt== "2"),]
data$BSJ<- factor(data$BSJ, levels = c('FSM', 'NIC', 'NNC'), labels=c('FSM', 'NIC', 'NNC'), ordered=TRUE)
data$ISO<- factor(data$ISO, levels = c('FSM', 'NIC', 'NNC'), labels=c('FSM', 'NIC', 'NNC'), ordered=TRUE)
data$ReadCnt <- factor(data$ReadCnt, levels = c('1', '2', '3'), labels=c('Read count \u2265 1', 'Read count \u2265 2', 'Read count \u2265 3'))
xlab <- "BSJ"
ylab <- "FSJ"
ls=c("0", "6,000", "12,000", "18,000", "24,000") # HEK293
bs=seq(0, 24000, 6000) # HEK293
ggplot(data, aes(x=BSJ, y=ISO)) +
	geom_tile(aes(fill=IsoCnt)) +
	geom_text(aes(label = format(as.numeric(IsoCnt), big.mark=",", scientific=F)), size=label_size) +
	scale_colour_gradient2(low = low, high = high, space = "Lab", mid=mid, midpoint=midpoint, na.value = "grey50", guide = "colourbar", aesthetics = "fill", limits=c(0,max), labels=ls, breaks=bs) +
	# scale_fill_continuous(high = "red", low = "red") + 
	# scale_fill_distiller(palette = "Reds") +
    ggtitle(plot_title) +
    xlab(xlab) + 
    ylab(ylab) +
    labs(fill=legend_title) + 
	facet_wrap(. ~ ReadCnt, scales = "free", nrow=3) +
    theme(text=element_text(size=font_size, family=font_family), strip.text=element_text(hjust=0.5, size=title_size), plot.title = element_text(hjust = 0.5, size=title_size), axis.text.x = element_text(color=font_color, size=font_size, family=font_family, angle=90, hjust=1), axis.text.y = element_text(color=font_color, size=font_size, family=font_family)) # put title in the center
ggsave(fig, plot=last_plot(), dpi=300, height=7, width=8)