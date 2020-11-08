# # plot scatter and line plot
# # created by Yan Gao, 06/15/2017
# library(ggplot2)
# library(MASS)
# library(ggthemes)
# library(LSD)
# argv<-commandArgs(trailingOnly = TRUE)

# data<-read.table(argv[1], header=T)
# fig<-argv[2]

# # theme
# theme_jack <- function (base_size = 12, base_family = "") {
#  theme_gray(base_size = base_size, base_family = base_family) %+replace%
#    theme(
#      panel.background = element_rect(fill=NULL),
#      #plot.background = element_rect(fill=NULL,size = 0.5,linetype = “dash”,colour = “black”),# element_line(linetype = “dash”),
#      panel.grid.minor.y = element_line(size=0),
#      panel.grid.minor.x = element_line(size=0),
#      panel.grid.major = element_line(colour = "grey80",size=.2)#,
#      #axis.line = element_line(colour = NULL)
#    )
# }
# theme_set(theme_jack())

# # Get density of points in 2 dimensions.
# # @param x A numeric vector.
# # @param y A numeric vector.
# # @param n Create a square n by n grid to compute density.
# # @return The density within each square.
# get_density <- function(x, y, ...) {
#   # Calculate 2d density over a grid
#   dens <- kde2d(x,y)
#   # create a new data frame of that 2d density grid
#   # (needs checking that I haven't stuffed up the order here of z?)
#   gr <- data.frame(with(dens, expand.grid(x,y)), as.vector(dens$z))
#   names(gr) <- c("xgr", "ygr", "zgr")

#   # Fit a model
#   mod <- loess(zgr~xgr*ygr, data=gr)
#   return(predict(mod, newdata=data.frame(xgr=x, ygr=y)))
# }

# # parameters
# # color
# color_palatte<-"Set1" # http://mkweb.bcgsc.ca/brewer/swatches/brewer-palettes-swatches.pdf
# # font
# font_size=20
# font_family="sans"
# # text
# plot_title="High-confidence BSJ read count in replicates"
# legend_title="Legend"
# xlab="Rep#1 read count"
# ylab="Rep#2 read count"

# scal_f <- function(x){log2(x+1)} # {x}
# data$Count1<-scal_f(data$Count1)
# data$Count2<-scal_f(data$Count2)
# # get density
# s12 <- data[which(data$Sample == 'rep1_rep2'),]
# s12$Density <- get_density(s12$Count1, s12$Count2)
# # s13$Density <- get_density(s13$Count1, s13$Count2)
# # s23$Density <- get_density(s23$Count1, s23$Count2)

# new_data<-rbind(s12) #,s13,s23)
# new_data$Sample<- factor(new_data$Sample, levels = c('rep1_rep2'), labels= c('Rep#1 vs Rep#2'), ordered=TRUE)
# scal_f <- function(x){log2(x+1)} # {x}
# # ylim<-15
# # xlim<-15

# ggplot(new_data, aes(x=Count1, y=Count2, color=Density)) + geom_point() +
#   scale_colour_gradientn(colours = colorpalette('heat', 5)) + 
#   ggtitle(plot_title) +
#   xlab(xlab) +
#   ylab(ylab) +
#   # ylim(c(0,ylim)) + 
#   # xlim(c(0,xlim)) + 
#   scale_y_continuous(breaks=c(0, 3.45943, 6.65821, 9.96722, 12.28800), labels=c("0","10", "100", "1,000", "5,000"), limits=c(0,12.28800)) +
#   scale_x_continuous(breaks=c(0, 3.45943, 6.65821, 9.96722, 12.28800), labels=c("0","10", "100", "1,000", "5,000"), limits=c(0,12.28800)) +
#   theme(text=element_text(size=font_size, family=font_family), plot.title = element_text(hjust = 0.5),
#     legend.position= "right") +
#   facet_wrap(. ~ Sample, scales = "free", nrow=1)

# ggsave(fig, plot=last_plot(), dpi=300, width=7, height=6)
library(ggplot2)
library(MASS)
library(ggthemes)
library(LSD)

argv<-commandArgs(trailingOnly = TRUE)
data<-read.table(argv[1], header=T)
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

# Get density of points in 2 dimensions.
# @param x A numeric vector.
# @param y A numeric vector.
# @param n Create a square n by n grid to compute density.
# @return The density within each square.
get_density <- function(x, y, b) {
  # Calculate 2d density over a grid
  dens <- kde2d(x,y, 1)
  # create a new data frame of that 2d density grid
  # (needs checking that I haven't stuffed up the order here of z?)
  gr <- data.frame(with(dens, expand.grid(x,y)), as.vector(dens$z))
  names(gr) <- c("xgr", "ygr", "zgr")

  # Fit a model
  mod <- loess(zgr~xgr*ygr, data=gr)
  return(predict(mod, newdata=data.frame(xgr=x, ygr=y)))
}

# parameters
# color
color_palatte<-"Set1" # http://mkweb.bcgsc.ca/brewer/swatches/brewer-palettes-swatches.pdf
# font
font_size=20
title_size=25
font_color="black"
font_family="sans"
# text
plot_title="High-confidence BSJ read count in replicates"
legend_title="Legend"
# xlab="Rep#1 read count"
# ylab="Rep#2 read count"
xlab="" #"Blood" #"Brain"
ylab="" #Testis" #Skeletal Muscle"

# removed Pancreas, Colon, Thyroid
raw_names=c('nano_43', 'nano_44', 'nano_45', 'nano_58', 'nano_59', 'nano_60')
names=c('Nano-43', 'Nano-44', 'Nano-45', 'Nano-58', 'Nano-59', 'Nano-60')
names=c('Rep1_1', 'Rep1_2', 'Rep1_3', 'Rep2_1', 'Rep2_2', 'Rep2_3')
data$Sample1<- factor(data$Sample1, levels = raw_names, labels= names, ordered=TRUE)
data$Sample2<- factor(data$Sample2, levels = raw_names, labels= names, ordered=TRUE)

scal_f <- function(x){log2(x+1)} # {x}
data$Count1<-scal_f(data$Count1)
data$Count2<-scal_f(data$Count2)

# calculate density for each pair
newdata <- data[which(data$Sample1== data$Sample2),]
newdata$Density <- 0
for (name1 in names) {
  for (name2 in names) {
    if (name1 != name2) {
      subdata = data[which(data$Sample1==name1 & data$Sample2 == name2),]
      subdata$Density <- get_density(subdata$Count1, subdata$Count2)
      newdata <- rbind(newdata, subdata)
    }
  }
}
ls=c("0.00", "0.05", "0.10", "0.15", "0.20")
bs=seq(0, 0.2, 0.05)
ggplot(newdata, aes(x=Count1, y=Count2, color=Density)) +
  geom_point(size=0.8) + 
  scale_colour_gradientn(colours = colorpalette('heat', 5), limits= c(-0.01,0.2), breaks=bs, labels=ls) + 
  ggtitle(plot_title) +
  xlab(xlab) +
  ylab(ylab) +
  scale_y_continuous(breaks=c(0, 6.65821, 12.28800), labels=c("0", "100", "5,000"), limits=c(0,12.28800)) +
  scale_x_continuous(breaks=c(0, 6.65821, 12.28800), labels=c("0", "100", "5,000"), limits=c(0,12.28800)) +
  # facet_wrap(. ~ Sample, scales = "free", nrow=1) +
  facet_grid(Sample1 ~ Sample2) +
  theme(text=element_text(size=font_size, family=font_family), plot.title = element_text(hjust = 0.5, size=title_size), legend.position= "right", panel.spacing.x=unit(0.5, "lines") , panel.spacing.y=unit(0.5,"lines"), axis.text.x = element_text(color=font_color, angle=45, hjust=1), axis.text.y = element_text(color=font_color)) 

ggsave(fig, plot=last_plot(), dpi=300, width=10.5, height=9.5, limitsize=FALSE)