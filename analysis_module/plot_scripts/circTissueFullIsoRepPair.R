# plot scatter and line plot
# created by Yan Gao, 06/15/2017
library(ggplot2)
library(MASS)
# library(ggthemes)
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
font_size=15
title_size=25
font_family="sans"
font_color="black"
# text
plot_title="Pairwise comparison of circRNA isoform ratios across 12 tissues"
legend_title="Legend"
xlab="" #"Blood" #"Brain"
ylab="" #Testis" #Skeletal Muscle"
# xlab="Brain - log2(read count + 1)"
# ylab="Skeletal Muscle - log2(read count + 1)"

# scal_f <- function(x){x} #{log2(x+1)} # {x}
# data$Count1<-scal_f(data$Ratio1)
# data$Count2<-scal_f(data$Ratio2)
# get density
# data$Density <- get_density(data$Ratio1, data$Ratio2)

# removed Pancreas, Colon, Thyroid
raw_names=c('Adipose', 'Adrenal', 'Blood', 'Brain', 'Heart', 'Kidney', 'Liver', 'Lung', 'Prostate', 'SkeletalMuscle', 'SmoothMuscle', 'Testis')
names=c('Adipose', 'Adrenal', 'Blood', 'Brain', 'Heart', 'Kidney', 'Liver', 'Lung', 'Prostate', 'Skeletal\nMuscle', 'Smooth\nMuscle', 'Testis')
data$Tissue1 <- factor(data$Tissue1, levels = raw_names, labels= names, ordered=TRUE)
data$Tissue2 <- factor(data$Tissue2, levels = raw_names, labels= names, ordered=TRUE)

# calculate density for each pair
newdata <- data[which(data$Tissue1 == data$Tissue2),]
newdata$Density <- 0
for (name1 in names) {
  for (name2 in names) {
    if (name1 != name2) {
      subdata = data[which(data$Tissue1==name1 & data$Tissue2 == name2),]
      subdata$Density <- get_density(subdata$Ratio1, subdata$Ratio2)
      newdata <- rbind(newdata, subdata)
    }
  }
}

ylim<-1.0
xlim<-1.0

ggplot(newdata, aes(x=Ratio1, y=Ratio2, color=Density)) +
  geom_point(size=0.8) +
  scale_colour_gradientn(colours = colorpalette('heat', 5)) + 
  ggtitle(plot_title) +
  xlab(xlab) +
  ylab(ylab) +
  # ylim(c(0,ylim)) + 
  scale_y_continuous(breaks=c(0, 0.5, 1.0), labels=c("0","0.5", "1")) +
  # xlim(c(0,xlim)) + 
  scale_x_continuous(breaks=c(0, 0.5, 1.0), labels=c("0", "0.5", "1")) +
  # facet_wrap(. ~ Sample, scales = "free", nrow=1) +
  facet_grid(Tissue1 ~ Tissue2) +
  theme(text=element_text(size=font_size, family=font_family), plot.title = element_text(hjust = 0.5, size=title_size), legend.position= "right", panel.spacing.x=unit(0.3, "lines") , panel.spacing.y=unit(0.3,"lines"), axis.text=element_text(color=font_color)) 

ggsave(fig, plot=last_plot(), dpi=300, width=13.5, height=12.75, limitsize=FALSE)

