library(ggplot2)
library(MASS)
library(ggthemes)
library(LSD)

argv<-commandArgs(trailingOnly = TRUE)
data<-read.table(argv[1], header=T)
fig<-argv[2]

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

# Get density of points in 2 dimensions.
# @param x A numeric vector.
# @param y A numeric vector.
# @param n Create a square n by n grid to compute density.
# @return The density within each square.
get_density <- function(x, y, ...) {
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
font_family="sans"
font_color="black"
# text
plot_title="Read count \u2265 2\nHigh-confidence BSJ read count"
legend_title="Legend"
xlab="isoCirc long-read count"
ylab="Short-read count"

scal_f <- function(x){log2(x+1)} # {x}
data$ShortCnt<-scal_f(data$ShortCnt)
data$LongCnt<-scal_f(data$LongCnt)
# get density
data$Density <- get_density(data$LongCnt, data$ShortCnt)
spear_cor = round(cor(data$ShortCnt, data$LongCnt, method = 'spearman'), digits = 2)
spearman_text = expression(paste('Spearman ', rho, ' = 0.5'))
ylim<-12
xlim<-12

ls=c("0.00", "0.05", "0.10", "0.15")
bs=seq(0, 0.15, 0.05)
ggplot(data, aes(x=LongCnt, y=ShortCnt, color=Density)) + 
  geom_point() +
  geom_text(x=1.584962500721156,y=11, label=spearman_text, family=font_family, color=font_color, size=8, hjust=0) +
  scale_colour_gradientn(colours = colorpalette('heat', 5), limits= c(-0.01,0.15), breaks=bs, labels=ls) + 
  ggtitle(plot_title) +
  xlab(xlab) +
  ylab(ylab) +
  scale_y_continuous(breaks=c(0, 1.584962500721156, 3.4594316186372973, 6.658211482751795, 9.967226258835993, 12.288000889707574), labels=c("0", "2", "10", "100", "1,000", "5,000"), limits=c(0,12.288000889707574)) +
  scale_x_continuous(breaks=c(0, 1.584962500721156, 3.4594316186372973, 6.658211482751795, 9.967226258835993, 12.288000889707574), labels=c("0", "2", "10", "100", "1,000", "5,000"), limits=c(0,12.288000889707574)) +
  theme(text=element_text(size=font_size, family=font_family), plot.title = element_text(hjust = 0.5), legend.position= "right", axis.text=element_text(color="black")) +
ggsave(fig, plot=last_plot(), dpi=300, width=9, height=8)