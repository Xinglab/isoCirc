library(ggplot2)
library(ggrepel)


argv<-commandArgs(trailingOnly = TRUE)

data<-read.table(argv[1], header=T)
fig<-argv[2]

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
# font size
# font_size=25
font_size=12
font_family="sans"
font_color="black"
text_font_size=4
plot_title="Fraction of Alu pairs"
legend_title=""
xlab=""
ylab="Fraction"

data <- data[which((data$MinCount == '2' | data$MinCount == '3') & (data$FlankLen == '1000' | data$FlankLen == '2000')),]

# data$Cate <- factor(data$Cate, levels = c(1,2), labels= c('Detected in 1 replicate', 'Detected in 2 replicates'), ordered=T)
data$Cate <- factor(data$Cate, levels = c('Convergent', 'Divergent', 'Both', 'Other'), labels= c('Convergent', 'Divergent', 'Both', 'None'),ordered=T)
data$RepAnno <- factor(data$RepAnno, levels = c('1Novel', '2Novel', '1Known', '2Known', 'Negative'), ordered=T)
# data$RepAnno <- factor(data$RepAnno, levels = c('1Novel', '2Novel', '1Known', '2Known'), labels = c('Novel, Detected in 1 replicate', 'Novel, Detected in 2 replicates', 'Known, Detected in 1 replicate', 'Known, Detected in 2 replicates'), ordered=T)
data$MinCount <- factor(data$MinCount, levels = c('1', '2', '3'), labels=c(expression('Read count \u2265 1'), 'Read count \u2265 2', 'Read count \u2265 3'), ordered=T)
# data$MinCount <- factor(data$MinCount, levels = c('1', '2', '3'), labels=c('Read count ≥ 1', 'Read count ≥ 2', 'Read count ≥ 3'), ordered=T)
data$FlankLen <- factor(data$FlankLen, levels = c('500', '1000', '2000'), labels=c('Flanking window size: 500 bp', 'Flanking window size: 1,000 bp', 'Flanking window size: 2,000 bp'), ordered=T)

ggplot(data, aes(x=RepAnno, y=Frac, fill=Cate)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label=format(as.numeric(Count), big.mark=",", scientific=F)), position = position_stack(vjust=0.5), size = text_font_size, family=font_family) + 
    # geom_text_repel(aes(label = format(as.numeric(Count), big.mark=",", scientific=F)), position=position_stack(vjust=0.5), size = text_font_size, family=font_family, direction='y')+
    ggtitle(plot_title) +
    xlab(xlab) + 
    ylab(ylab) +
    labs(fill=legend_title) +
    # facet_wrap(. ~ MinReadCount, scales = "free", nrow=1) +
    facet_grid(FlankLen ~ MinCount) + 
    theme(text=element_text(size=font_size, family=font_family), 
          plot.title = element_text(hjust = 0.5), axis.text=element_text(color=font_color),
          legend.position="bottom"
          # axis.text.x = element_text(size=font_size, color=font_color, family=font_family, angle=0, hjust=0)
          ) + # put title in the center
    guides(fill=guide_legend(ncol=4)) + 
    scale_fill_brewer(palette = color_palatte)
ggsave(fig, plot=last_plot(), dpi=300, width=7, height=7)
# ggsave(fig, plot=last_plot(), dpi=300, width=16, height=16)