library(ggplot2)

argv<-commandArgs(trailingOnly = TRUE)

data<-read.table(argv[1], header=T)
fig<-argv[2]

# parameters
# color
color_palatte<-"Set1" # http://mkweb.bcgsc.ca/brewer/swatches/brewer-palettes-swatches.pdf
# font size
font_size=20
font_family="sans"
font_color="black"
text_font_size=8
plot_title="Fraction of BSJs"
legend_title="Category"
xlab="Number of replicates"
ylab="Fraction"

# data$Cate <- factor(data$Cate, levels = c(1,2), labels= c('Detected in 1 replicate', 'Detected in 2 replicates'), ordered=T)
data$Cate <- factor(data$Cate, levels = c(1,2), ordered=T)
data$Anno <- factor(data$Anno, levels = c('circBase', 'MiOncoCirc', 'Both', 'Novel'), ordered=T)

ggplot(data, aes(x=Cate, y=Fraction, fill=Anno)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label=format(as.numeric(Cnt), big.mark=",", scientific=F)), position = position_stack(vjust=0.5), size = text_font_size, family=font_family) + 
    ggtitle(plot_title) +
    xlab(xlab) + 
    ylab(ylab) +
    labs(fill=legend_title) + 
    theme(text=element_text(size=font_size, family=font_family), 
          plot.title = element_text(hjust = 0.5), axis.text=element_text(color=font_color)
          # axis.text.x = element_text(size=font_size, color=font_color, family=font_family, angle=0, hjust=0)
          ) + # put title in the center
    scale_fill_brewer(palette = color_palatte)
ggsave(fig, plot=last_plot(), dpi=300, width=6, height=7)