library(ggplot2)

argv<-commandArgs(trailingOnly = TRUE)

data<-read.table(argv[1], header=T)
fig<-argv[2]

# parameters
# color
color_palatte<-"Set1" # http://mkweb.bcgsc.ca/brewer/swatches/brewer-palettes-swatches.pdf
# font size
font_size=14
font_color="black"
font_family="sans"
text_font_size=10
plot_title="Consensus sequence" #Fraction of BSJs"
legend_title="" #Category"
xlab=""
ylab="Length (bp)"

raw_name <- c('Nano-43', 'Nano-44', 'Nano-45', 'Nano-58', 'Nano-59', 
  'Nano-60', 'Adipose', 'Adrenal', 'Blood', 'Brain', 'Heart', 'Kidney', 'Liver', 'Lung', 
  'Prostate', 'SkeletalMuscle', 'SmoothMuscle', 'Testis')
new_name <- c('HEK293 rep1_1', 'HEK293 rep1_2', 'HEK293 rep1_3', 'HEK293 rep2_1', 'HEK293 rep2_2', 
  'HEK293 rep2_3', 'Adipose', 'Adrenal', 'Blood', 'Brain', 'Heart', 'Kidney', 'Liver', 'Lung', 
  'Prostate', 'Skeletal Muscle', 'Smooth Muscle', 'Testis')
data$Name <- factor(data$Name, levels = rev(raw_name), 	labels = rev(new_name), ordered=T)
# data$Type <- factor(data$Type, levels = c('readLen', 'consLen'), label=c('Read', 'Consensus sequence'))

scal_f<-function(x){log2(x)}
# ylims=10000 
# ylims=500
# data$Len<-data[which(Len<=ylims),]
# ylims=6000 #300 
# ys=c(0, 2000, 4000, 6000)
# ylabs = c("60", "500", "4,000", "20,000", "200,000")
# ys = c(5.906890595608519, 8.965784284662087, 11.965784284662087, 14.287712379549449, 17.609640474436812)
# ylims = c(5.906890595608519, 17.609640474436812)

ylabs = c("30", "60",  "100", "300", "1,000", "2,000")
ys = c(4.906890595608519, 5.906890595608519, 6.643856189774724, 8.228818690495881, 9.965784284662087, 10.965784284662087)
ylims = c(4.906890595608519, 10.965784284662087)

ggplot(data, aes(x=Name, y=scal_f(Len), fill=Name)) +
    # geom_bar(stat = "identity", width=0.7, fill="red") +
    geom_violin(show.legend = FALSE) +
    # geom_boxplot(show.legend = FALSE)+
    coord_flip() +  
    scale_y_continuous(breaks=ys, labels=ylabs, limits=ylims) +  
    ggtitle(plot_title) +
    xlab(xlab) + 
    ylab(ylab) +
    labs(fill=legend_title) + 
    theme(text=element_text(size=font_size, family=font_family), 
          plot.title = element_text(hjust = 0.5), axis.text=element_text(color="black"),
          axis.text.x = element_text(color=font_color, size=font_size, family=font_family),
          axis.text.y = element_text(color=font_color, size=font_size, family=font_family)) #, angle=60, hjust=1))
   	# facet_wrap(. ~ Type, scales = "free", nrow=1) +
    # guides(fill = guide_legend(reverse = TRUE))
    #+ scale_fill_brewer(palette = color_palatte)
ggsave(fig, plot=last_plot(), dpi=300, width=8, height=7)
# ggsave(fig, plot=last_plot(), dpi=300, width=12, height=5)