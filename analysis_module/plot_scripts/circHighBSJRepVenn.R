
library(VennDiagram)
library(RColorBrewer)

argv<-commandArgs(trailingOnly = TRUE)
data<-read.table(argv[1], header=T)
fig<-argv[2]

# parameters
# color
color_palatte<-"Set1" # http://mkweb.bcgsc.ca/brewer/swatches/brewer-palettes-swatches.pdf
fill_color_template=brewer.pal(n = 8, name = color_palatte)
n = 2
fill_color = fill_color_template[1:n]
# font size
font_size=1
font_face="bold"
font_family="sans"
plot_title="" #High-confidence BSJs (read count >=2) \nof 3 replicates"
title_font_size=1.5
# legend_title="Legend"

data$Type <- factor(data$Type, levels = c('rep1', 'rep2'), labels= c('Rep#1', 'Rep#2'), ordered=TRUE)

t1="Rep#1"
t2="Rep#2"
l1=as.character(data[data$Type==t1,]$Value)
l2=as.character(data[data$Type==t2,]$Value)

p <- venn.diagram(x = list(Replicate_1=l1, Replicate_2=l2),
            category = c('Rep#1', 'Rep#2'),
            main=plot_title,
            main.cex=title_font_size, main.fontfamily=font_family, main.fontface = font_face,
            cex = font_size, fontfamily = font_family, fontface = font_face,
            cat.cex = font_size, cat.fontfamily = font_family, cat.fontface = font_face,
            filename=NULL, resolution=1000, height=5000, width=4000,
            fill = fill_color,
            col = "transparent", alpha = 0.50, #cat.col = fill_color,
            ext.text = TRUE)# , cat.cex = 1, cat.pos = c(210,180, 330), cat.dist = c(0.3, 0.4,0.15), rotation.degree = 0, margin = 0.2)
tot_n <- 2^n - 1
idx <- sapply(p, function(i) grepl("text", i$name))
for(i in 1:tot_n){
      p[idx][[i]]$label <- format(as.numeric(p[idx][[i]]$label), big.mark=",", scientific=F)
      # p[idx][[i]]$label <- paste(format(as.numeric(p[idx][[i]]$label), big.mark=",", scientific=F), as.numeric(p[idx][[i]]$label)/100, "%")
}
pdf(fig)
grid.newpage()
grid.draw(p)
dev.off()
