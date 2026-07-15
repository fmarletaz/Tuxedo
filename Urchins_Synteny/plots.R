require(RColorBrewer)
require(pheatmap)
require(tidyverse)
require(wesanderson)
require(stringr)
require(ggalluvial)

source("~/Desktop/Urchins_Synteny//synteny-functions.R")
setwd("~/Desktop/Urchins_Synteny/local/")
#source("../synteny-functions.R")

colors=data.frame(row.names = c('A1','A2','B1','B2+C2','B3','C1','D','E','F','G','H','I','J1','J2','K','L','M','N','O1','O2','P', 'Q', 'R'),
                  col=c("A6CEE3","1F78B4","1B9E77","33A02C","B2DF8A","FB9A99","E31A1C","FFFF99","B15928","D95F02","7570B3","E7298A","66A61E","6A3D9A","CAB2D6","E6AB02","A6761D","8DD3C7","666666","FF7F00","FDBF6F","BEBADA","ffd92f"))

#by position by chromosomes
bmglob_rmglob_synt <- read.table("mglob_output_syn_with_clg-mglob_red_synt.txt", sep='\t', h=T)
bmglob_rmglob_synt.exp <- testEnrich(bmglob_rmglob_synt)
bmglob_rmglob_synt.exp %>% mutate(clg=na_if(clg,"")) -> bmglob_rmglob_synt.exp
bmglob_rmglob_synt.exp$clg <- as.factor(bmglob_rmglob_synt.exp$clg)
bmglob_rmglob_synt.exp$s1chr <- as.factor(bmglob_rmglob_synt.exp$s1chr)
bmglob_rmglob_synt.exp$s2chr <- as.factor(bmglob_rmglob_synt.exp$s2chr)

color_mapping <- setNames(paste0("#", colors$col), rownames(colors))
color_mapping["NA"] <- "#CCCCCC"
  
# Create directory if it doesn't exist
dir.create("bluevd_red/Chromosome_mapping", recursive = TRUE, showWarnings = FALSE)

# Count syntenic blocks for each chromosome pair
chr_counts <- bmglob_rmglob_synt.exp %>%
  group_by(s1chr, s2chr) %>%
  summarise(count = n(), .groups = 'drop') %>%
  arrange(desc(count))

# Find best match for each s1chr (most syntenic blocks)
best_pairs <- chr_counts %>%
  group_by(s1chr) %>%
  slice_max(count, n = 1) %>%
  ungroup()

# Create individual plots for each best chromosome pair
for(i in 1:nrow(best_pairs)) {
  chr1 <- best_pairs$s1chr[i]
  chr2 <- best_pairs$s2chr[i]
  block_count <- best_pairs$count[i]
  
  # Subset data for this chromosome pair
  chr_data <- bmglob_rmglob_synt.exp %>%
    filter(s1chr == chr1 & s2chr == chr2) %>%
    arrange(s1gp, s2gp)
  
  if(nrow(chr_data) >= 1) {
    # Create filename
    filename <- paste0("bluevd_red/Chromosome_mapping/", chr1, "_vs_", chr2, ".png")
    
    # Open PNG device with larger dimensions
    png(filename, width = 1200, height = 800, res = 150)
    
    # Set margins - less space on right, more for plot
    par(mar=c(5, 4, 4, 6) + 0.1)
    
    # Get colors for each CLG
    point_colors <- color_mapping[as.character(chr_data$clg)]
    
    # Create dotplot with CLG colors
    plot(chr_data$s1gp, chr_data$s2gp, 
         xlab = paste(chr1, "position (bp)"), 
         ylab = paste(chr2, "position (bp)"),
         main = paste(chr1, "vs", chr2, "(", block_count, "syntenic blocks)"),
         pch = 19, cex = 1.2, col = point_colors)
    
    # Add smaller legend positioned outside
    unique_clgs <- unique(chr_data$clg[!is.na(chr_data$clg)])
    if(length(unique_clgs) > 0) {
      legend("topright", 
             inset = c(-0.15, 0),  # Position outside plot area
             legend = unique_clgs, 
             fill = color_mapping[as.character(unique_clgs)],
             title = "CLG",
             cex = 0.5,  # Much smaller legend text
             xpd = TRUE)
    }
    
    # Close device
    dev.off()
    
    print(paste("Saved:", filename))
  }
}

print(paste("Total chromosome pairs plotted:", nrow(best_pairs)))
#by index 
bmglob_rmglob_synt=read.table("mglob_output_syn_with_clg-mglob_red_synt.txt",sep='\t',h=T)
bmglob_rmglob_synt.exp<-testEnrich(bmglob_rmglob_synt)
bmglob_rmglob_synt.exp %>% mutate(clg=na_if(clg,"")) -> bmglob_rmglob_synt.exp
bmglob_rmglob_synt.exp$clg<-as.factor(bmglob_rmglob_synt.exp$clg)
par(mar=c(8, 8, 8, 8) + 0.1) 
plotSyntCLGrs(ordSynt(bmglob_rmglob_synt.exp),'M. globulus (red)', 'M. globulus (blue)')
write.table(bmglob_rmglob_synt.exp,file="bmglob_rmglob_syntr.txt",sep='\t',quote=F,row.names=F)

#setwd("~/Dropbox/Genomes/Synteny/synt-forge/Urchins/")
hleu_mglob_synt=read.table("Holothuria_holleu_clean-mglob_output_syn_synt.txt",sep='\t',h=T)
hleu_mglob_synt.exp<-testEnrich(hleu_mglob_synt)
hleu_mglob_synt.exp %>% mutate(clg=na_if(clg,"")) -> hleu_mglob_synt.exp
hleu_mglob_synt.exp$clg<-as.factor(hleu_mglob_synt.exp$clg)
par(mar=c(8, 8, 8, 8) + 0.1) 
plotSyntCLGrs(ordSynt(hleu_mglob_synt.exp),'M. globulus', 'H. leucospilota')
write.table(hleu_mglob_synt.exp,file="hleu_mglob_syntr.txt",sep='\t',quote=F,row.names=F)

mglob_pliv_synt=read.table("mglob_output_syn_with_clg-pliv_output_synt.txt",sep='\t',h=T)
mglob_pliv_synt.exp<-testEnrich(mglob_pliv_synt)
mglob_pliv_synt.exp %>% mutate(clg=na_if(clg,"")) -> mglob_pliv_synt.exp
mglob_pliv_synt.exp$clg<-as.factor(mglob_pliv_synt.exp$clg)
par(mar=c(8, 8, 8, 8) + 0.1) 
plotSyntCLGrs(ordSynt(mglob_pliv_synt.exp),'P. lividus','M. globulus')
write.table(mglob_pliv_synt.exp,file="mglob_pliv_syntr.txt",sep='\t',quote=F,row.names=F)

pliv_spurp_synt=read.table("pliv_output_with_clg-spur_synt_cl_synt.txt",sep='\t',h=T)
pliv_spurp_synt.exp<-testEnrich(pliv_spurp_synt)
pliv_spurp_synt.exp %>% mutate(clg=na_if(clg,"")) -> pliv_spurp_synt.exp
pliv_spurp_synt.exp$clg<-as.factor(pliv_spurp_synt.exp$clg)
par(mar=c(8, 8, 8, 8) + 0.1) 
plotSyntCLGrs(ordSynt(pliv_spurp_synt.exp),'S. Purpuratus','P. lividus')
write.table(pliv_spurp_synt.exp,file="pliv_spurp_syntr.txt",sep='\t',quote=F,row.names=F)

spurp_lpict_synt=read.table("spur-lpic_cleaned_synt.txt",sep='\t',h=T)
spurp_lpict_synt.exp<-testEnrich(spurp_lpict_synt)
spurp_lpict_synt.exp %>% mutate(clg=na_if(clg,"")) -> spurp_lpict_synt.exp
spurp_lpict_synt.exp$clg<-as.factor(spurp_lpict_synt.exp$clg)
par(mar=c(8, 8, 8, 8) + 0.1) 
plotSyntCLGrs(ordSynt(spurp_lpict_synt.exp),'L. pictus','S. Purpuratus')
write.table(spurp_lpict_synt.exp,file="spurp_lpict_syntr.txt",sep='\t',quote=F,row.names=F)

lpict_lvar_synt=read.table("lpic_cleaned_with_clg-lvar_output_synt.txt",sep='\t',h=T)
lpict_lvar_synt.exp<-testEnrich(lpict_lvar_synt)
lpict_lvar_synt.exp %>% mutate(clg=na_if(clg,"")) -> lpict_lvar_synt.exp
lpict_lvar_synt.exp$clg<-as.factor(lpict_lvar_synt.exp$clg)
par(mar=c(8, 8, 8, 8) + 0.1) 
plotSyntCLGrs(ordSynt(lpict_lvar_synt.exp),'L. variegatus','L. pictus')
write.table(lpict_lvar_synt.exp,file="lpict_lvar_syntr.txt",sep='\t',quote=F,row.names=F)

lvar_lpict_synt=read.table("lvar_output_with_clg-lpic_cleaned_synt.txt",sep='\t',h=T)
lvar_lpict_synt.exp<-testEnrich(lvar_lpict_synt)
lvar_lpict_synt.exp %>% mutate(clg=na_if(clg,"")) -> lvar_lpict_synt.exp
lvar_lpict_synt.exp$clg<-as.factor(lvar_lpict_synt.exp$clg)
par(mar=c(8, 8, 8, 8) + 0.1) 
plotSyntCLGrs(ordSynt(lvar_lpict_synt.exp),'L. pictus','L. variegatus')
write.table(lvar_lpict_synt.exp,file="lvar_lpict_syntr.txt",sep='\t',quote=F,row.names=F)

spurp_lvar_synt=read.table("spur-lvar_output_synt.txt",sep='\t',h=T)
spurp_lvar_synt.exp<-testEnrich(spurp_lvar_synt)
spurp_lvar_synt.exp %>% mutate(clg=na_if(clg,"")) -> spurp_lvar_synt.exp
spurp_lvar_synt.exp$clg<-as.factor(spurp_lvar_synt.exp$clg)
par(mar=c(8, 8, 8, 8) + 0.1) 
plotSyntCLGrs(ordSynt(spurp_lvar_synt.exp),'L. variegatus','S. purpuratus')
write.table(spurp_lvar_synt.exp,file="spurp_lvar_syntr.txt",sep='\t',quote=F,row.names=F)

