
## R script for tree annotation using ggtree
library(reshape)  
library(magrittr) 
library(ape)
library(phangorn)
library(ggtree)
library(ggrepel)

source("./facet_plot.R")

brewer_qualitative <- c("#0000ff","#ff0000","#483215","#008900","#7244c4","#e65a11","#000000","#e6e528","#ff00ee","#6e0000","#00c7dd","#d4b455","#8f008d","#736b00","#7d8cbf","#c2a91b","#374a12")


## Files
raxml_file <- "./RAxML_bipartitions.RAXML_TREE_ANNOT"
sample_info <- "choleraesuis_format.txt"

## Read tree and sample info
raxml_tree <- read.tree(raxml_file)
info <- read.table(sample_info,sep="\t",header=T)

## Midroot tree                    
raxml_tree <- midpoint(raxml_tree) 

## Print first tree with node labels for improving visualization.              
p0 <- ggtree(raxml_tree) + geom_text2(aes(subset=!isTip,label=node),hjust=-.3) + geom_tiplab()  
pdf(file="tree_nodes.pdf",width=40,height=50)                                  
	print(p0)                                                                  
dev.off()                                                                      

## Remove nodes spoiling the image because of big branch  
clade_remote <- c("SRR1122712","SRR5081321","SRR1969188") 
raxml_tree_trim <- drop.tip(raxml_tree,clade_remote)      

## Group clades by ST                                                     
groupInfo <- split(as.character(info$sample),as.character(info$ST))   
raxml_tree_trim_groupInfo <- groupOTU(raxml_tree_trim,groupInfo)          

## Bootstrap manipulation. Only show bootstrap more than 80.
p <- ggtree(raxml_tree_trim_groupInfo)
bootstrap_df <- p$data
bootstrap_df <- bootstrap_df[!bootstrap_df$isTip,]
bootstrap_df$label <- as.numeric(bootstrap_df$label)
bootstrap_df <- bootstrap_df[bootstrap_df$label > 80,]
bootstrap_df$.panel <- factor("Tree")

## Cladogram with branch annotation
p1 <- ggtree(raxml_tree_trim_groupInfo,branch.length="none") %<+% info + 
  #geom_text(aes(x=branch,label=branch.length,vjust=-.5)) +
  geom_label_repel(data=bootstrap_df,aes(label=label)) +
  geom_tiplab(size=6,align=TRUE) + 
  #geom_treescale(width=15000,color="red",offset=1.5,linesize = 1.2) +
  geom_tippoint(aes(color=ST),size=7,alpha=0.25) +
  geom_text(aes(label=Source_Type),hjust=1.4,vjust=-0.4,size=3) +
  scale_color_manual("ST",values=brewer_qualitative) + 
  geom_text(aes(label=Country),hjust=1.4,vjust=1.4,size=3) +
  theme(legend.position="right")

## Phylogram with branch annotation 
p2 <- ggtree(raxml_tree_trim_groupInfo) %<+% info +              
   #geom_text(aes(x=branch,label=branch.length,vjust=-.5)) +            
   geom_text2(aes(x=branch,subset=!isTip,label=label,size=4,vjust=-.5)) +         
   geom_tiplab(size=6,align=TRUE) +                                     
   #geom_treescale(width=15000,color="red",offset=1.5,linesize = 1.2) + 
   geom_tippoint(aes(color=ST),size=7,alpha=0.25) +                     
   geom_text(aes(label=Source_Type),hjust=1.4,vjust=-0.4,size=3) +      
   theme(legend.position="right")                                       

## Phylogram  without annotation for side table.   
p3 <- ggtree(raxml_tree_trim_groupInfo,aes(color=group),size=1.5) %<+% info +                                     
	#geom_text2(aes(label=label, subset = !is.na(as.numeric(label)) & as.numeric(label) > 80),size=6,vjust=-.5,hjust=-.4) +	
	#geom_text(data=bootstrap_df,aes(label=label),size=4,vjust=-.5) +             
	geom_label_repel(data=bootstrap_df,aes(label=label),size=6) + 
    geom_tiplab(size=6,align=TRUE) +                                        
    scale_color_brewer("ST",palette="Spectral") +
    scale_size(guide=FALSE) +
    #xlim(0,0.036) +
    #theme_tree2() +
    theme(legend.position="bottom",legend.text = element_text(size=25),legend.title=element_text(size=25),legend.box="horizontal")                                          

## Table preparation for side table
d3 <- data.frame(sample=info$sample,country=info$Country,date=info$Collection_Year,source=info$Source_Type)
# Melt table
d3 <- melt(d3,id=c("sample"))

# Create "fake" continous scale for dimension the table. (It fails with discrete values)
# d3$variable <- gsub("st",0.001,d3$variable)  
# d3$variable <- gsub("country",0.005,d3$variable)
# d3$variable <- gsub("date",0.01,d3$variable) 
# d3$variable <- gsub("source",0.015,d3$variable)   
# d3$variable <- as.numeric(d3$variable)

d3$variable <- gsub("country",0,d3$variable)   
d3$variable <- gsub("date",0.005,d3$variable)       
d3$variable <- gsub("source",0.01,d3$variable)    
d3$variable <- as.numeric(d3$variable)             

## Create facet grid table next to the tree
# Custom face_plot, space="free_x" added.
# inward ajust for visualization of the labels
p4 <- facet_plot(p3,panel="data",data=d3,grid_space="free_x",geom=geom_text,aes(x=variable,label=value),size=8,vjust = "inward", hjust = "inward")

## Tree painting into files
pdf(file="tree_cladogram_annot_tipsalign.pdf",width=40,height=50)
  print(p1)
dev.off()

pdf(file="tree_phylogram_tipsalign_annot.pdf",width=40,height=70)
  print(p2)
dev.off()

pdf(file="tree_cladogram_annot_test_notable.pdf",width=40,height=70) 
   print(p3)                                                       
dev.off()                                                         

pdf(file="tree_cladogram_annot_test_table.pdf",width=40,height=70)  
    print(p4)                                                          
dev.off()                                                             
