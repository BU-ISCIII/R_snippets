# R script for variant stats in AFumigatus_EMellado Service

# 1) Read files in directories
library(plyr)
library(reshape)
source("functions.R")

## Read variants files
a1163_files <- print(list.files("A1163",pattern="^[0-9a-zA-Z/-]+.vcf.table",full.names=T))
af293_files <- print(list.files("Af293",pattern="^[0-9a-zA-Z/-]+.vcf.table",full.names=T))
pheno <- read.table("pheno.csv",sep="\t",header=T)

variants_A1163 <- read_files_variants(a1163_files,pheno)
variants_Af293 <- read_files_variants(af293_files,pheno)

## Read annotation files
a1163_annot_files <- print(list.files("A1163",pattern="*annot.vcf.table",full.names=T))
af293_annot_files <- print(list.files("Af293",pattern="*annot.vcf.table",full.names=T))
variants_annot_A1163 <- read_files_annot(a1163_annot_files,pheno)
variants_annot_Af293 <- read_files_annot(af293_annot_files,pheno)

### TABLE 1: VARIANTS NUMBER VS DIFFERENT REFS ###
table1 <- NULL

## Number of variants
## A1163
number_variants_a1163 <- ddply(variants_A1163,.(Sample,Phenotype),summarize,Variants_vs_A1163=length(variant))
## Af293                                                                                                                        
number_variants_af293 <- ddply(variants_Af293,.(Sample,Phenotype),summarize,Variants_vs_Af293=length(variant)) 

table1 <- cbind(number_variants_a1163,number_variants_af293[,3])
colnames(table1) <- c("Sample","Phenotype","Variants_vs_A1163","Variants_vs_Af293")
write.table(table1,file="table1_numberVariants_vs_references.txt",row.names=F,sep="\t")

### TABLE 2: Variant Annotation statistics ###
## Table 2.1: A1163
table2_1 <- NULL

variants_annot_stats_a1163 <- ddply(variants_annot_A1163,.(Sample),summarize,Variants_Annotated=length(unique(location)),Number_Genes=length(gene),Number_transcripts=length(feature))

consequence_list <- strsplit(as.character(variants_annot_A1163$consequence),",")
variants_annot_A1163$consequence <- sapply(consequence_list,"[[",1)
variants_annot_conseq_a1163 <- ddply(variants_annot_A1163,.(Sample,consequence),summarize,number=length(location))
variants_annot_conseq_a1163 <- recast(variants_annot_conseq_a1163,Sample~consequence)

table2_1 <- cbind(number_variants_a1163,variants_annot_stats_a1163[,2:ncol(variants_annot_stats_a1163)],variants_annot_conseq_a1163[,2:ncol(variants_annot_conseq_a1163)])
table2_1 <- t(table2_1)
write.table(table2_1,file="table2_1_A1163_consequence_stats.txt",row.names=T,col.names=F,sep="\t")

## Table 2.2: Af293
table2_2 <- NULL

variants_annot_stats_af293 <- ddply(variants_annot_Af293,.(Sample),summarize,Variants_Annotated=length(unique(location)),Number_Genes=length(gene),Number_transcripts=length(feature))

consequence_list <- strsplit(as.character(variants_annot_Af293$consequence),",")
variants_annot_Af293$consequence <- sapply(consequence_list,"[[",1)
variants_annot_conseq_af293 <- ddply(variants_annot_Af293,.(Sample,consequence),summarize,number=length(location))
variants_annot_conseq_af293 <- recast(variants_annot_conseq_af293,Sample~consequence)

table2_2 <- cbind(number_variants_af293,variants_annot_stats_af293[,2:ncol(variants_annot_stats_af293)],variants_annot_conseq_af293[,2:ncol(variants_annot_conseq_af293)])
table2_2 <- t(table2_2)
table2_2[is.na(table2_2)] <- 0
write.table(table2_2,file="table2_2_Af293_consequence_stats.txt",row.names=T,col.names=F,sep="\t")


#### TABLE 3 : Number of variants unique intragroup 
## TABLE 3.1 A1163

number_variants_unique_pop <- ddply(variants_A1163,.(Phenotype),summarize,number=length(unique(variant)))
variants_unique_A1163_wt <- unique(variants_A1163$variant[variants_A1163$Phenotype=="WT"])
variants_unique_A1163_wt <- gsub("_",":",variants_unique_A1163_wt)
variants_unique_A1163_L1 <- unique(variants_A1163$variant[variants_A1163$Phenotype=="L1"])
variants_unique_A1163_L1 <- gsub("_",":",variants_unique_A1163_L1)
variants_unique_A1163_L2 <- unique(variants_A1163$variant[variants_A1163$Phenotype=="L2"])
variants_unique_A1163_L2 <- gsub("_",":",variants_unique_A1163_L2)

variants_unique_annot_A1163 <- variants_annot_A1163[variants_annot_A1163$location %in% variants_unique_A1163_wt & variants_annot_A1163$location %in% variants_unique_A1163_L1 & variants_annot_A1163$location %in% variants_unique_A1163_L2, ]
variants_unique_A1163_conseq <- ddply(variants_unique_annot_A1163,.(Phenotype,consequence),summarize,number=length(location))
variants_unique_A1163_conseq <- recast(variants_unique_A1163_conseq,Phenotype~consequence)

table3_1 <- cbind(number_variants_unique_pop,variants_unique_A1163_conseq[,2:ncol(variants_unique_A1163_conseq)])
table3_1 <- t(table3_1)
write.table(table3_1,file="table3_1_number_variants_A1163_unique_between_populations.txt",row.names=T,col.names=F,sep="\t")

## TABLE 3.2 Af293
number_variants_unique_pop <- ddply(variants_Af293,.(Phenotype),summarize,number=length(unique(variant)))
variants_unique_Af293_wt <- unique(variants_Af293$variant[variants_Af293$Phenotype=="WT"])
variants_unique_Af293_wt <- gsub("_",":",variants_unique_Af293_wt)
variants_unique_Af293_L1 <- unique(variants_Af293$variant[variants_Af293$Phenotype=="L1"])
variants_unique_Af293_L1 <- gsub("_",":",variants_unique_Af293_L1)
variants_unique_Af293_L2 <- unique(variants_Af293$variant[variants_Af293$Phenotype=="L2"])
variants_unique_Af293_L2 <- gsub("_",":",variants_unique_Af293_L2)

variants_unique_annot_Af293 <- variants_annot_Af293[variants_annot_Af293$location %in% variants_unique_Af293_wt & variants_annot_Af293$location %in% variants_unique_Af293_L1 & variants_annot_Af293$location %in% variants_unique_Af293_L2, ]
variants_unique_Af293_conseq <- ddply(variants_unique_annot_Af293,.(Phenotype,consequence),summarize,number=length(location))
variants_unique_Af293_conseq <- recast(variants_unique_Af293_conseq,Phenotype~consequence)

table3_2 <- cbind(number_variants_unique_pop,variants_unique_Af293_conseq[,2:ncol(variants_unique_Af293_conseq)])
table3_2 <- t(table3_2)
write.table(table3_2,file="table3_2_number_variants_Af293_unique_between_populations.txt",row.names=T,col.names=F,sep="\t")
            

### TABLE 4:
#A1163
# Variants unique per Sample. They are not in the rest of Samples.
variants_in_Samples <- ddply(variants_A1163,.(variant),summarize,number_Samples=length(Sample)) 
variants_diff <- variants_in_Samples$variant[variants_in_Samples$number_Samples == 1]
variants_diff_all <- variants_A1163[variants_A1163$variant %in% variants_diff,]
variant_diff_Sample <- ddply(variants_diff_all,.(Sample),summarize,number_variants_diff_Sample=length(variant))      

variants_diff <- gsub("_",":",variants_diff)
variants_diff_annot <- variants_annot_A1163[variants_annot_A1163$location %in% variants_diff,]
variants_diff_annot <- ddply(variants_diff_annot,.(Sample,consequence),summarize,number_variants_diff_Sample=length(location))
variants_diff_annot <- recast(variants_diff_annot,Sample~consequence)

table4_1 <- cbind(variant_diff_Sample,variants_diff_annot[,2:ncol(variants_diff_annot)])
table4_1 <- t(table4_1)
table4_1[is.na(table4_1)] <- 0
write.table(table4_1,file="table4_1_number_variants_A1163_unique_per_Sample.txt",row.names=T,col.names=F,sep="\t")

#Af293
# Variants unique per Sample. They are not in the rest of Samples.                                                              
variants_in_Samples <- ddply(variants_Af293,.(variant),summarize,number_Samples=length(Sample)) 
variants_diff <- variants_in_Samples$variant[variants_in_Samples$number_Samples == 1]
variants_diff_all <- variants_Af293[variants_Af293$variant %in% variants_diff,]
variant_diff_Sample <- ddply(variants_diff_all,.(Sample),summarize,number_variants_diff_Sample=length(variant))      

variants_diff <- gsub("_",":",variants_diff)
variants_diff_annot <- variants_annot_Af293[variants_annot_Af293$location %in% variants_diff,]
variants_diff_annot <- ddply(variants_diff_annot,.(Sample,consequence),summarize,number_variants_diff_Sample=length(location))
variants_diff_annot <- recast(variants_diff_annot,Sample~consequence)

table4_2 <- cbind(variant_diff_Sample,variants_diff_annot[,2:ncol(variants_diff_annot)])
table4_2 <- t(table4_2)
table4_2[is.na(table4_2)] <- 0
write.table(table4_2,file="table4_2_number_variants_Af293_unique_per_Sample.txt",row.names=T,col.names=F,sep="\t")

#######
## TABLE 5
#######
## Variants in L1 but not in L2 and WT
#####
## A1163
####
genotypes_A1163 <- read.table("A1163/snpma.table",header=F)
samples <- read.table("samples_id.txt")
samples <- samples$V1
samples <- sort(samples)
col_names1 <- c("CHROM","POS","REF","ALT","FILTER")
col_names2 <- c("GT","SDP","RD","AD")
col_names2 <- apply(as.matrix(samples),1,function(x) paste(x,col_names2,sep="_"))
col_names2 <- as.vector(col_names2)
col_names <- c(col_names1,col_names2)
colnames(genotypes_A1163) <- col_names
genotypes_A1163$position <- paste(genotypes_A1163$CHROM,genotypes_A1163$POS,sep=":")
genotypes_A1163 <- genotypes_A1163[genotypes_A1163$ALT != ".",]

diff_L1_L2 <- setdiff(variants_unique_A1163_L1,variants_unique_A1163_L2)
diff_L1_L2_WT <- setdiff(diff_L1_L2,variants_unique_A1163_wt)
diff_L1_L2_WT <- gsub(":","_",diff_L1_L2_WT)

variants_ONLY_L1 <- variants_A1163[variants_A1163$variant %in% diff_L1_L2_WT,]

variants_in_samples_L1 <- ddply(variants_ONLY_L1,.(variant),summarize,number_samples=length(Sample))
number_samples_variants_ONLY_L1 <- ddply(variants_in_samples_L1,.(number_samples),summarize,number_variants=length(variant))
write.table(number_samples_variants_ONLY_L1,file="table5_1_distribution_variants_samples_ONLY_L1_A1163.txt",row.names=T,col.names=F,sep="\t")

diff_L1_L2_WT <- gsub("_",":",diff_L1_L2_WT)
variants_annot_A1163_nosample <- variants_annot_A1163[,1:14]
variants_annot_A1163_nosample <- unique(variants_annot_A1163_nosample)
variants_annot_ONLY_L1 <- variants_annot_A1163_nosample[variants_annot_A1163_nosample$location %in% diff_L1_L2_WT,]
variants_annot_ONLY_L1_fil <- variants_annot_ONLY_L1[variants_annot_ONLY_L1$consequence == "stop_gained" | variants_annot_ONLY_L1$consequence == "missense_variant",]
genotypes_A1163_variants_annot_ONLY_L1 <- merge(genotypes_A1163,variants_annot_ONLY_L1_fil,by.x="position",by.y="location",all.y=T)
write.table(genotypes_A1163_variants_annot_ONLY_L1,file="table5_2_only_L1_missense_stop_gained_A1163.txt",sep="\t",row.names=F)

#######
## Af293
########
genotypes_Af293 <- read.table("Af293/snpma.table",header=F)
colnames(genotypes_Af293) <- col_names
genotypes_Af293$position <- paste(genotypes_Af293$CHROM,genotypes_Af293$POS,sep=":")
genotypes_Af293 <- genotypes_Af293[genotypes_Af293$ALT != ".",]


diff_L1_L2 <- setdiff(variants_unique_Af293_L1,variants_unique_Af293_L2)
diff_L1_L2_WT <- setdiff(diff_L1_L2,variants_unique_Af293_wt)
diff_L1_L2_WT <- gsub(":","_",diff_L1_L2_WT)

variants_ONLY_L1 <- variants_Af293[variants_Af293$variant %in% diff_L1_L2_WT,]

variants_in_samples_L1 <- ddply(variants_ONLY_L1,.(variant),summarize,number_samples=length(Sample))
number_samples_variants_ONLY_L1 <- ddply(variants_in_samples_L1,.(number_samples),summarize,number_variants=length(variant))
write.table(number_samples_variants_ONLY_L1,file="table5_3_distribution_variants_samples_ONLY_L1_Af293.txt",row.names=T,col.names=F,sep="\t")

diff_L1_L2_WT <- gsub("_",":",diff_L1_L2_WT)
variants_annot_Af293_nosample <- variants_annot_Af293[,1:14]
variants_annot_Af293_nosample <- unique(variants_annot_Af293_nosample)
variants_annot_ONLY_L1 <- variants_annot_Af293_nosample[variants_annot_Af293_nosample$location %in% diff_L1_L2_WT,]
variants_annot_ONLY_L1_fil <- variants_annot_ONLY_L1[variants_annot_ONLY_L1$consequence == "stop_gained" | variants_annot_ONLY_L1$consequence == "missense_variant",]
genotypes_Af293_variants_annot_ONLY_L1 <- merge(genotypes_Af293,variants_annot_ONLY_L1_fil,by.x="position",by.y="location",all.y=T)
write.table(genotypes_Af293_variants_annot_ONLY_L1,file="table5_4_only_L1_missense_stop_gained_Af293.txt",sep="\t",row.names=F)


## Variants in L2 but not in L1 and WT

####
## A1163
####
diff_L2_L1 <- setdiff(variants_unique_A1163_L2,variants_unique_A1163_L1)
diff_L2_L1_WT <- setdiff(diff_L2_L1,variants_unique_A1163_wt)
diff_L2_L1_WT <- gsub(":","_",diff_L2_L1_WT)

variants_ONLY_L2 <- variants_A1163[variants_A1163$variant %in% diff_L2_L1_WT,]

variants_in_samples_L2 <- ddply(variants_ONLY_L2,.(variant),summarize,number_samples=length(Sample))
number_samples_variants_ONLY_L2 <- ddply(variants_in_samples_L2,.(number_samples),summarize,number_variants=length(variant))
write.table(number_samples_variants_ONLY_L2,file="table5_5_distribution_variants_samples_ONLY_L2_A1163.txt",row.names=T,col.names=F,sep="\t")

diff_L2_L1_WT <- gsub("_",":",diff_L2_L1_WT)
variants_annot_A1163_nosample <- variants_annot_A1163[,1:14]
variants_annot_A1163_nosample <- unique(variants_annot_A1163_nosample)
variants_annot_ONLY_L2 <- variants_annot_A1163_nosample[variants_annot_A1163_nosample$location %in% diff_L2_L1_WT,]
variants_annot_ONLY_L2_fil <- variants_annot_ONLY_L2[variants_annot_ONLY_L2$consequence == "stop_gained" | variants_annot_ONLY_L2$consequence == "missense_variant",]
genotypes_A1163_variants_annot_ONLY_L2 <- merge(genotypes_A1163,variants_annot_ONLY_L2_fil,by.x="position",by.y="location",all.y=T)
write.table(genotypes_A1163_variants_annot_ONLY_L2,file="table5_6_only_L2_missense_stop_gained_A1163.txt",sep="\t",row.names=F)

####
## Af293
####
diff_L2_L1 <- setdiff(variants_unique_Af293_L2,variants_unique_Af293_L1)
diff_L2_L1_WT <- setdiff(diff_L2_L1,variants_unique_Af293_wt)
diff_L2_L1_WT <- gsub(":","_",diff_L2_L1_WT)

variants_ONLY_L2 <- variants_Af293[variants_Af293$variant %in% diff_L2_L1_WT,]

variants_in_samples_L2 <- ddply(variants_ONLY_L2,.(variant),summarize,number_samples=length(Sample))
number_samples_variants_ONLY_L2 <- ddply(variants_in_samples_L2,.(number_samples),summarize,number_variants=length(variant))
write.table(number_samples_variants_ONLY_L2,file="table5_7_distribution_variants_samples_ONLY_L2_Af293.txt",row.names=T,col.names=F,sep="\t")

diff_L2_L1_WT <- gsub("_",":",diff_L2_L1_WT)
variants_annot_Af293_nosample <- variants_annot_Af293[,1:14]
variants_annot_Af293_nosample <- unique(variants_annot_Af293_nosample)
variants_annot_ONLY_L2 <- variants_annot_Af293_nosample[variants_annot_Af293_nosample$location %in% diff_L2_L1_WT,]
variants_annot_ONLY_L2_fil <- variants_annot_ONLY_L2[variants_annot_ONLY_L2$consequence == "stop_gained" | variants_annot_ONLY_L2$consequence == "missense_variant",]
genotypes_Af293_variants_annot_ONLY_L2 <- merge(genotypes_Af293,variants_annot_ONLY_L2_fil,by.x="position",by.y="location",all.y=T)
write.table(genotypes_Af293_variants_annot_ONLY_L2,file="table5_8_only_L2_missense_stop_gained_Af293.txt",sep="\t",row.names=F)


## Variants in WT but not in L1 and L2
####
## A1163
####
diff_WT_L1 <- setdiff(variants_unique_A1163_wt,variants_unique_A1163_L1)
diff_WT_L1_L2 <- setdiff(diff_WT_L1,variants_unique_A1163_L2)
diff_WT_L1_L2 <- gsub(":","_",diff_WT_L1_L2)

variants_ONLY_WT <- variants_A1163[variants_A1163$variant %in% diff_WT_L1_L2,]

variants_in_samples_WT <- ddply(variants_ONLY_WT,.(variant),summarize,number_samples=length(Sample))
number_samples_variants_ONLY_WT <- ddply(variants_in_samples_WT,.(number_samples),summarize,number_variants=length(variant))
write.table(number_samples_variants_ONLY_WT,file="table5_9_distribution_variants_samples_ONLY_WT_A1163.txt",row.names=T,col.names=F,sep="\t")

diff_WT_L1_L2 <- gsub("_",":",diff_WT_L1_L2)
variants_annot_A1163_nosample <- variants_annot_A1163[,1:14]
variants_annot_A1163_nosample <- unique(variants_annot_A1163_nosample)
variants_annot_ONLY_WT <- variants_annot_A1163_nosample[variants_annot_A1163_nosample$location %in% diff_WT_L1_L2,]
variants_annot_ONLY_WT_fil <- variants_annot_ONLY_WT[variants_annot_ONLY_WT$consequence == "stop_gained" | variants_annot_ONLY_WT$consequence == "missense_variant",]
genotypes_A1163_variants_annot_ONLY_WT <- merge(genotypes_A1163,variants_annot_ONLY_WT_fil,by.x="position",by.y="location",all.y=T)
write.table(genotypes_A1163_variants_annot_ONLY_WT,file="table5_10_only_WT_missense_stop_gained_A1163.txt",sep="\t",row.names=F)

####
## Af293
####
diff_WT_L1 <- setdiff(variants_unique_Af293_wt,variants_unique_Af293_L1)
diff_WT_L1_L2 <- setdiff(diff_WT_L1,variants_unique_Af293_L2)
diff_WT_L1_L2 <- gsub(":","_",diff_WT_L1_L2)

variants_ONLY_WT <- variants_Af293[variants_Af293$variant %in% diff_WT_L1_L2,]

variants_in_samples_WT <- ddply(variants_ONLY_WT,.(variant),summarize,number_samples=length(Sample))
number_samples_variants_ONLY_WT <- ddply(variants_in_samples_WT,.(number_samples),summarize,number_variants=length(variant))
write.table(number_samples_variants_ONLY_WT,file="table5_11_distribution_variants_samples_ONLY_WT_Af293.txt",row.names=T,col.names=F,sep="\t")

diff_WT_L1_L2 <- gsub("_",":",diff_WT_L1_L2)
variants_annot_Af293_nosample <- variants_annot_Af293[,1:14]
variants_annot_Af293_nosample <- unique(variants_annot_Af293_nosample)
variants_annot_ONLY_WT <- variants_annot_Af293_nosample[variants_annot_Af293_nosample$location %in% diff_WT_L1_L2,]
variants_annot_ONLY_WT_fil <- variants_annot_ONLY_WT[variants_annot_ONLY_WT$consequence == "stop_gained" | variants_annot_ONLY_WT$consequence == "missense_variant",]
genotypes_Af293_variants_annot_ONLY_WT <- merge(genotypes_Af293,variants_annot_ONLY_WT_fil,by.x="position",by.y="location",all.y=T)
write.table(genotypes_Af293_variants_annot_ONLY_WT,file="table5_12_only_WT_missense_stop_gained_Af293.txt",sep="\t",row.names=F)


#######
## Extra
####

# #A1163
# ## Number of variants shared by the two populations. 26952 variants
# unique_variants_wt <- unique(variants_A1163$variant[variants_A1163$Phenotype=="WT"])
# length(unique(variants_A1163$variant[variants_A1163$variant %in% unique_variants_wt & variants_A1163$Phenotype == "I-group"]))


# ## Af293                                                                                                                                                                                                                                                         
                                                                                                                                 
# ## Number of variants shared by the two populations.
# unique_variants_wt <- unique(variants_Af293$variant[variants_Af293$Phenotype=="WT"])                                           
# length(unique(variants_Af293$variant[variants_Af293$variant %in% unique_variants_wt & variants_Af293$Phenotype == "I-group"])) 
    
# ### TABLE 6:
# ## Number of variants in common among each populations

# #A1163
# # WT. 4399 common variants in the two Samples. 36105 different.
# variants_wt <- variants_A1163[variants_A1163$Phenotype == "WT",]
# variants_in_Samples <- ddply(variants_wt,.(variant),summarize,number_Samples=length(Sample))
# length(variants_in_Samples$variant[variants_in_Samples$number_Samples == 2])
# length(variants_in_Samples$variant[variants_in_Samples$number_Samples == 1])

# #Af293 
# variants_wt <- variants_Af293[variants_Af293$Phenotype == "WT",]                                                               
# variants_in_Samples <- ddply(variants_wt,.(variant),summarize,number_Samples=length(Sample))                                    
# length(variants_in_Samples$variant[variants_in_Samples$number_Samples == 2])                                                    
# length(variants_in_Samples$variant[variants_in_Samples$number_Samples == 1])                                                    


# ## Number of variants in common among each populations                                                                          
# # I-group. 2780 variants in common between the 15 Samples
# variants_igroup <- variants_A1163[variants_A1163$Phenotype == "I-group",]
# variants_in_Samples <- ddply(variants_igroup,.(variant),summarize,number_Samples=length(Sample))
# length(variants_in_Samples$variant[variants_in_Samples$number_Samples == 15])
                                                                                                                                                                                                                                                                                                                                                                                                                               
# # I-group. 
# variants_igroup <- variants_Af293[variants_Af293$Phenotype == "I-group",]                                                      
# variants_in_Samples <- ddply(variants_igroup,.(variant),summarize,number_Samples=length(Sample))                                
# length(variants_in_Samples$variant[variants_in_Samples$number_Samples == 15])                                                   
                                                                                                                                 
