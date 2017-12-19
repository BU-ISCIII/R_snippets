 ## Functions
 
 read_files_variants <- function(files,pheno){                
 	df_bind <- NULL                                          
 	for (f in files){                             
 		df <- read.table(f,sep="\t")                    
 		colnames(df) <- c("chr","pos","ref","alt")      
 		df$variant <- paste(df$chr,df$pos,sep="_")
 		sample <- gsub("^([0-9a-zA-Z]+)/([0-9a-zA-Z]+).*$","\\2",f) 
        df$Sample <- sample 
 		df$reference <- gsub("^([0-9a-zA-Z]+)/([0-9a-zA-Z]+).*$","\\1",f)
 		df$Phenotype <- pheno$resistance[pheno$sample == sample]
 		df_bind <- rbind(df_bind,df)                    
 	}                                                   
     return(df_bind)                                     
 }                                                       

read_files_annot <- function(files,pheno){                
 	df_bind <- NULL                                          
 	for (f in files){                             
 		df <- read.table(f,sep="\t")                    
 		colnames(df) <- c("variation","location","allele","gene","feature","feature_type","consequence","cDNA_position","CDS_position","protein_position","amino_acids","codons","existing_variation","extra")      
 		sample <- gsub("^([0-9a-zA-Z]+)/([0-9a-zA-Z]+).*$","\\2",f) 
        df$Sample <- sample 
 		df$reference <- gsub("^([0-9a-zA-Z]+)/([0-9a-zA-Z]+).*$","\\1",f)
 		df$Phenotype <- pheno$resistance[pheno$sample == sample]
 		df_bind <- rbind(df_bind,df)                    
 	}                                                   
     return(df_bind)                                     
 }
