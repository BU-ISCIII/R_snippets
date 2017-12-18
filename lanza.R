source("script_comp.R")

all_variants <- data(dbname="/home/saram/Documentos/Desarrollo/virtual_tumors/20130618_virtualTumor_freq_10_EXOME/20130716_0921_Report_Variants/20130716_0921_mutations.sqlite",
					 gatk_c="/home/saram/Documentos/Desarrollo/virtual_tumors/gatk_samtools_counts/gatk_virtualFreq_10_2/all_variants_S_gatk_virtualFreq",
					 samtools="/home/saram/Documentos/Desarrollo/virtual_tumors/gatk_samtools_counts/samtools_virtualFreq_10_2/all_variants_samtools_virtualFreq"
					)
cosmic <- cosmic_parse("/home/saram/Documentos/Desarrollo/virtual_tumors/gatk_samtools_counts/cosmic_10")

venn(all_variants)

sumDep <- summary_depth(all_variants)
sumFreq <- summary_freq(all_variants)

virtual_freq(all_variants,cosmic) 
fp_count_f <- fp_count_freq(all_variants,cosmic)
write.csv(fp_count_f,file="fp_count.csv")


sumDepC <- do.call(rbind,sumDep[[2]])
sumDepT <- do.call(rbind,sumDep[[1]])

write.csv(sumDepC,file="sumDepC.csv")
write.csv(sumDepT,file="sumDepT.csv")

sumFreqC <- do.call(rbind,sumFreq[[2]])
sumFreqT <- do.call(rbind,sumFreq[[1]])

write.csv(sumFreqC,file="sumFreqC.csv")
write.csv(sumFreqT,file="sumFreqT.csv")

got(all_variants,cosmic)

