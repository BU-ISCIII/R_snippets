library(VennDiagram)
library(RSQLite)
library(ggplot2)
library(plyr)
library(scales)
library(gridExtra)
library(reshape2)

data <- function(dbname,gatk_c,samtools){ #,tgatkcontrol,tsamtools
	# Conexión con la base de datos.
	drv <- dbDriver("SQLite")
	con <- dbConnect(drv,dbname=dbname)

	# Obtener tablas para comparar.

	variantJoint <- dbGetQuery(con,"select var_sample_id,chr,pos,ref,alt,gt_control,gq_control,gt_tumor,gq_tumor from variant")
	variantStrelka <- dbGetQuery(con,"select var_sample_id,chr,pos,ref,alt,gt_control,gq_control,gt_tumor,gq_tumor from variantStrelka")
	variantGatk <- dbGetQuery(con,"select var_sample_id,chr,pos,ref,alt,gt_control,gq_control,gt_tumor,gq_tumor from variantGatk")
	gatkControl <- read.table(gatk_c,sep="\t",header=T)
	variantSamtools <- read.table(samtools,sep="\t",header=T)

	## Ñapa monumental para poder hacer el gráfico boxplot, no usar para estudio por separado de muestras ni nada parecido, las filas no corresponden
	## con el valor de freqC.
	# newrow <- rep(NA, 7)
	# recm <- matrix(rep(newrow,9),nr=9,byrow=TRUE)
	# gatkControlm <- as.matrix(gatkControl)
	# gatkControl <- data.frame(rbind(gatkControlm,recm))

	## Añadimos una columna con los campos que definen a la variante. Teoricamente debería ser con ref y alt también,
	# pero como samtools para ALT usa una nomenclatura diferente, lo dejamos en chr-pos únicamente.
	variantJoint$merged <- do.call(paste, c(variantJoint[c("chr", "pos")], sep = "_"))
	variantStrelka$merged <- do.call(paste, c(variantStrelka[c("chr","pos")], sep = "_"))
	variantGatk$merged <- do.call(paste, c(variantGatk[c("chr","pos")], sep = "_"))
	variantSamtools$merged <- do.call(paste, c(variantSamtools[c("chr","pos")], sep = "_"))

	##
	# Código para calcular meter columna de cobertura..
	##
	variantJoint$DepthControl <- as.numeric(variantJoint$gt_control) + as.numeric(variantJoint$gq_control)
	variantJoint$DepthTumor <- as.numeric(variantJoint$gt_tumor) + as.numeric(variantJoint$gq_tumor)

	variantStrelka$DepthControl <- as.numeric(variantStrelka$gt_control) + as.numeric(variantStrelka$gq_control)
	variantStrelka$DepthTumor <- as.numeric(variantStrelka$gt_tumor) + as.numeric(variantStrelka$gq_tumor)

	# variantGatk$DepthControl <- NA
	gatkControl$DepthControl <- as.numeric(gatkControl$gt_control) + as.numeric(gatkControl$gq_control)
	variantGatk$DepthControl <- gatkControl$DepthControl
	variantGatk$DepthTumor <- as.numeric(variantGatk$gt_tumor) + as.numeric(variantGatk$gq_tumor)

	variantSamtools$DepthControl <- as.numeric(variantSamtools$gt_control) + as.numeric(variantSamtools$gq_control)
	variantSamtools$DepthTumor <- as.numeric(variantSamtools$gt_tumor) + as.numeric(variantSamtools$gq_tumor)

	##
	# Código para añadir columna de vaf.
	##
	variantJoint$freqT <- as.numeric(variantJoint$gq_tumor)/(as.numeric(variantJoint$gq_tumor)+as.numeric(variantJoint$gt_tumor))
	variantJoint$freqT[variantJoint$freqT == "Inf"] <- NA
	variantJoint$freqC <- as.numeric(variantJoint$gq_control)/(as.numeric(variantJoint$gq_control)+as.numeric(variantJoint$gt_control))
 	variantJoint$freqC[variantJoint$freqC == "Inf"] <- NA

	variantStrelka$freqT <- as.numeric(variantStrelka$gq_tumor)/(as.numeric(variantStrelka$gq_tumor)+as.numeric(variantStrelka$gt_tumor))
	variantStrelka$freqT[variantStrelka$freqT == "Inf"] <- NA
	variantStrelka$freqC <- as.numeric(variantStrelka$gq_control)/(as.numeric(variantStrelka$gq_control)+as.numeric(variantStrelka$gt_control))
	variantStrelka$freqC[variantStrelka$freqC == "Inf"] <- NA
	
	variantGatk$freqT <- as.numeric(variantGatk$gq_tumor)/(as.numeric(variantGatk$gq_tumor)+as.numeric(variantGatk$gt_tumor))
	variantGatk$freqT[variantGatk$freT == "Inf"] <- NA
	gatkControl$freqC <- as.numeric(gatkControl$gq_control)/(as.numeric(gatkControl$gq_control)+as.numeric(gatkControl$gt_control))
	gatkControl$freqC[gatkControl$freqC == "Inf"] <- NA
	# variantGatk$freqC <- NA
	## OJO ÑAPA para el boxplot únicamente.
	variantGatk$freqC <- gatkControl$freqC
	################################
 	
 	variantSamtools$freqT <- as.numeric(variantSamtools$gq_tumor)/(as.numeric(variantSamtools$gq_tumor)+as.numeric(variantSamtools$gt_tumor))
	variantSamtools$freqT[variantSamtools$freT == "Inf"] <- NA
	variantSamtools$freqC <- as.numeric(variantSamtools$gq_control)/(as.numeric(variantSamtools$gq_control)+as.numeric(variantSamtools$gt_control))
	variantSamtools$freqC[variantSamtools$freqC == "Inf"] <- NA

	## Añado una columna que identifique a cada software
	variantJoint$software <- "JointSNVmix"
	variantStrelka$software <- "Strelka"
	variantGatk$software <- "UnifiedGenotyper"
	variantSamtools$software <- "Samtools"

	## Información de distribución de frequencias de snps de hapmap.
	snps124 <- read.table("/home/saram/Documentos/Desarrollo/LowFrequencyValidation/snps124_T_bwa.txt",sep="\t",header=T)
	snps124$merged <- do.call(paste, c(snps124[c("chr", "pos", "ref", "alt")], sep = "_"))
	snps124$freqT <- (snps124$gq_tumor/(snps124$gq_tumor+snps124$gt_tumor))
	snps124$freqT <- (snps124$gq_tumor/(snps124$gq_tumor+snps124$gt_tumor))
	snps124$freqC <- (snps124$gq_control/(snps124$gq_control+snps124$gt_control))
	snps124$freqC <- (snps124$gq_control/(snps124$gq_control+snps124$gt_control))
	snps124$DepthControl <- as.numeric(snps124$gt_control) + as.numeric(snps124$gq_control)
	snps124$DepthTumor <- as.numeric(snps124$gt_tumor) + as.numeric(snps124$gq_tumor)
	snps124$software <- "Hapmap SNP"
	
	## Juntamos todo en una única tabla.
	all_variants <- rbind(variantJoint,variantStrelka,snps124, variantSamtools, variantGatk) 
	
	p <- lapply(all_variants$software,function(x){
		if(x == "JointSNVmix" || x == "Strelka"){
		 	y <- "joint"
	   	}else{
	   		y <- "independent"
	   	}
	  })
	
	all_variants$type <- unlist(p)

	all_variants$het[all_variants$freqT > 0.3] <- "> 0.3"
	all_variants$het[all_variants$freqT <= 0.3] <- "<= 0.3"

	return(all_variants)

}

cosmic_parse <- function(cosmic_dir){
	cosmic_files <- list.files(path=cosmic_dir,full.names=TRUE,pattern=".ok-T.txt")
	cosmic <- data.frame(chr =character(), start =integer(), end =integer(), mut=character(),strand=character(),ref=character(),alt=character(),freq=integer())

	## Creo una tabla que una todos los ficheros .bed, añadiendo una columna a cada una que se corresponda con la frequencia
	# de las mutaciones que contiene.
	for (i in cosmic_files){
		## Separo la ruta de los ficheros para separar vaf y meterlo como columna.
		tmp <- strsplit(i,"[._]",perl=TRUE)
		## Leo el fichero y le cambio los nombres a las columnas.
		cosmic_tmp <- read.table(i,sep="\t")
		names(cosmic_tmp) <- paste(c("chr","pos","ref","alt","gt_tumor","gq_tumor","pileup"))
		cosmic_tmp$pileup <- NULL
		## Meto dos columnas nuevas haciendo un split de la columna mut, para tener ref y alt por separado.
		# cosmic_tmp <- with(cosmic_tmp, cbind(cosmic_tmp, colsplit(cosmic_tmp$mut, pattern = ">", names = c('ref', 'alt'))))
		## Meto una columna que indique en qué frequencia se introdujo esa mutación.
		cosmic_tmp$freqTT <- as.integer(tmp[[1]][8])
	
		## Se va creando una única tabla.
		cosmic <- rbind(cosmic,cosmic_tmp)
	}
	## Añadimos la columna de software y de identificación de la mutación para los gráficos. Dividimos la frequencia 
	## entre 100 para que coincida con la medida del resto de métodos. Aquí si incluimos alt en la comparación
	## también habría que meterlo.
	cosmic$software <- "True Mutations"
	cosmic$merged <- do.call(paste, c(cosmic[c("chr", "pos")], sep = "_"))
	cosmic$freqTT <- cosmic$freqT/100
	cosmic$freqT <- as.numeric(cosmic$gq_tumor)/(as.numeric(cosmic$gt_tumor) + as.numeric(cosmic$gq_tumor))
	cosmic$DepthTumor <- cosmic$gt_tumor + cosmic$gq_tumor
	cosmic_good <- subset(cosmic,gq_tumor != 0)

	return(cosmic_good)

}

######################
## Diagrama de Venn ##
######################

venn <- function(varTable){

venn.diagram(
	x = list(
		Strelka = varTable$merged[varTable$software == "Strelka"],
		JointSNVmix = varTable$merged[varTable$software == "JointSNVmix"],
		GATK = varTable$merged[varTable$software == "UnifiedGenotyper"],
		Samtools = varTable$merged[varTable$software == "Samtools"]
		),
	filename = "VennJointStrelkagatk.tiff",
	lwd = 2,
	height= 4250,
	width= 4370,
	fill = c("cornflowerblue", "darkorchid1","orange","darkgreen"), 
	alpha = 0.5,
	label.col = "white",
	cex = 2,
	fontfamily = "serif",
	fontface = "bold",
	cat.col = c("cornflowerblue", "darkorchid1","orange","darkgreen"), #
	cat.cex = 2,
	cat.fontfamily = "serif",
	cat.fontface = "bold",
	cat.dist = c(0.08,0.08,0.08,0.08),#
	cat.pos = 0
	);

} 

######################################
## Análisis individual de cobertura ##
######################################
summary_depth <- function(varTable){

	# mean, sd, min y max para las variantes llamadas por cada uno de los software.
	summaryDepthT <- by(varTable,as.factor(varTable$software),summarise,mean_depthT=mean(DepthTumor,na.rm=TRUE),sd_depthT = sd(DepthTumor,na.rm=TRUE),min_depthT = min(DepthTumor,na.rm=TRUE),max_depthT = max(DepthTumor,na.rm=TRUE))
	summaryDepthC <- by(varTable,as.factor(varTable$software),summarise,mean_depthC=mean(DepthControl,na.rm=TRUE),sd_depthC = sd(DepthControl,na.rm=TRUE),min_depthC = min(DepthControl,na.rm=TRUE),max_depthC = max(DepthControl,na.rm=TRUE))
	summaryDepth <- list(summaryDepthT,summaryDepthC)
	# Histograma para cobertura para todos los metodos.
	pdf("histDepth_bwa.pdf")
	 print(ggplot(varTable, aes(x=DepthTumor)) + geom_histogram(aes(y=..density..,fill=..count..)) + xlim(c(0,200)) + scale_fill_continuous(low="green",high="red", name="# mutation") + geom_density() + labs(title="Depth distribution on Tumor sample") + facet_wrap(~ software))
	 print(ggplot(varTable, aes(x=DepthControl)) + geom_histogram(aes(y=..density..,fill=..count..)) + xlim(c(0,200)) + scale_fill_continuous(low="green",high="red", name="# mutation") + geom_density() + labs(title="Depth distribution on Tumor sample") + facet_wrap(~ software))
	dev.off()
	return(summaryDepth)
}



#######################################
## Análisis individual para freq ##
#######################################
summary_freq <- function(varTable){
	# mean, sd, min y max para las variantes llamadas por cada uno de los software.
	summaryFreqT <- by(varTable,as.factor(varTable$software),summarise,mean_freqT=mean(freqT,na.rm=TRUE),sd_freqT = sd(freqT,na.rm=TRUE),min_freqT = min(freqT,na.rm=TRUE),max_freqT = max(freqT,na.rm=TRUE))
	summaryFreqC <- by(varTable,as.factor(varTable$software),summarise,mean_freqC=mean(freqC,na.rm=TRUE),sd_freqC = sd(freqC,na.rm=TRUE),min_freqC = min(freqC,na.rm=TRUE),max_freqC = max(freqC,na.rm=TRUE))

	summaryFreq <- list(summaryFreqT,summaryFreqC)

	png("boxplot.png",width=880)
	plot1 <- ggplot(varTable, aes(x=software, y=freqT, fill=type)) + scale_fill_manual(values=c("#698DCD","#69CD79")) + guides(fill=FALSE) + geom_boxplot(alpha=0.5) + stat_summary(fun.y=mean, geom="point", shape=3, size=3) + scale_x_discrete(limits=c("Samtools","UnifiedGenotyper","Strelka","JointSNVmix")) + labs(list(title="Boxplot for mutation frequencies in tumor sample", y="Frequency"))
	plot2 <- ggplot(varTable, aes(x=software, y=freqC, fill=type)) + scale_fill_manual(values=c("#698DCD","#69CD79")) + guides(fill=FALSE) + geom_boxplot(alpha=0.5) + stat_summary(fun.y=mean, geom="point", shape=3, size=3) + scale_x_discrete(limits=c("Samtools","UnifiedGenotyper","Strelka","JointSNVmix")) + labs(list(title="Boxplot for mutation frequencies in control sample", y="Frequency"))
	(plot3 <- grid.arrange(plot1, plot2, ncol=2))
	dev.off()

	# print(t.test(variantStrelka$freqT,variantJoint$freqT))
	# print(t.test(variantSamtools$freqT,variantGatk$freqT))
	# print(t.test(varTable$freqT[all_variants1$type == "joint"],all_variants1$freqT[all_variants1$type == "independent"]))

	#########################################################################
	## Histogramas para freq: todo junto sofwares y snps para una muestra. ##
	##########################################################################

	cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
	pdf("histFreq_bwa.pdf")
	print(ggplot(varTable, aes(x=freqT, fill=software, colour=software)) + scale_fill_brewer(palette="Accent") + scale_colour_brewer(palette="Accent") + geom_density(alpha=0.3) + labs(list(title="Density Function for mutation frequencies", x="Frequencies"))) 
	dev.off()

	# #############################################################################################
	# ## Gráfico de barras para frecuencia, menor de 0,3 y mayor de 0,3 para todos los software  ##
	# #############################################################################################

	## Importante para calcular porcentajes para cada una de las longitudes.
	df <- ddply(varTable, .(software,het), summarise, p = length(freqT))
	df <- ddply(df, .(software), transform, p2 = p/sum(p))
	
	pdf("bars_freq_software.pdf")
	 print(ggplot(df, aes(x = software,fill=het)) +  geom_bar(aes(y = p2), position="dodge", stat="identity") + scale_y_continuous(labels = percent_format()))
	dev.off()

	return(summaryFreq)
}

fp <- function(x,file_tp,file_tn){
	tp<-0
	fp<-0
	fn<-0
	tn<-0
	
	for (i in 1:nrow(x)){
		if (x[i,"merged"] %in% file_tp$merged) {
			tp=tp+1
		}else{
			fp=fp+1
		}
	}		
	fn <- nrow(file_tp) - tp
	tn <- nrow(file_tn) - fp
	s <- tp / (tp+fn)
	e <- tn / (tn+fp)

	result <- data.frame(TP=tp,FP=fp,FN=fn,TN=tn,Total_tp=nrow(file_tp),Total_tn=nrow(file_tn),Specifity=e,Sensitivity=s)

	return(result)
}

# #######################################################################################
# ####  Tumores virtuales: Análisis tumor con mutaciones en 6 frecuencias diferentes ####
# #######################################################################################

virtual_freq <- function(all_variants,cosmic){
	## Obtengo subsets de las tablas para quedarme sólo con lo que necesito para hacer el gráfico.
	## Las uno todas en una misma tabla.
	all_variants <- subset(all_variants,select=c(merged,freqT,software)) 
	cosmic <- subset(cosmic,select=c(merged,freqT,software))
	
	virtual_all <- rbind(all_variants,cosmic)
	
	## Gráfico de puntos que representa la frequencia de las mutaciones que llama cada software.
	pdf("virtual_freq.pdf")
	print(ggplot(virtual_all, aes(x=software, y=freqT)) + geom_point(aes(colour=factor(software))) + scale_fill_brewer(palette="Accent") + theme(axis.text.x = element_text(size = 7.5,angle=75, vjust=0.5), strip.text.x = element_text(size=6.5))+ scale_x_discrete(limits=c("True Mutations","Hapmap SNP","Samtools","UnifiedGenotyper","Strelka","JointSNVmix")) + guides(colour=FALSE) + labs(list(title="Mutation VAF values detected by each method", x="Methods", y="AAF")))
	dev.off()

}

###############################################
#### Cálculo de sensibilidad y especificidad ##
###############################################


fp_count_freq <- function(varTable,cosmic){
	control0 <- subset(varTable,var_sample_id == "124_virtual")
	control5 <- subset(varTable,var_sample_id == "124_0.05")
	control10 <- subset(varTable,var_sample_id == "124_0.1")

	cosmic_tp <- cosmic
	cosmic_tn <- data.frame(chr=character(),start=character(),end=character(),mut=character(),strand=character(),ref=character(),alt=character(),freqT=integer(),software=(character),merged=character())
	
	cosmic_tp_5 <- subset(cosmic, freqTT != 0.4)
	cosmic_tn_5 <- subset(cosmic, freqTT == 0.4)

	cosmic_tp_10 <- subset(cosmic, freqTT != 0.2)
	cosmic_tn_10 <- subset(cosmic, freqTT == 0.2)

	fp_count_c0_by <- by(control0,list(as.factor(control0$software)),fp,file_tp=cosmic_tp,file_tn=cosmic_tn)
	fp_count_c0 <- do.call(rbind,fp_count_c0_by)
	fp_count_c0$fp_control <- 0
	fp_count_c0$sample <- rownames(fp_count_c0)
	rownames(fp_count_c0) <- c(1:length(rownames(fp_count_c0)))

	fp_count_c5_by <- by(control5,list(as.factor(control5$software)),fp,file_tp=cosmic_tp_5,file_tn=cosmic_tn_5)
	fp_count_c5 <- do.call(rbind,fp_count_c5_by)
	fp_count_c5$fp_control <- 0.05
	fp_count_c5$sample <- rownames(fp_count_c5)
	rownames(fp_count_c5) <- c(1:length(rownames(fp_count_c5)))

	fp_count_c10_by <- by(control10,list(as.factor(control10$software)),fp,file_tp=cosmic_tp_10,file_tn=cosmic_tn_10)
	fp_count_c10 <- do.call(rbind,fp_count_c10_by)
	fp_count_c10$fp_control <- 0.10
	fp_count_c10$sample <- rownames(fp_count_c10)
	rownames(fp_count_c10) <- c(1:length(rownames(fp_count_c10)))

	fp_count_freq <- rbind(fp_count_c0,fp_count_c5,fp_count_c10)

	return(fp_count_freq)
}

fp_count_cob <- function(all_variants,cosmic){

}

got <- function(varTable,cosmic){
	varTable <- subset(varTable,software != "True Mutations")

	control0 <- subset(varTable,var_sample_id == "124_virtual")
	control5 <- subset(varTable,var_sample_id == "124_0.05")
	control10 <- subset(varTable,var_sample_id == "124_0.1")

	cosmic_tp <- cosmic
	cosmic_tn <- data.frame(chr=character(),start=character(),end=character(),mut=character(),strand=character(),ref=character(),alt=character(),freqT=integer(),software=(character),merged=character())
	
	cosmic_tp_5 <- subset(cosmic, freqTT != 0.2)
	cosmic_tn_5 <- subset(cosmic, freqTT == 0.2)

	cosmic_tp_10 <- subset(cosmic, freqTT != 0.4)
	cosmic_tn_10 <- subset(cosmic, freqTT == 0.4)


	got_by <- by(control0,list(as.factor(control0$software)),got_notgot,file_tp=cosmic_tp)
	got <- do.call(rbind,got_by)

	pdf("point_freq_depth.pdf",width=15,height=15)
	print(ggplot(got, aes(x=DepthTumor, y=freqT)) + geom_point(aes(colour=factor(got_it))) + stat_smooth(aes(group=got_it),method="glm",formula= y~ x,size=1,fullrange=FALSE,se=FALSE) + scale_fill_brewer(palette="Accent")+ ylim(c(0,0.50))+xlim(c(0,300)) + labs(list(title="AAF vs DC for mutations reported and not reported", x="Depth of Coverage", y="AAF")) + facet_wrap(~software))
	dev.off()

}

got_notgot <- function(x,file_tp){
	got <- data.frame(var_sample_id=character(),chr=character(),pos=character(),ref=character(),alt=character(),gt_control=character(),gq_control=character(),gt_tumor=character(),gq_tumor=character(),merged=character(),DepthControl=character(),DepthTumor=character(),freqT=character(),freqC=character(),software=character(),type=character(),het=character())
	not_got <- data.frame(var_sample_id=character(),chr=character(),pos=character(),ref=character(),alt=character(),gt_control=character(),gq_control=character(),gt_tumor=character(),gq_tumor=character(),merged=character(),DepthControl=character(),DepthTumor=character(),freqT=character(),freqC=character(),software=character(),type=character(),het=character())

	for (i in 1:nrow(file_tp)){
		if (file_tp[i,"merged"] %in% x$merged) {
			got <- rbind(got,file_tp[i,])
		}else{
			not_got <- rbind(not_got,file_tp[i,])
		}
	}

	got$got_it <- "got"
	got$software <- x$software[1]
	not_got$got_it <- "not_got_it"
	not_got$software <- x$software[1]
	got_q <- rbind(got,not_got)

	return(got_q)	
}


