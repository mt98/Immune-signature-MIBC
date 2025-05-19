##Script for WGCNA to identify genes for a signature for immunotherapy in muscle invasive bladder cancer

##Load libraries
library("WGCNA")
library("cluster")
library(AnnotationDbi)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()
library(stringr)
library("clusterProfiler")


##Used the PURE01 cohort
###Clinical data
PURE1_clinical<-read_excel("databasePURELM_07.13.20_survival data.xlsx", sheet = "Patients in our PURE01 cohort",col_names=T)
PURE1_clinical<-as.data.frame(PURE1_clinical)
rownames(PURE1_clinical)<-as.character(PURE1_clinical$PURE01.ID)
###censoring at 24 months
censorship<-24
PURE1_clinical_event_censorshipo  <- ifelse(  PURE1_clinical$Time.to.recurrence <= censorship & PURE1_clinical$RELAPSE == 1, 1 ,0)
PURE1_clinical_time_censorshipo   <- ifelse( PURE1_clinical_event_censorshipo == 0 &PURE1_clinical$Time.to.recurrence >= censorship, censorship , PURE1_clinical$Time.to.recurrence )   
PURE1_clinical               <- cbind( PURE1_clinical, PURE1_clinical_time_censorshipo, PURE1_clinical_event_censorshipo )
colnames(PURE1_clinical   )[ (ncol(PURE1_clinical  )-1):ncol(PURE1_clinical  ) ] <- c( "censored_time", "censored_status" )



###Look at association with PD-L1, immune infiltration and estimate scores


datTraits<-read_excel("PURE01_supplementary_information.20210903.xlsx", sheet = "Table S1",col_names=T,skip=1)
datTraits<-as.data.frame(datTraits)
datTraits1<-datTraits[,c(6)]
names(datTraits1)<-datTraits$Pre.therapy.RNAseq.file


datTraits<-read_excel("PURE01_supplementary_information.20210903.xlsx", sheet = "Table S4",col_names=T,skip=1)
datTraits<-as.data.frame(datTraits)
 datTraits2<-datTraits[,c(2:5)]
  rownames(datTraits2)<-datTraits$id

datTraits3<-merge(datTraits1,datTraits2,by="row.names")
rownames(datTraits3)<-datTraits3$Row.names
datTraits3<-datTraits3[,-c(1)]


datTraits<-read_excel("PURE01_supplementary_information.20210903.xlsx", sheet = "Table S6",col_names=T,skip=1)
datTraits<-as.data.frame(datTraits)
 datTraits4<-datTraits[,c(4:13)]
  rownames(datTraits4)<-datTraits$id





##look at only the 82 pretreatment biopsies of patients
PURE1_expression<-PURE_expression[,intersect(rownames(PURE_data),colnames(PURE_expression))]
PURE1_expression<-log2(PURE1_expression+1)


###look at only mRNA based genes, 20,000
library(biomaRt)
ensembl = useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")

res <- getBM(attributes=c("ensembl_gene_id","gene_biotype"),filters = c("ensembl_gene_id","biotype"), values=list(rownames(PURE1_expression),"protein_coding"), mart=ensembl)
PURE1_expression<-PURE1_expression[res$ensembl_gene_id,] 


###Filter out genes that are lowly or not expressed
###Filter out genes which are in the bottom 50th percentile of the network
variance<-apply((PURE1_expression),1,mad)
PURE1_filtered<-PURE1_expression[which(variance>quantile(variance,0.50)),] ##9707 genes


###build a network using WGCNA

###Remove outlier samples based on hierarchical clustering
datExpr0<-as.data.frame(t(PURE1_filtered))
sampleTree = hclust(dist(datExpr0), method = "average");
pdf("removal.pdf")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)
dev.off()
##Build network
###Inorder to develop a weighed correlation network we need to find the optimal value of the soft threshold



powers = c(c(1:10), seq(from = 12, to=20, by=2))

# Call the network topology analysis function

##Draw an abline for connectivity to make sure it is not zero
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5,networkType="signed hybrid",corFnc = "bicor", corOptions = list (maxPOutliers = 0.05))
pdf("R2.pdf")
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,col="red")
dev.off()

# Mean connectivity as a function of the soft-thresholding power 
cex1 = 0.9;
pdf("connectivity.pdf")
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n", main = paste("Mean connectivity")) 
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
abline(h=0.0, col="red")
dev.off()


softPower = 4
adjacency = adjacency(datExpr0, power = softPower,type="signed hybrid", corFnc = "bicor",corOptions = list (maxPOutliers = 0.05))

# Turn adjacency into topological overlap

TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM

geneTree = hclust(as.dist(dissTOM), method = "average")

# Plot the resulting clustering tree (dendrogram)
pdf("genetree.pdf")
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity")
dev.off()

# We like large modules, so we set the minimum module size relatively high:

minModuleSize = 30
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                deepSplit = 2, pamRespectsDendro = FALSE,
                minClusterSize = minModuleSize)

dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)


pdf("dendro.pdf")
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()

# Calculate eigengenes

MEList = moduleEigengenes(datExpr0, colors = dynamicColors)
MEs = MEList$eigengenes

# Calculate dissimilarity of module eigengenes

MEDiss = 1-bicor(MEs)

# Cluster module eigengenes

METree = hclust(as.dist(MEDiss), method = "average")

# Plot the result
MEDissThres = 0.25
pdf("clustering.pdf")
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
abline(h=MEDissThres, col = "red")
dev.off()


merge = mergeCloseModules(datExpr0, dynamicColors, cutHeight = MEDissThres, verbose = 3,corFnc = "bicor",corOptions = list (maxPOutliers = 0.05))

mergedColors = merge$colors

mergedMEs = merge$newMEs

pdf("geneDendro.pdf")
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off ()


moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs

nGenes = ncol(datExpr0)
nSamples = nrow(datExpr0)

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr0, moduleColors)$eigengenes
MEs = orderMEs(MEs0)


##Look at correlation of MEs with  PD-L1, immune infiltration and estimate scores

moduleTraitCor = (cor(MEs, datTraits3[rownames(MEs),], use = "p"))
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
#Displaying correlations and its p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)


#Displaying the correlation values in a heatmap plot
pdf("moduletrait1.pdf")
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits3),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               setStdMargins = FALSE,
               cex.lab.x = 0.8,cex.lab.y=0.8,
               zlim = c(-1,1),
               main = paste("Module-Immune-biomarker relationships"))
dev.off()


moduleTraitCor1 = (cor(MEs, datTraits4[rownames(MEs),], use = "p"))
moduleTraitPvalue1 = corPvalueStudent(moduleTraitCor1, nSamples)

textMatrix1 =  paste(signif(moduleTraitCor1, 2), "\n(",
                    signif(moduleTraitPvalue1, 1), ")", sep = "")
dim(textMatrix1) = dim(moduleTraitCor1)

pdf("moduletrait2.pdf")
par(mar = c(6, 8.5, 3, 3))

#Displaying the correlation values in a heatmap plot
labeledHeatmap(Matrix = (moduleTraitCor1),
               xLabels = names(datTraits4),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               setStdMargins = FALSE,
               cex.lab.x = 0.8,cex.lab.y=0.8,
               zlim = c(-1,1),
               main = paste("Module-Immune cell relationships"))
dev.off()


###Association of module eigen genes with survival

association<-merge(PURE1_clinical,MEs,by="row.names")

###survival analysis
cox_gene_train   <- t(sapply( colnames(MEs),function(x){
res.cox<- coxph(Surv(censored_time, censored_status) ~ association[,x], data = association)
y<-summary(res.cox)
p.value<-signif(y$wald["pvalue"], digits=2)
 HR <-signif(y$coef[2], digits=2)
 HR.confint.lower <- signif(y$conf.int[,"lower .95"], 2)
 HR.confint.upper <- signif(y$conf.int[,"upper .95"],2)

                          res<-c(HR, HR.confint.lower, HR.confint.upper , p.value)
                          names(res)<-c("HR" , "HR.confint.lower", "HR.confint.upper",
                                        "p.value")
                          return(res)
 }))

###Look at module membership and connectivity of genes in different modules

modNames = names(MEs)

geneModuleMembership = as.data.frame(bicor(datExpr0, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
p<-intramodularConnectivity.fromExpr(datExpr0, moduleColors, 
              corFnc = "bicor", corOptions = "use = 'p'",
              weights = NULL,
              networkType = "signed hybrid", 
              scaleByMax = FALSE,
              ignoreColors = if (is.numeric(colors)) 0 else "grey",
              getWholeNetworkConnectivity = TRUE)
rownames(p)<-colnames(datExpr0)   

###Chose pink module based on inflammatory genes and correlation with immune score and PD-L1 along with association with survival

selectedgenes3<-names(datExpr0)[moduleColors=="pink"]

go<-enrichGO(
selectedgenes3,
'org.Hs.eg.db',
  keyType = "ENSEMBL",
  ont = "BP",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH")

write.table(summary(go),file="pinksigned.csv",sep=",",quote=F)

##Identify the most  connected or hub genes

pink<-geneModuleMembership[,which(grepl("Pink",colnames(geneModuleMembership)))] 
names(pink)<-rownames(geneModuleMembership)

selectedgenes_pink<-names(datExpr0)[moduleColors=="pink"]
pink<-pink[selectedgenes_pink]
significant_pink<-pink[abs(pink)>=(0.8)]
p_pink<-MMPvalue[names(significant_pink),2]
pink_connectivity<-p[(selectedgenes_pink),2]
names(pink_connectivity)<-(selectedgenes_pink)

connectivity_pink<-pink_connectivity[pink_connectivity>=quantile(pink_connectivity,0.90)]

length(intersect(names(connectivity_pink),names(significant_pink)))

pink_hub<-(intersect(names(connectivity_pink),names(significant_pink)))

hubpink <- getBM(attributes=c("hgnc_symbol","ensembl_gene_id"),filters = c("ensembl_gene_id"), values=(pink_hub), mart=ensembl)

###median expression of pink hub gnes
median<-apply(datExpr0[,pink_hub],1,median)

association_pink<-merge(PURE1_clinical,median,by="row.names")
res.cox<- coxph(Surv(censored_time, censored_status) ~ y, data = association_pink)
summary(res.cox)
####Look at stratification based on the median
stratify<-quantile(median,0.50)
category_pink<-ifelse(association_pink$y>=stratify,"high","low")
association_pink$category_pink<-category_pink

res.cox<- coxph(Surv(censored_time, censored_status) ~ category_pink, data = association_pink)
summary(res.cox)

###plot it 

fit<-survfit(Surv(censored_time,censored_status) ~ category_pink, data = association_pink)


par(pty="s")
pdf("survival_pink.pdf")
ggsurvplot(fit,censor=TRUE,pval=TRUE,linetype=c(1,3), risk.table.x.text = FALSE,tables.theme = clean_theme(),
                      risk.table = TRUE,palette=c("red","blue"),ylab="Recurrence free survival",break.time.by = 5,
                      xlab="Time in months", font.x = "bold", font.y = "bold", font.tickslab = "bold",
                      font.legend = "bold")
dev.off()



###path response


###Look at this signature in ABACUS

clinicaldata_ABACUS<-read.table(file="ABACUS metadata.csv",sep=",",header=T, stringsAsFactors=F)
rownames(clinicaldata_ABACUS)<-clinicaldata_ABACUS$ID
clinicaldata_ABACUS<-clinicaldata_ABACUS[clinicaldata_ABACUS$VISIT=="PRE",]
clinicaldata_ABACUS$PDL1<-ifelse(clinicaldata_ABACUS$PDL1_IC=="IC2+","positive","negative")

clinicaldata_ABACUS <- clinicaldata_ABACUS %>%
  mutate(response = case_when(clinicaldata_ABACUS$PCR=="Yes" ~ "CR",
 clinicaldata_ABACUS$MPR=="Yes"~ "PR",
clinicaldata_ABACUS$PCR=="No"&clinicaldata_ABACUS$MPR=="No" ~ "NR"
    ),
    response.2grp = case_when(clinicaldata_ABACUS$PCR=="Yes"~ "CR.PR",
 clinicaldata_ABACUS$MPR=="Yes" ~ "CR.PR",
clinicaldata_ABACUS$PCR=="No"&clinicaldata_ABACUS$MPR=="No"~ "NR"
    )
  )

###make it for the development of the ABACUS response
clinicaldata_ABACUS_compact<-clinicaldata_ABACUS[,c(1,33,34,35)]
rownames(clinicaldata_ABACUS_compact)<-clinicaldata_ABACUS_compact$ID
expression_ABACUS<-read.delim(file="ABACUS_n148_FPKM_19936_coding_genes.unique_gene_symbols.txt",sep="\t",header=T)
expression_ABACUS1<-expression_ABACUS[,4:151]
rownames(expression_ABACUS1)<-expression_ABACUS$gene_name
expression_ABACUS1<-log2(expression_ABACUS1+1)
colnames(expression_ABACUS1)<- sub("X","",colnames(expression_ABACUS1))





####pink module
pink_ABACUS<-expression_ABACUS1[,]
mean_pink_ABACUS<-apply(pink_ABACUS,2,mean)
category<-ifelse(median_pink_ABACUS>=quantile(mean_pink_ABACUS,0.50),"high","low")

ABACUS_category<-merge(clinicaldata_ABACUS_compact,category,by="row.names")
rownames(ABACUS_category)<-ABACUS_category$Row.names
colnames(ABACUS_category)[ncol(ABACUS_category)]<-"category"

ABACUS_category_high<-ABACUS_category[ABACUS_category$category=="high",]
table(ABACUS_category_high[,4]) 

ABACUS_category_low<-ABACUS_category[ABACUS_category$category=="low",]
table(ABACUS_category_low[,4]) 


table(ABACUS_category_high[,5]) 
table(ABACUS_category_high[,5]) 











