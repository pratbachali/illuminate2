load("~/Dropbox (AMPEL BioSolutions)/Bioinformatics/studies_sle/adult/GSE72747/T24_T12_T0_LIMMA/final_results_de/GSE72747_samples_T0_T12_T24_affy_ba.RData")
ls()
sig_uniq <- GSE72747.affy.ba.sig
#sig_uniq_UP <- sig_uniq[grep("UP", sig_uniq$fc_direction),]
#sig_uniq_DOWN <- sig_uniq[grep("DOWN", sig_uniq$fc_direction),]
sig_uniq <- sig_uniq[ order(sig_uniq$LIMMA.pVal.adj), ] #12,953
sig_uniq <- sig_uniq[ !duplicated(sig_uniq$geneSymbol), ] # 6975 


# Load Affy CDF
load("~/Dropbox (AMPEL BioSolutions)/Bioinformatics/studies_sle/adult/GSE72747/esets/eset_affy_gcrma_all_samples_ampel.RData")
ls()
eset_affy <- eset_gcrma_ampel
#rm(eset)
patient_data <- pData(eset_affy)
table(eset_affy$cohort) # 10 T0, 10 T12

load("~/Dropbox (AMPEL BioSolutions)/Bioinformatics/studies_sle/adult/GSE72747/esets/eset_brainarray_gcrma_all_samples_ampel.RData")
ls()
eset_ba <- eset_brainarray_gcrma_ampel
#rm(eset)
pDATA <- pData(eset_ba)
table(eset_ba$cohort)

cohort = as.character(eset_affy$cohort)
eset_affy$friendlyName <- paste( cohort, substr( sampleNames(eset_affy),7,10),eset_affy$SLEDAI, sep="." )

colnames <- eset_affy$friendlyName
colnames(eset_affy) <- colnames
exprs_heatmap <- matrix(nrow=nrow(sig_uniq),ncol=ncol(eset_affy))


##################################################################
# GENERATE THE EXPRESSION VALUES FOR THE HEATMAP
for(i in 1:nrow(sig_uniq)){
  if (grepl("^Affy.",sig_uniq[i,"CDF"])) {
    index <- which(fData(eset_affy)[,"probe"] == as.character(sig_uniq[i,"probe"]))
    if (length(index)>0) {
      exprs_heatmap[i,] = exprs(eset_affy)[index,]
      print( paste("Mapping entry ",i,": ", fData(eset_affy)[index,"geneSymbol"], " to Affy eset", sep="") )
    } else {
      print( paste("Mapping entry ",i,": NOT FOUND IN AFFY ESET", sep="") )
    }
  } else {
    index <- which(fData(eset_ba)[,"probe"] == as.character(sig_uniq[i,"probe"]))
    if (length(index)>0) {
      exprs_heatmap[i,] = exprs(eset_ba)[index,]
      print( paste("Mapping entry ",i,": ", fData(eset_ba)[index,"geneSymbol"], " to BrainArray eset", sep="") )
    } else {
      print( paste("Mapping entry ",i,": NOT FOUND IN BRAINARRAY ESET", sep="") )
    }
  }
}  
#Identify row numbers which contain probes not found in the esets
index <- which(is.na(exprs_heatmap[,1]))

# Rename the columns to patient arrays
colnames(exprs_heatmap) <- colnames(eset_affy)

# Find Affy CDF entries
index_affy <- which(grepl("^Affy.",sig_uniq$CDF)) # 2331
exprs_affy <- exprs_heatmap[index_affy,]
rownames(exprs_affy) <- sig_uniq[ index_affy, "probe" ]
fData_affy <- fData(eset_affy)
fData_affy$probe <- as.character(fData_affy$probe)
exprs_affy <- data.frame(exprs_affy)
exprs_affy$probe <- rownames(exprs_affy)
exprs_affy <- merge( exprs_affy, fData_affy[,1:2] )
rownames(exprs_affy)  <- exprs_affy$geneSymbol


colnames(exprs_affy)
exprs_affy <- exprs_affy[,c(-1,-ncol(exprs_affy))]

# STOP! Adjust the columns to remove according to this data set
#exprs_affy <- exprs_affy[,c(-1,-11)]

# Find BrainArray entries
index_ba <- which(!grepl("^Affy.",sig_uniq$CDF)) # 48
exprs_ba <- exprs_heatmap[index_ba,]
rownames(exprs_ba) <- sig_uniq[ index_ba, "probe" ]
fData_ba <- fData(eset_ba)
fData_ba$probe <- as.character(fData_ba$probe)
exprs_ba <- data.frame(exprs_ba)
exprs_ba$probe <- rownames(exprs_ba)
exprs_ba <- merge( exprs_ba, fData_ba[,1:2] )
#rownames(exprs_ba)  <- exprs_ba$probe
rownames(exprs_ba)  <- exprs_ba$geneSymbol

exprs_ba <- exprs_ba[,c(-1,-ncol(exprs_ba))]
exprs_heatmap <- rbind( exprs_affy, exprs_ba )
exprs_heatmap <- data.matrix(exprs_heatmap)

#index <- which(colnames(exprs_heatmap)=="T0.9396.17" | colnames(exprs_heatmap)=="T12.9397.2" | colnames(exprs_heatmap)=="T0.9402.14" | colnames(exprs_heatmap)=="T12.9403.2" | colnames(exprs_heatmap)=="T0.9411.17" | colnames(exprs_heatmap)=="T12.9412.2" | colnames(exprs_heatmap)=="T0.9417.14" | colnames(exprs_heatmap)=="T12.9418.4")
#exprs_heatmap_one <- exprs_heatmap[,index]
#index <- which(colnames(exprs_heatmap)=="T0.9399.13" | colnames(exprs_heatmap)=="T12.9400.8" | colnames(exprs_heatmap)=="T0.9408.19" | colnames(exprs_heatmap)=="T12.9409.14" | colnames(exprs_heatmap)=="T0.9420.16" | colnames(exprs_heatmap)=="T12.9421.14")
#exprs_heatmap_two <- exprs_heatmap[,index]

setwd("~/Dropbox (AMPEL BioSolutions)/Bioinformatics/studies_sle/adult/GSVA/adult/gsva_lists")

Increases_T12_T24 <- read.delim("GSE72747_Timecourse_forGSVA_CellSubsets_Increased_T12_T24_9Jan2018.txt")


for(i in 1:ncol(Increases_T12_T24)){
  testlist1 <- as.character(Increases_T12_T24[,1])
  testlist1 <- testlist1[1:34]
  testlist1 <- list('Monocytes'=c(testlist1))
  testlist2 <- as.character(Increases_T12_T24[,2])
  testlist2 <- testlist2[1:14]
  testlist2 <- list('Myeloid'=c(testlist2))
  testlist3 <- as.character(Increases_T12_T24[,3])
  testlist3 <- testlist3[1:8]
  testlist3 <- list('LDG'=c(testlist3))
  testlist4 <- as.character(Increases_T12_T24[,4])
  testlist4 <- testlist4[1:13]
  testlist4 <- list('Hematopoietic'=c(testlist4))
  testlist5 <- as.character(Increases_T12_T24[,5])
  testlist5 <- testlist5[1:3]
  testlist5 <- list('Ery_Platelets'=c(testlist5))
}

testList <- c(testlist1,testlist2,testlist3,testlist4,testlist5)


rm(testlist1,testlist2,testlist3,testlist4,testlist5)

n <- names(testList)
uniqueList <- lapply(testList, unique)

makeSet <- function(geneIds, n) {
  GeneSet(geneIds, geneIdType=SymbolIdentifier(), setName=n)
}

library(dplyr)
library(GSEABase)
library(GSVA)
gsList <- gsc <- mapply(makeSet, uniqueList[], n)
newlist = gsList %>% lapply(function(x) trimws(x@geneIds))
gsc <- GeneSetCollection(gsList)

#### AFTER CREATING THE GENESET COLLECTION OBJECT NOW WE CAN PERFORM GSVA FUNCTION AND 
#### GET EXPRESSION SET WHICH IS TRANSFORMED FROM GENE BY SAMPLE TO GENESET BY SAMPLE

eset <- exprs_heatmap

### emtab2713_pbmc #####
gse72747_gsva <- gsva(eset, newlist, method = c("gsva"), rnaseq=FALSE,
                      min.sz=1, max.sz=Inf, verbose=TRUE)$es.obs


setwd("~/Dropbox (AMPEL BioSolutions)/Bioinformatics/studies_sle/adult/GSVA/adult/Time-series/GSE72747/Enrichment_Scores_Lists")
write.csv(gse72747_gsva,file="GSE72747_WB_MixedSex_MixedRac_10T0_10T12_10T24_Increased_CellSubsets_Enrichment_scores.csv")


########Heatmap gsva#################
exprs_hm <- gse72747_gsva
index <- which(is.na(rowSums(exprs_hm)))
if (length(index) > 0) { exprs_hm <- exprs_hm[-index,] }

# Here we'll use Euclidean distances
distance <- dist(t(exprs_hm),method="euclidean")
hc <- hclust(distance)

dist2 = dist2 <- function(x, ...)
  dist(x,method="euclidean")
label.size = 1.2
samples = colnames(exprs_hm)

hmcol <- colorRampPalette(c("blue","grey","red"))(256)


library(d3heatmap)
#row.names(exprs_hm) <- exprs_hm$geneSet




d3heatmap(data.matrix(exprs_hm), colors = hmcol, scale = "row", dendrogram = "both", hclustfun = hclust, distfun = dist2,
          k_row = 3, k_col = 3, height = 220, width = 700, show_grid = FALSE)




