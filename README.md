# GOenrichment_GOstats
Testing GO enrichment with the R package GOstats
###################################
source("https://bioconductor.org/biocLite.R")
biocLite("GOstats")
biocLite("biomaRt")
biocLite("GO.db")
biocLite("topGO")
biocLite("GOSemSim")

library(GO.db)
library(topGO)
library(GOSemSim)
library(biomaRt)

mart = useMart(biomart="plants_mart",host="plants.ensembl.org", dataset = "zmays_eg_gene")
univ.geneID <- getBM(attributes=c("ensembl_gene_id", "entrezgene"), mart = mart) # 40481
univ.geneID = arrange(univ.geneID, ensembl_gene_id)
head(univ.geneID3)
## remove genes with no corresponding Entrez Gene ID
univ.geneID2 <- univ.geneID[!is.na(univ.geneID[,2]),] # 14142 
## remove duplicated Entrez Gene ID
univ.geneID3 <- univ.geneID2[ !duplicated(univ.geneID2[,2]),] # 13630
##Get GO terms
univ.geneID4<-getBM(attributes=c("entrezgene","goslim_goa_accession","name_1006","namespace_1003","go_linkage_type"),filters='entrezgene',mart=mart,values=univ.geneID3$entrezgene)
#listAttributes(mart)
#grep("name", listAttributes(mart)[,1], value = TRUE)

###Make dataframe for GOStats
goframeData <- data.frame(go_id = univ.geneID4$goslim_goa_accession, Evidence = univ.geneID4$go_linkage_type, gene_id = univ.geneID4$entrezgene,stringsAsFactors=F)
head(goframeData)

## Make my geneID
f2_ra1_1 = filter(f2_df, series == "mut_series", size == "1mm", q2 == "ra1", significant == "yes")
f2_ra1_1_cut = f2_ra1_1[, 1:9]
colnames(f2_ra1_1_cut)[1] <- "ensembl_gene_id"
f2_ra1_1_cut_A = arrange (f2_ra1_1_cut, ensembl_gene_id)
univ.geneID3_A = arrange (univ.geneID3, ensembl_gene_id)
f2_ra1_1_GO <- merge(f2_ra1_1_cut_A, univ.geneID3_A, by ="ensembl_gene_id") # 620
f2_ra1_1_GO

## GO enrichment
library("GOstats")
library("GOSemSim")
##Prepare GO to gene mappings
library(GSEABase)
goFrame=GOFrame(goframeData,organism="Zea mays")
goAllFrame=GOAllFrame(goFrame)
gsc <- GeneSetCollection(goAllFrame, setType = GOCollection())
params <- GSEAGOHyperGParams(name="Zea mays GO", geneSetCollection=gsc, geneIds = f2_ra1_1_GO[,10], universeGeneIds = univ.geneID4$entrezgene, ontology = "BP", pvalueCutoff = 0.06, conditional = TRUE, testDirection = "over")
BP <- hyperGTest(params)
summary(BP)[,c(1,2,7)]
