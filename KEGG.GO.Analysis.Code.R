library(clusterProfiler)

#genelist is a text list of ncbi-protein IDs of Ostrinia furnacalis ("ofu") orthologs for outliers as determined by MMSeq2 
genelistOF<-read.table('genelist.txt', header=T, sep="\t")

gene1<-as.character(genelistOF$V1)


#KEGG enrichment

kk1<-enrichKEGG(
gene1,
organism = "ofu",
keyType = "ncbi-proteinid",
pvalueCutoff = 0.1,
pAdjustMethod = "BH",
minGSSize = 5,
maxGSSize = 100,
qvalueCutoff = 0.1,
use_internal_data = FALSE)

#GO enrichment

#geneDMout is a text list of Uniprot IDs for outliers of Drosophila melanogaster ("DM") orthologs as determined by MMSeq2 
geneDMout<-read.table('genelMout.txt', header=T, sep="\t")
gene4<-as.character(geneDMout$V1)


#geneDMUni is a text list of Uniprot IDs for all Drosophila melanogaster ("DM") orthologs as determined by MMSeq2 "Gene universe"
geneDMUni<-read.table('Universe.txt', header=F, sep="\t")
uni3<-as.character(geneDMUni<$V1)


library(org.Dm.eg.db)

ego1 <- enrichGO(gene        = gene4,
               OrgDb         = org.Dm.eg.db,
               ont           = "ALL", keyType = "UNIPROT",
               pAdjustMethod = "fdr", universe=uni3,
               pvalueCutoff  = 0.10, 
               minGSSize = 10, 
               maxGSSize = 500, 
               qvalueCutoff  = 0.10,
               readable      = TRUE)

