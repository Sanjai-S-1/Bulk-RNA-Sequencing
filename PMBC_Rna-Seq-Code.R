#Transcriptome Analysis of Peripheral Blood Mononuclear Cells in Pulmonary Sarcoidosis

#Importing Essential Libraries

library(tidyverse) # For data manipulation
library(dplyr) # For data wrangling
library(GEOquery) # Access and analyze Gene Expression Omnibus (GEO) data
library(edgeR) # RNA-Seq differential expression analysis
library(matrixStats) # statistical operations on matrices.
library(cowplot) # Combine multiple plots
library(DT) # for making interactive tables
library(plotly) # for making interactive plots
library(gt) # Data representation in interactive table
library(limma) # venerable package for differential gene expression using linear modeling
library(plotly) # For plots
library(RColorBrewer) # colors to make heatmaps
library(ComplexHeatmap) # For heatmap
library(GSEABase) #functions and methods for Gene Set Enrichment Analysis
library(Biobase) #base functions for bioconductor; required by GSEABase
library(GSVA) #Gene Set Variation Analysis
library(gprofiler2) #tools for accessing the GO enrichment results using g:Profiler web resources
library(clusterProfiler) # provides a suite of tools for functional enrichment analysis
library(msigdbr) # access to msigdb collections directly within R
library(enrichplot) # For making the GSEA enrichment plots
library(ggplot2) # For plots
library(ggforce) # Extends ggplot2
library(gplots) # For plots
library(ggrepel) # Avoids label overlap in ggplot2

# Reading Counts File
# Counts File contains Genes in rows and the Samples in columns

gene_counts <- read.csv("counts.txt", sep = "\t",  row.names = 1) 

# Getting Metadata from GEO database
# Metadata contains the additional information about the samples which is helpful for the analysis

geo_id <- "GSE192829"
gse <- getGEO(geo_id, GSEMatrix = T)
phenodata <- pData(phenoData(gse[[1]]))
head(phenodata)

targets <- phenodata[,c(1,2,8,44,48)] %>% as.tibble()
colnames(targets) <- c("title", "geo_accesion", "condition", "age", "sex")
targets$title <- c("H1","H2","H3","H4","H5","H6","H7","H8","H9","H10","H11", "PS1", "PS2", "PS3","PS4","PS5", "PS6", "PS7", "PS8")
colnames(gene_counts) <- c(targets$title)
targets$condition <- gsub(" control", "",targets$condition)
targets$condition <- gsub(" ", "_",targets$condition)
sampleLabels <- targets$title

# Log CPM transformation of Counts data 

myDGEList <- DGEList(gene_counts)
log2.cpm <- cpm(myDGEList, log=TRUE)

log2.cpm.df <- as_tibble(log2.cpm, rownames = "geneID")
colnames(log2.cpm.df) <- c("geneID", sampleLabels)

cpm <- cpm(myDGEList)
keepers <- rowSums(cpm>1)>=8 # Onlt the genes which have CPM greater than 1 in at least 8 samples are kept
myDGEList.filtered <- myDGEList[keepers,]


log2.cpm.filtered <- cpm(myDGEList.filtered, log=TRUE)
log2.cpm.filtered.df <- as_tibble(log2.cpm.filtered, rownames = "geneID")
colnames(log2.cpm.filtered.df) <- c("geneID", sampleLabels)

# Normalisation of Data

myDGEList.filtered.norm <- calcNormFactors(myDGEList.filtered, method = "TMM")
log2.cpm.filtered.norm <- cpm(myDGEList.filtered.norm, log=TRUE)
log2.cpm.filtered.norm.df <- as_tibble(log2.cpm.filtered.norm, rownames = "geneID")
colnames(log2.cpm.filtered.norm.df) <- c("geneID", sampleLabels)

##############################################################################

# Performing Principle Component Analysis [PCA]

# Assigning the groups
colnames(targets)[3] <- "group"
group <- targets$group
group <- factor(group)


pca.res <- prcomp(t(log2.cpm.filtered.norm), scale.=F, retx=T)
pc.var <- pca.res$sdev^2 # sdev^2 captures these eigenvalues from the PCA result
pc.per <- round(pc.var/sum(pc.var)*100, 1) 
pca.res.df <- as_tibble(pca.res$x) %>%
  add_column(sample = sampleLabels,
             group = group)



ggplot(pca.res.df, aes(x=PC1, y=PC2)) +
  coord_fixed(xlim=c(-40,60), ylim=c(-40,40)) +
  geom_point(aes(color = group, shape = group), size=4)  +
  #geom_label() +
  geom_mark_ellipse(aes(x=PC1, y=PC2, color = group)) +
  #stat_ellipse() + 
  xlab(paste0("PC1 (",pc.per[1],"%",")")) + 
  ylab(paste0("PC2 (",pc.per[2],"%",")")) +
  labs(title="PCA plot") +
  theme_bw() + 
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16))

#########################################################################################

# Identification of differential expressed genes

group <- factor(targets$group)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)


v.DEGList.filtered.norm <- voom(myDGEList.filtered.norm, design, plot = FALSE)
fit <- lmFit(v.DEGList.filtered.norm, design)
contrast.matrix <- makeContrasts(infection = Pulmonary_sarcoidosis - Healthy,
                                 levels=design)

fits <- contrasts.fit(fit, contrast.matrix)
ebFit <- eBayes(fits)
myTopHits <- topTable(ebFit, adjust ="BH", coef=1, number=40000, sort.by="logFC")
myTopHits.df <- myTopHits %>%
  as_tibble(rownames = "geneID")

#########################################################################################

# Violin plot of differential expressed genes

# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2fc respectively positive or negative)<br /><br /><br />
myTopHits$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP"
myTopHits$diffexpressed[myTopHits$logFC > 1 & myTopHits$P.Value < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN
myTopHits$diffexpressed[myTopHits$logFC < -1 & myTopHits$P.Value < 0.05]<- "DOWN"


myTopHits$gene_symbol <- rownames(myTopHits)
myTopHits$delabel <- NA

new_df <- myTopHits[order(abs(myTopHits$logFC) > 2 & myTopHits$P.Value < 0.05, decreasing = T),]
new_df[1:70, "delabel"] <- head(new_df$gene_symbol,70)

df <- new_df

ggplot(data = df, aes(x = logFC , y = -log10(P.Value ), col = diffexpressed, label = delabel)) +
  geom_vline(xintercept = c(-1, 1), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
  geom_point(size = 2) + 
  scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00"), # to set the colours of our variable  
                     labels = c("Downregulated", "Not significant", "Upregulated")) + # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
  coord_cartesian(ylim = c(0, 26), xlim = c(-12, 8)) + # since some genes can have minuslog10padj of inf, we set these limits
  labs(color = 'Severe', #legend_title, 
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) + 
  scale_x_continuous(breaks = seq(-12, 8, 2)) + # to customise the breaks in the x axis
  scale_y_continuous(breaks = seq(0,26,2)) +
  ggtitle('Differential expressed genes in Healthy vs Diseased individuals') + # Plot title 
  geom_text_repel(max.overlaps = Inf) # To show all labels 


##################################################################################

# Getting the list of differencially expressed genes

results <- decideTests(ebFit, method="global", adjust.method="BH", p.value=0.05, lfc=1) # determine differentially expressed genes
colnames(v.DEGList.filtered.norm$E) <- sampleLabels
diffGenes <- v.DEGList.filtered.norm$E[results[,1] !=0,]
diffGenes.df <- as_tibble(diffGenes, rownames = "geneID")

################################################################################

# Heat Map 

myheatcolors1 <- bluered(100)

clustRows <- hclust(as.dist(1-cor(t(diffGenes), method="pearson")), method="complete") 
clustColumns <- hclust(as.dist(1-cor(diffGenes, method="spearman")), method="complete")

plot(clustColumns)

colanno <- HeatmapAnnotation(group = targets$group, sex = targets$sex,
                             show_annotation_name = T, show_legend = T)
columnsplit <- as.vector(targets$group)
columnsplit <- gsub("_", " ", columnsplit)

scaled_mat = t(scale(t(diffGenes)))

h1 <- Heatmap(scaled_mat,
        cluster_rows = clustRows, cluster_columns = F,
        row_names_side = "left", col = myheatcolors1,
        show_column_dend = F, column_split = columnsplit, row_split = 2, row_gap = unit(c(2, 1), "mm"),
        column_gap = unit(c(2, 1), "mm"), border = T,
        row_title = NULL, name = "Z-score", row_names_gp = gpar(fontsize = 2),
        show_row_names = F, show_row_dend = F, show_column_names = T)


log2_val <- as.data.frame(myTopHits[,1], row.names = rownames(myTopHits)) 

log2_val <-  as.data.frame(log2_val[rownames(diffGenes),], row.names = rownames(diffGenes)) 
colnames(log2_val)[1] <- "logFC"
  

colours <- bluered(100)

h2 <- Heatmap(as.matrix(log2_val) , 
              cluster_rows = clustRows, name="logFC", col = colours, show_row_names = F)

h = h1+h2
h

##################################################################################

# GSEA analysis

gost.res_up <- gost(rownames(myModule_up), organism = "hsapiens", correction_method = "fdr")
gostplot(gost.res_up, interactive = T, capped = F)
gost.res_down <- gost(rownames(myModule_down), organism = "hsapiens", correction_method = "fdr")
gostplot(gost.res_down, interactive = T, capped = F)

c2cp <- read.gmt("c2.cp.kegg.v2023.1.Hs.symbols.gmt") 

hs_gsea_c2 <- msigdbr(species = "Homo sapiens", 
                      category = "C2") %>% 
  dplyr::select(gs_name, gene_symbol)  

mydata.df.sub <- dplyr::select(myTopHits.df, geneID, logFC)
mydata.gsea <- mydata.df.sub$logFC
names(mydata.gsea) <- as.character(mydata.df.sub$geneID)
mydata.gsea <- sort(mydata.gsea, decreasing = TRUE)

# run GSEA using the 'GSEA' function from clusterProfiler
myGSEA.res <- GSEA(mydata.gsea, TERM2GENE=hs_gsea_c2, verbose=FALSE)
myGSEA.df <- as_tibble(myGSEA.res@result)

# view results as an interactive table
datatable(myGSEA.df, 
          extensions = c('KeyTable', "FixedHeader"), 
          caption = 'Signatures enriched in leishmaniasis',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(3:10), digits=2)

# create enrichment plots using the enrichplot package
gseaplot2(myGSEA.res, 
          geneSetID = 1, 
          pvalue_table = FALSE,
          title = myGSEA.res$Description[1]) 

# add a variable to this result that matches enrichment direction with phenotype
myGSEA.df <- myGSEA.df %>%
  mutate(phenotype = case_when(
    NES > 0 ~ "disease",
    NES < 0 ~ "healthy"))

# Bubble plot
ggplot(myGSEA.df[50:100,], aes(x=phenotype, y=ID)) + 
  geom_point(aes(size=setSize, color = NES, alpha=-log10(p.adjust))) +
  scale_color_gradient(low="blue", high="red") +
  theme_bw()


################################################################################

# Gene ontology Analysis

genes_to_test <- rownames(myTopHits[abs(myTopHits$logFC) > 1.5,])

ego_ORA_90cutoff <- enrichGO(gene     = genes_to_test ,
                             OrgDb         = "org.Hs.eg.db",
                             keyType       = 'SYMBOL',
                             ont           = "ALL",
                             pAdjustMethod = "BH",
                             pvalueCutoff  = 0.01,
                             qvalueCutoff  = 0.05)

barplot(ego_ORA_90cutoff, split = "ONTOLOGY")+facet_grid(ONTOLOGY~., scale = "free")+ggtitle("Barplot for GO_ORA_90cuttoff_PC")

mytophits <- dplyr::filter(mytophits, by = abs(log)  )

library(org.Hs.eg.db)
hs <- org.Hs.eg.db
my.symbols <- c("ANKRD62P1-PARP4P3")
select(hs, 
       keys = my.symbols,
       columns = c("ENTREZID", "SYMBOL"),
       keytype = "SYMBOL")

###################################################################################

# KEGG PATHWAY ANALYSIS

res <- mytophits[,-3:-7] 

library(AnnotationDbi)
library(org.Hs.eg.db)


res$Entrez = mapIds(org.Hs.eg.db, 
                    keys = res$geneID,
                    keytype = "SYMBOL",
                    column = "ENTREZID")

head(res)


library(pathview)
library(gage)
library(gageData)

res <- res[complete.cases(res), ]

foldchanges <- res$logFC
names(foldchanges) <- res$Entrez

head(foldchanges)

data("kegg.sets.hs")
data("sigmet.idx.hs") 
kegg.sets.hs <- kegg.sets.hs[sigmet.idx.hs]

keggres <- gage(exprs = foldchanges,
                gsets = kegg.sets.hs,
                same.dir = T)

head(keggres$greater)

head(keggres$less)

library(pathview)

pathview(gene.data = df, 
         pathway.id = "03010",
         species = "hsa",
         out.suffix = "downregulated",
         kegg.native = T)

# The End
# Thank You
