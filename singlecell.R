
# ----------------------------------------------- 30.06.2025 --------------------------------------------- #

install.packages("Seurat")
install.packages("dplyr")
install.packages("patchwork")
install.packages("ggplot2")
install.packages("hdf5r") #bu paketi yükleyerek .h5 formatındaki 10X Genomics veri dosyasını okuyabilir hale geldin.
install.packages("BiocManager")
devtools::install_github("immunogenomics/presto")
install.packages("devtools")

#install sadece bir kere indirilir bilgisayara dosya yüklüyor
#library hafızadaki paketleri çağırma her projede vb kullanılır.

BiocManager::install("clusterProfiler")
BiocManager::install("limma")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("celldex")
BiocManager::install("SingleR")
BiocManager::install("scrapper")
BiocManager::install("presto")

#org.Hs.eg.db ->organism - homo sapiens - entire genome - data base
#org.Mm.eg.db -> Mus musculus (fare)


library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(hdf5r)
library(BiocManager)
library(clusterProfiler)
library(celldex)
library(org.Hs.eg.db)
library(SingleR)
library(scrapper)
library(presto)

# Load the PBMC dataset
pbmc.data <- Read10X_h5("/Users/elif/Desktop/Staj/10k_PBMC.h5") #PBMC (peripheral blood mononuclear cells) verisini okuduk.
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc10k", min.cells = 3, min.features = 200)
# min.cells= bir gen en az 3 hücrede görülsün görülmüyorsa sil. min.Features bir hücrede gen en az 200 kere görülsün görülmüyorsa sil. 
pbmc


# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
# benim bi objem var ona percent mt ekle dedik peki percent.mt nedir-> her satırın başında "^MT-" arıyor bu da mitokondriyel gen demek
#Eğer mitokondriyel gen belli bir seviyenin üstündeyse o hücre ya patlamış ya ölmüş demektir
#mt patterninini içeren yüzde görülmesini hesapla



# Show QC metrics for the first 5 cells in the control group
head(pbmc@meta.data, 5)



# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.0001)
#Bu şemadaki her nokta bir hücreye tekabül eder. İlk başta 
#nFeature ->içersinde kaç gen var 
# nCountRNA= ne kadar gen okumuşuz
# percent.mt = mitokondriyel genin yüzdesi
#Nasıl filtreleyeceğiz? ilk ikisini alttan yatay kesip üst tarafı alacağız. üçüncüyü alttan kesip alt tarafı alacağız.



# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt") + 
  theme(legend.position="none")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
  theme(legend.position="none")
plot1 + plot2



#Korelasyon: -1 ve 1 arasında değişiyor:
#0 'da hiç yok demek'a yakınsa bağlantı yok korelasyon yok
#1e yakında 1i artarken diğeride artar
#-1e yakında bir tanesi artarken diğeri azalır



pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15)
#pbmc verimin üstüne yaz: subset yap pbmc
#nFeature_RNA > 200 & nFeature_RNA < 5000 kısmı isteğe bağlı yazılıyor


VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.001)


pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
#Hiçbir şey yazmak zorunda değiliz içinde aslında içine zaten default)
#Bir sonrakinde verbose=True dene kodun nasıl çalıştığını gösteriyor.



pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc) + 
  theme(legend.position="top")
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE) + 
  theme(legend.position="none")
plot1 + plot2

#Her nokta = 1 Gen bu tabloda
#Average Expression = bizim verimizde hangi veriler(genler) daha çok görülmüş (sadece sayı olarak karşılaştırma olarak değil)
#Nokta ne kadar sağda ve yukarıda ise o kadar anlamlı 

all.genes <- rownames(pbmc)
#all genes siye bir obje atayım pbmc deki isimlerin hepsini içine koyduk.
pbmc <- ScaleData(pbmc, features = all.genes, verbose = FALSE)
#ölçeklendirme yaptık.

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc), verbose = FALSE)
#princible component analysis(PCA): Veriyi kaç tane hücre varsa o kadar boyut var gibi düşün. ^D evrende yaşadığımız için biz maks 3e indirgeyebiliyoruz.


DimPlot(pbmc, reduction = "pca")
#Bu aslında 10000 boyutlu bir olay ama bilgisayar bunu 2 ya da 3 boyuta ingirgiyor. 
#Her nokta bir hücre demek. Veri nasıl dağıldığını inceliyoruz.

ElbowPlot(pbmc)
#Princeble Components arttıkça gürültü sayısı artıyor. Yani veriye hata ekliyoruz. (Kolun dirseğini bulma)
#Bir veride kaç Princible component kullanılmalı -> Soru gelir

pbmc <- FindNeighbors(pbmc, dims = 1:10, verbose = FALSE)
pbmc <- FindClusters(pbmc, resolution = 0.5, verbose = FALSE)
pbmc <- RunUMAP(pbmc, dims = 1:10, verbose = FALSE)
DimPlot(pbmc, reduction = "umap")
#dims=1:10 dediysek her yerde 1:10 olmalı

#Hücreleri(veriyi) kümelere ayırıyor. KNN (Key Nearest Neighbour)
#Örnekte 13'e bölmüş. Kümedeki hücreler kendi arasında çok benzemeli yanındakilerle bir o kadar benzememeli.
#Bu kümelere ayırma sebebimiz hücre türü atamak istememiz. Bir kümenin Gen ekspresyonunun diğerlerinden farkına göre. DEG'lere bakıyoruz.
#Bir kümede diğerlerine kıyasla daha fazla sentezlenen genlere bakıyoruz.
#Marker genler var. Örneğin OCA1 geni bir kümede görünüyorsa bu küme b/nöron/kanser vb hücresidir.

pbmc <- RunTSNE(pbmc, dims = 1:10, verbose = FALSE)
DimPlot(pbmc, reduction = "tsne")

#eskiden tsne daha fazla kullanılmasına rağmen sonradan matematiksel olarak umap'ın daha doğru olduğu bulundu.


DimPlot(pbmc, reduction = "umap", label = TRUE)

# ----------------------------------------------- 01.07.2025 -------------------------------------------- # 

# Hücre tipi ataması

refHPCA <- HumanPrimaryCellAtlasData()
sce <- GetAssayData(object = pbmc, assay = "RNA", slot = "data")
prediction_HPCA_main <- SingleR(test = sce, ref = refHPCA, clusters = Idents(pbmc), labels = refHPCA$label.main)
prediction_HPCA_fine <- SingleR(test = sce, ref = refHPCA, clusters = Idents(pbmc), labels = refHPCA$label.fine)

pbmc$HPCA_main <- prediction_HPCA_main$labels[Idents(pbmc)]
pbmc$HPCA_fine <- prediction_HPCA_fine$labels[Idents(pbmc)]

DimPlot(pbmc, group.by = "HPCA_main", label = TRUE, pt.size = 0.5) + NoLegend()
# main = Genel türleri atama 
DimPlot(pbmc, group.by = "HPCA_fine", label = TRUE, pt.size = 0.5) + NoLegend()
# fine = alt türleri de atama yolu


# ----------------------------------------------- 02.07.2025 --------------------------------------------- #

cluster.all.markers <- FindAllMarkers(pbmc, logfc.threshold = 0.25, only.pos = TRUE, min.pct = 0.25)
degs <- cluster.all.markers$gene #ENSEMBLID Ver GeneID yerine



convert_to_ensembl <- function(ensg_ids) {
  ensembl_ids <- mapIds(org.Hs.eg.db, keys = ensg_ids, column = "ENSEMBL", keytype = "SYMBOL", multiVals = "first")
  return(ensembl_ids)
}

degs_ensembl <- convert_to_ensembl(degs)
# console kısmına degs_ensembl yazıyoruz




ego <- enrichGO(gene = degs, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
                ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05)


dotplot(ego) + ggtitle("GO Biological Processes") #hadi görselleştir
#Bu görselde hiçbir şey karşılaştırmadık aslında bu aşamada yapmamalıydık



convert_to_entrez <- function(ensg_ids) {
  entrez_ids <- mapIds(org.Hs.eg.db, keys = ensg_ids, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
  return(entrez_ids)
}
degs_entrez <- convert_to_entrez(degs)
degs_entrez <- degs_entrez[!is.na(degs_entrez)]
print(degs_entrez)
ekegg <- enrichKEGG(gene = degs_entrez, organism = "hsa", keyType = "kegg",
                    pvalueCutoff = 0.05, qvalueCutoff = 0.05)



#KEGG: kyoto encloypedia genes and genomes
#KEGG görselleştirdik:
if (!is.null(ekegg)) {
  dotplot(ekegg) + ggtitle("KEGG Descriptions")
} else {
  print("No Descriptions were enriched.")
}



clusters <- levels(Idents(pbmc))
deg_list <- list()
for (cluster in clusters) {
  deg_list[[cluster]] <- FindMarkers(pbmc, ident.1 = cluster, ident.2 = NULL,
                                     logfc.threshold = 0.25, min.pct = 0.1)
  write.csv(deg_list[[cluster]], paste0("DEGs_cluster_", cluster, ".csv"))
}

#ident1 = cluster, ident2= null ===== burada clusterdaki her kümeyi geriye kalan tüm kümelerle karşılaştırıyor
# ident1 = cluster1, ident2= cluster2 ===== 1.kümeyi 2.kümeyle kıyaslayıp farkına bakıyor



#ALTTAKİLERİN HEPSİ BERABER OKUNUYOR:
extract_significant_genes <- function(deg_results, p_val_thresh = 0.05, logfc_thresh = 0.25) {
  significant_genes <- deg_results[deg_results$p_val_adj < p_val_thresh & abs(deg_results$avg_log2FC) > logfc_thresh, ]
  gene_symbols <- rownames(significant_genes)
  entrez_ids <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  return(entrez_ids$ENTREZID)}

#ASIL İŞ:
#Sonuçlara GO ve KEGG uyguladık
perform_go_kegg_analysis <- function(gene_list, cluster_name) {
  go_results <- enrichGO(gene = gene_list, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05)
  kegg_results <- enrichKEGG(gene = gene_list, organism = "hsa", pAdjustMethod = "BH")
  write.csv(as.data.frame(go_results), paste0("GO_enrichment_cluster_", cluster_name, ".csv"))
  write.csv(as.data.frame(kegg_results), paste0("KEGG_enrichment_cluster_", cluster_name, ".csv"))
  return(list(go = go_results, kegg = kegg_results))}

#Enrichment -> Zenginleştirme

go_kegg_results <- list()

for (cluster in names(deg_list)) {
  significant_genes <- extract_significant_genes(deg_list[[cluster]])
  if (length(significant_genes) > 0) {
    go_kegg_results[[cluster]] <- perform_go_kegg_analysis(significant_genes, cluster)
  } else {
    cat("No significant genes found for cluster", cluster, "\n")
  } }  


#GÖRSELLEŞTİRME:
#https://metascape.org/gp/index.html#/main/step1

for (cluster in names(go_kegg_results)) {
  cat("Visualizing and saving plots for cluster", cluster, "\n")
  go_results <- go_kegg_results[[cluster]]$go
  kegg_results <- go_kegg_results[[cluster]]$kegg
  if (!is.null(go_results) && nrow(go_results) > 0) {
    p <- barplot(go_results, showCategory = 10, title = paste("GO Enrichment for Cluster", cluster))
    print(p)}}







