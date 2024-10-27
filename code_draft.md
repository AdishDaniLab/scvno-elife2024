

# code for analysis of single cell RNAseq of mouse vomeronasal neuroepithelium 
#### 10X sequencing requirements for 3' gene expression analysis V3.1 chemistry  

| Read |	Read 1 |	i7 Index |	i5 Index| 	Read 2|
-------|---------|-----------|----------|---------|
|Purpose| 	Cell barcode & UMI| 	Sample Index| 	N/A| 	Insert|
|Length |	28 |	8 |	0| 	91|

```shell
parallel --jobs 4 trim_galore --hardtrim5 28 -o ../trimmed_fastq ::: *_R1.fastq.gz
parallel --jobs 4 trim_galore --hardtrim5 91 -o ../trimmed_fastq ::: *_R2.fastq.gz
```

#### Alignment of trimmed reads to the mm10reference genome:
Mouse mm10 reference genome and gene annotation file(gtf) was downloaded from ensembl. Xntrpc was removed as it is a read through transcript and interferes with the alignemnt of Trcp2 reads. More information on read through transcripts can be found [here](https://www.ensembl.info/2019/02/11/annotating-readthrough-transcription-in-ensembl/)

We built\-4 a custom reference based on modified mm10 and used it for my analysis. As per 10x, alignment using cellranger need to be done separately for each sample in different well to avoid barcode clash

```shell
cellranger count --id=10x_M2 --transcriptome=/home/nandan/vno/scrna/Mus_musculus_10Xref_Xntrpc_removed --fastqs=/home/nandan/vno/scrna/Male_Female_rawdata/trimmed_fastq --sample=10x_M2
cellranger count --id=10x_F1 --transcriptome=/home/nandan/vno/scrna/Mus_musculus_10Xref_Xntrpc_removed --fastqs=/home/nandan/vno/scrna/Male_Female_rawdata/trimmed_fastq --sample=10x_F1
cellranger count --id=10x_F2 --transcriptome=/home/nandan/vno/scrna/Mus_musculus_10Xref_Xntrpc_removed --fastqs=/home/nandan/vno/scrna/Male_Female_rawdata/trimmed_fastq --sample=10x_F2
cellranger count --id=10x_M1 --transcriptome=/home/nandan/vno/scrna/Mus_musculus_10Xref_Xntrpc_removed --fastqs=/home/nandan/vno/scrna/Male_Female_rawdata/trimmed_fastq --sample=10x_M1
```


contents of female_aggr.scv
```shell
library_id,molecule_h5
10x_F1,/home/nandan/vno/scrna/male_female_run/10x_F1/outs/molecule_info.h5
10x_F2,/home/nandan/vno/scrna/male_female_run/10x_F2/outs/molecule_info.h5
```

contents of male_aggr.scv
```shell
library_id,molecule_h5
10x_M1,/home/nandan/vno/scrna/male_female_run/10x_M1/outs/molecule_info.h5
10x_M2,/home/nandan/vno/scrna/male_female_run/10x_M2/outs/molecule_info.h5
```
Aggregating male and female replicates by normalizing the depth of sequencing 
```shell
cellranger aggr --id=male_fenale --csv=all_aggr.csv --normalize=mapped
```
Contents of all_aggr.csv
```shell
library_id,molecule_h5
10x_F1,/home/nandan/vno/scrna/male_female_run/10x_F1/outs/molecule_info.h5
10x_F2,/home/nandan/vno/scrna/male_female_run/10x_F2/outs/molecule_info.h5
10x_M1,/home/nandan/vno/scrna/male_female_run/10x_M1/outs/molecule_info.h5
10x_M2,/home/nandan/vno/scrna/male_female_run/10x_M2/outs/molecule_info.h5
```
#seurat code for analysis:
```r
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(SoupX)

#Correcting expression matrix for ambient RNA contamination using SoupX 
male_female_soupX.input=load10X('/home/nandan/vno/scrna/male_female_run/male_female/outs/count/')
male_female_soupX.input=autoEstCont(male_female_soupX.input)
male_female_soupX.out=adjustCounts(male_female_soupX.input)
male_female_bc=CreateSeuratObject(male_female_soupX.out)
mito.genes <- grep(pattern = "^mt-", x = rownames(male_female_bc@assays[["RNA"]]), value = TRUE)
male_female_bc[["percent.mito"]] <- PercentageFeatureSet(male_female_bc, pattern = "^mt-")
male_female_bc_subset <- subset(x = male_female_bc, subset = nFeature_RNA > 200 & nFeature_RNA < 7000)
male_female_bc_subset <- NormalizeData(object = male_female_bc_subset, normalization.method = "LogNormalize", scale.factor = 10000)

male_female_bc_subset=FindVariableFeatures(object=male_female_bc_subset, mean.function=ExpMean, dispersion.function=LogVMR, x.low.cutoff =0.0125, x.high.cutoff=3, y.cutoff=0.5, nfeatures=2000)
male_female_bc_subset <- ScaleData(object = male_female_bc_subset, features=rownames(male_female_bc_subset@assays[["RNA"]]), vars.to.regress = c("nCount_RNA", "percent.mito"))
male_female_bc_subset <- RunPCA(object = male_female_bc_subset,  npcs = 50, verbose = FALSE)

male_female_bc_subset <- JackStraw(object = male_female_bc_subset, reduction = "pca", dims = 50, num.replicate = 100,  prop.freq = 0.1, verbose = FALSE)  
male_female_bc_subset <- ScoreJackStraw(object = male_female_bc_subset, dims = 1:50, reduction = "pca")
JackStrawPlot(object = male_female_bc_subset, dims = 1:50, reduction = "pca")
ElbowPlot(object = male_female_bc_subset) # used to determine the number of principal components
male_female_bc_subset <- FindNeighbors(male_female_bc_subset, reduction = "pca", dims = 1:37) #dims decided from ElbowPlot
male_female_bc_subset <- FindClusters(male_female_bc_subset, resolution = 0.3) #resolution decides number of clusters
male_female_bc_subset <- RunUMAP(object = male_female_bc_subset, dims = 1:50)
male_female_bc_subset.markers <- FindAllMarkers(object =male_female_bc_subset, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
```


Removing RBC contamination(based on Hbb-bs gene expression)
```r
hbbbs=names(which(male_female_bc_subset@assays[["RNA"]]$counts["Hbb-bs",]>0))
male_female_bc_subset_fil=subset(male_female_bc_subset,cells=hbbbs,invert=TRUE)
```
Removing dying cell clusters with mito genes as markers: cluster 12. 
Note: The cluster number may change based with run or machine. Run FindAllMarkers to identify the cluster enriched with mito-genes as markers
```r
male_female_bc_subset_fil_final=subset(male_female_bc_subset_fil,idents=12,invert=TRUE)
```


Renaming the seurat object and combining two Gnai2 clusters. Gnai2 cluster numbers a vary between and run and machine, identify the clusters using FindAllMarkers
```r
vno_scrna= male_female_bc_subset_fil_final
oldnames=levels(male_female_bc_subset_fil_final)
newnames=oldnames
newnames[1]="6" # renaming one of Gnai2 neurons clusters using other cluster ID
names(newnames) <- levels(vno_scrna)
vno_scrna <- RenameIdents(vno_scrna, newnames)
```

Solitiary chemosensory neurons and endothelial cells were manually labelled as clusters on UMAP using GUI 
```r
plot <- DimPlot(vno_scrna, reduction = "umap")
select.cells <- CellSelector(plot = plot)#cells were selected from GUI 
Idents(seurat_obj, cells = select.cells) <- "SubCells" #and subset based on these cells
Idents(vno_scrna, cells = select.cells) <- "SubCells"

select.cells <- CellSelector(plot = plot)#cells were selected from GUI 
Idents(seurat_obj, cells = select.cells) <- "SCC" #and subset based on these cells
Idents(vno_scrna, cells = select.cells) <- "SCC"
```

Adding male and female information to dataset (see aggr file to find the order of female and male)
```r
#getting counts of each male and female replicate
test=grep(pattern="1$", x=names(vno_scrna$nCount_RNA),value=TRUE)
length(test)
#[1] 2034
test=grep(pattern="2$", x=names(vno_scrna$nCount_RNA),value=TRUE)
length(test)
#[1] 2173
test=grep(pattern="3$", x=names(vno_scrna$nCount_RNA),value=TRUE)
length(test)
#[1] 2309
test=grep(pattern="4$", x=names(vno_scrna$nCount_RNA),value=TRUE)
length(test)
#[1] 2669

data = c(rep("female", 4207),rep("male",4978))
vno_scrna[["sex"]]=data
colors_dimplot=c("#D55E00","#917A00","#AB7100","#9240C7","#708200","#4B52CD","#C83189","#BB35B1","#0088C5","#00914C","#0091B1","#009472","#398A00","#C83658","#917A00","#AB7100","#9240C7","#708200")
cairo_ps("male_female_DimPlot.eps",height=5,width=10,pointsize=56)
DimPlot(vno_scrna,split.by="sex", cols=colors_dimplot)& NoLegend() + theme(panel.border = element_rect(color="black",linewidth=1.5))
dev.off()
saveRDS(vno_scrna, file = "vno_scrna.rds")
```


Clusterwise differential expression of male and female samples
```r

#scatter plot comapring average gene expression profile of cell clusters in male and female samples. Figure 1-figure supplement 1
vno_scrna[["celltype.sex"]] <- paste(Idents(vno_scrna), vno_scrna$sex, sep = "_")
test = AggregateExpression(vno_scrna, group.by = c("celltype.sex"), return.seurat = TRUE)
genes.to.label = c("Eif2s3y", "Ddx3y", "Uty", "Kdm5d")

clus.ids = names(table(Idents(vno_scrna)))
n = length(clus.ids)
for (i in clus.ids[3:n]) {
    i = 1
    i = i + 1
    cairo_ps(paste0(i, "_male vs female .eps"), height = 3, width = 3)
    # FeatureScatter(gnao1_neurons, feature1='Vmn2r20',
    # feature2='Vmn2r22',cols=rep('black',4),pt.size=2)&NoLegend()
    CellScatter(test, paste0("g", i, "-male"), paste0("g", i, "-female"), highlight = genes.to.label,
        pt.size = 1) & NoLegend()
    # LabelPoints(plot = p1, points = genes.to.label, repel = TRUE )
    dev.off()

}

#writing the source data to XLSX file - Figure 1-figure supplement 1-source data 1
library("data.table")
library(openxlsx)
options(java.parameters = "-Xmx10000m")
dummy = vno_scrna
Idents(dummy) <- "celltype.sex"

clus.ids = names(table(Idents(vno_scrna)))
n = length(clus.ids)
wb <- createWorkbook()
header_style <- createStyle(halign = "center", textDecoration = "bold")
for (i in clus.ids[3:n]) {
    t = FindMarkers(dummy, ident.1 = paste0(i, "_male"), ident.2 = paste0(i, "_female"))
    t = setorder(t, p_val_adj)
    t = cbind(Gene = row.names(t), t)

    print(i)
    addWorksheet(wb, paste0(i, " male vs female"))
    writeData(wb, paste0(i, " male vs female"), t, headerStyle = header_style)
    freezePane(wb, paste0(i, " male vs female"), firstRow = TRUE)
    setColWidths(wb, paste0(i, " male vs female"), cols = 1:ncol(t), widths = "auto")
}

saveWorkbook(wb, file = "Figure1-figure supplement 1- Source Data1.xlsx",
    overwrite = TRUE)

```
UMAP of cell types in vomeronasal neuroepithelium - Figure 1A
```r
colors_dimplot = c("#D55E00", "#917A00", "#AB7100", "#9240C7", "#708200", "#4B52CD",
    "#C83189", "#BB35B1", "#0088C5", "#00914C", "#0091B1", "#009472", "#398A00",
    "#C83658", "#917A00", "#AB7100", "#9240C7", "#708200")
cairo_ps("Figure1A.eps", height = 7, width = 7, pointsize = 56)
DimPlot(vno_scrna, cols = colors_dimplot) & NoLegend() & NoAxes()
dev.off()
```

Dot plot of clusters markers - Figure 1B 
```r
vno_scrna.markers = FindAllMarkers(vno_scrna)
top10 <- vno_scrna.markers %>% 
    group_by(cluster) %>% 
    top_n(10, avg_log2FC)

cairo_ps("Figure 1B.eps", height = 85, width = 20, pointsize = 24)
Seurat::DotPlot(vno_scrna, features = unique(top10$gene), dot.scale = 16, cols = "RdBu") +
    theme(axis.text.x = element_text(angle = 75, size = 38, vjust = 0.9, hjust = 1.1),
        axis.text.y = element_text(size = 38), legend.key.size = unit(1.5, "in"),
        legend.title = element_text(size = 38), legend.text = element_text(size = 38)) +
    coord_flip()
dev.off()
```
Generating feature plot of known neuronal and non-neuronal marker genes - Figure 1C
```r
cairo_ps("Figure 1C.eps", height = 5, width = 7, pointsize = 24)
FeaturePlot(male_female_bc_subset_fil_final, features = c("Gap43", "Neurod1", "Omp",
    "Gnao1", "Gnai2", "Cbr2"), ncol = 3) & scale_colour_gradientn(colours = (brewer.pal(n = 9,
    name = "RdPu")[2:9])) & NoLegend() & NoAxes()
dev.off()
```

Generating feature plot of markers for non-neuronal cell types - Figure 2A (top panel)
```r
cairo_ps("Figure2A.eps", height = 2.5, width = 17.5, pointsize = 24)
FeaturePlot(male_female_bc_subset_fil_final, features = c("H2-Aa", "C1qb", "Ly6d",
    "Ppic", "Plp1", "Mgp", "Trpm5"), ncol = 7) & scale_colour_gradientn(colours = (brewer.pal(n = 9,
    name = "RdPu")[2:9])) & NoLegend() & NoAxes()
dev.off()
```

Analysis of macrophages - Figure 2B and 2C
```r
# UMAP of macrophages
mf_vno = subset(vno_scrna, idents = c(7, 9, 11)) # cluster number may vary with run and machine parameters, identify the macrophage clusters based on the markers
mf_cols = c("#BB35B1", "#00914C", "#009472")
cairo_ps("Figure 2B.eps", height = 7, width = 3, pointsize = 56, family = "Arial")
DimPlot(mf_vno) & NoLegend() & NoAxes()
dev.off()

features_to_Plot = c("Tmem119", "Aif1", "P2ry12", "Cx3cr1", "Trem2", "Ccl3", "Ccl4",
    "Ccl7", "C1qa", "C1qb", "C1qc", "Il1a", "Il1b", "Adgre1", "Siglec1", "Cd68",
    "Napsa", "Lsp1", "H2-Aa")
cairo_ps("Figure 2C.eps", height = 3, width = 7, pointsize = 56, family = "Arial")
Seurat::DotPlot(mf_vno, features = features_to_Plot) & scale_colour_gradientn(colours = (brewer.pal(n = 9,
    name = "RdPu"))) & theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
mf_vno = FindVariableFeatures(object = mf_vno, mean.function = ExpMean, dispersion.function = LogVMR,
    x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, nfeatures = 2000)
mf_vno.markers = FindAllMarkers(mf_vno)
mf_vno.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 30) %>%
    ungroup() -> top30
write.table(top30, file = "Supplementary file 3-top30markers_macrophage_comparison.tsv", sep = "\t") 
```

Creating neurons object:
```r
vno_neurons=subset(vno_scrna,idents=c(2,4,6,10,12)) #cluster number may vary with run and machine parameters, identify the neuronal clusters based on the markers as described in the manuscript
vno_neurons=FindVariableFeatures(object=vno_neurons, mean.function=ExpMean, dispersion.function=LogVMR, x.low.cutoff =0.0125, x.high.cutoff=3, y.cutoff=0.5, nfeatures=2000)
vno_neurons <- ScaleData(object = vno_neurons, features=rownames(vno_neurons@assays[["RNA"]]), vars.to.regress = c("nCount_RNA", "percent.mito"))
vno_neurons <- RunPCA(object = vno_neurons,  npcs = 50, verbose = FALSE)
ElbowPlot(object = vno_neurons,ndims=70)
vno_neurons <- FindNeighbors(vno_neurons, reduction = "pca", dims = 1:37)
vno_neurons <- FindClusters(vno_neurons, resolution = 0.4) #resolution decides number of clusters, 0.4 gives 13 clusters
vno_neurons <- RunTSNE(object = vno_neurons, dims.use = 1:37, do.fast = TRUE)
vno_neurons <- RunUMAP(object = vno_neurons, dims = 1:50)
saveRDS(vno_neurons, file = "vno_neurons.rds")
```

Generating DimPlot/UMAP for neurons - Figure 3-figure supplement 1A
```r

cairo_ps("neurons_dimplot.eps",height=8,width=5,pointsize=48)
DimPlot(vno_neurons)
dev.off()
```

Feature plot for Ascl1, Gap43, Gnao1, Gnai2, Neurod1, Neurog1 and Omp - Figure 3-figure supplement 1B
```r
cairo_ps("Figure 3-figure supplement 1B.eps",height=5.5,width=14.5,pointsize=24)
FeaturePlot(vno_neurons, features = c("Ascl1","Neurod1", "Neurog1", "Gap43", "Gnao1", "Gnai2", "Omp"),
    ncol = 4) & scale_colour_gradientn(colours = (brewer.pal(n = 9, name = "RdPu")[2:9])) &
    NoLegend() & NoAxes()
dev.off()
```
Dotplot of neuronal cluster markers - Figure 3-figure supplement 1C
```r
vno_neurons.markers <- FindAllMarkers(vno_neurons)
# dotplot for top10 cluster markers
top10 <-vno_neurons.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC)
cairo_ps("Figure 3-figure supplement 1C.eps", height = 15, width = 70, pointsize = 24)
Seurat::DotPlot(vno_neurons, features = unique(top10$gene), dot.scale = 16,
    cols = "RdBu") + theme(axis.text.x = element_text(angle = 75, size = 38, vjust = 0.9,
    hjust = 1.1), axis.text.y = element_text(size = 38), legend.key.size = unit(1.5,
    "in"), legend.title = element_text(size = 38), legend.text = element_text(size = 38))
dev.off()
```
Comparing the effect of VRs on clustering 
```r
vno_neurons_noVR <- vno_neurons_new
tesy1=VariableFeatures(vno_neurons_noVR)
write.table(tesy1, file="var_features_temp.tsv")
VariableFeatures(vno_neurons_noVR) <- grep("Vmn", tesy1, invert="TRUE", value="TRUE")
vno_neurons_noVR <- ScaleData(object = vno_neurons_noVR, features=rownames(vno_neurons_noVR@assays[["RNA"]]), vars.to.regress = c("nCount_RNA", "percent.mito"))
vno_neurons_noVR <- RunPCA(object = vno_neurons_noVR,  npcs = 50, verbose = FALSE)
p=ElbowPlot(object = vno_neurons_noVR,ndims=50)
ggsave(p,file="test.png")
vno_neurons_noVR <- FindNeighbors(vno_neurons_noVR, reduction = "pca", dims = 1:37)
vno_neurons_noVR <- FindClusters(vno_neurons_noVR, resolution = 0.4) #resolution decides number of clusters, 0.4 gives 13 clusters
vno_neurons_noVR <- RunTSNE(object = vno_neurons_noVR, dims.use = 1:37, do.fast = TRUE)
vno_neurons_noVR <- RunUMAP(object = vno_neurons_noVR, dims = 1:50)
saveRDS(vno_neurons_noVR,file="vno_neurons_noVR.rds")
vno_neurons_noVR=readRDS("vno_neurons_noVR.rds")
anchors.1 <- FindTransferAnchors(
    reference = vno_neurons_new,
    query = vno_neurons_noVR,
    reference.reduction = 'pca',
    features = rownames(x = vno_neurons_new[["pca"]]@feature.loadings),
    dims = 1:37,
    nn.method = "annoy",
    k.filter = NA,
    verbose = TRUE
  )
vno_neurons_noVR.1 <- MapQuery(anchorset = anchors.1, reference = vno_neurons_new, query = vno_neurons_noVR,
    refdata = list(celltype = "seurat_clusters"), reference.reduction = "pca", reduction.model = "umap")

cairo_ps("Figure 3-figure supplement 2B.eps", height = 4, width = 4, pointsize = 24)
DimPlot(vno_neurons_noVR.1,reduction="ref.umap", label=T) & NoLegend() & NoAxes()
dev.off()

```


Pseduotime analysis on neurons object
```r
new_names=c("Gi","Gi","Gi","Go","Gi","Go","Go","7","8","9","Gi","Gi","Go")
names(new_names) <- levels(vno_neurons)
vno_neurons=RenameIdents(vno_neurons, new_names)
library(slingshot)
library(RColorBrewer)
set.seed(1)
neurons.sce=as.SingleCellExperiment(vno_neurons)
cl1=neurons.sce$ident
sds <- slingshot(neurons.sce, clusterLabels=cl1, start.clus=9, omega=FALSE,approx_points=300,reducedDim='UMAP')
cairo_ps("Figure 3A.eps", height = 12, width = 12, pointsize = 24)
plot(reducedDims(neurons.sce)$UMAP, col = brewer.pal(9,"Set1")[cl1], pch = 16, cex = 0.5)
lines(SlingshotDataSet(sds), lwd = 4, col = "black")
dev.off()

```


Differential expression of immature Gnao1 vs Gnai2 neurons - Figure 3B
```r

#Figure 3B and its source data
immatureGnao1vsGnai2=FindMarkers(neurons,ident.1=8,ident.2=7)
immatureGnao1vsGnai2=cbind(gene=row.names(immatureGnao1vsGnai2),immatureGnao1vsGnai2)  
library(openxlsx)
wb <- createWorkbook()
header_style <- createStyle(halign = "center", textDecoration = "bold")
addWorksheet(wb,"immatureGnao1vsGnai2")
writeData(wb,"immatureGnao1vsGnai2",immatureGnao1vsGnai2,headerStyle=header_style)
freezePane(wb,"immatureGnao1vsGnai2",firstRow = TRUE)
setColWidths(wb, "immatureGnao1vsGnai2",cols=1:ncol(immatureGnao1vsGnai2), widths="auto")
addWorksheet(wb,"filtered")
filtered.de=immatureGnao1vsGnai2 %>% filter((avg_log2FC > 1 | avg_log2FC < -1) & (p_val_adj < 0.0005))
writeData(wb,"filtered",filtered.de,headerStyle=header_style)
freezePane(wb,"filtered",firstRow = TRUE)
setColWidths(wb, "filtered",cols=1:ncol(filtered.de), widths="auto")
saveWorkbook(wb, file= "immature Gnao1 vs Gnai2 differential expression.xlsx", overwrite=TRUE)


immatureGnao1vsGnai2$gene=row.names(immatureGnao1vsGnai2)
cairo_ps("Figure 3B.eps",height=20,width=20,pointsize = 24)
EnhancedVolcano(immatureGnao1vsGnai2,
    lab = immatureGnao1vsGnai2$gene,
    x = 'avg_log2FC',
    y = 'p_val_adj', 
    selectLab = c('Gnao1','Gnai2','Robo2','Meis2','Tfap2e','Lmo4','Creb5','Prrxl1','Foxo1','Shisa8','Cnpy1','Gm36028'),
    title="Immature Gnao1 vs Gnai2",
    subtitle="Cut-off for log2FC is >|2| and P value is 10e-6",
    xlab = bquote(~Log[2]~ 'fold change'),
    pointSize = 9.0,
    labSize = 12.0,  
    colAlpha = 0.7,
    legendPosition = 'bottom',
    legendLabSize = 12,
    legendIconSize = 6.0,
    colConnectors = 'black',
    drawConnectors = TRUE,
    lengthConnectors= unit(0.0001, 'npc')
)

dev.off()
```
Feature plot of Transcription factors -Figure 3C
```r 
cairo_ps("Transcription factors_featureplot.eps", height = 10, width = 18, pointsize = 24)
FeaturePlot(neurons, features = c("Robo2", "Tfap2e", "Meis2", "Creb5", "Prrxl1",
    "Shisa8", "Lmo4", "Foxo1"), ncol = 4) & scale_colour_gradientn(colours = (brewer.pal(n = 9,
    name = "RdPu")[2:9]))
dev.off()
```

UMAP of Gnao1 neurons and Heatmap of V2R and H2-Mv expression - Figures 4A-4E
```r
gnao1_neurons <- subset(vno_neurons, idents=c(3,5,6,12))
cairo_ps("Figure 4A.eps", height = 4, width = 3.33, pointsize = 24)
DimPlot(gnao1_neurons,label=T) & NoLegend() & NoAxes() & xlim(-14,-7)
dev.off()

cairo_ps("Figure 4B-D.eps", height = 4, width = 10, pointsize = 24)
FeaturePlot(gnao1_neurons, features=c("Vmn2r1", "Vmn2r2", "H2-M10.3"), ncol =3) & 
    scale_colour_gradientn(colours = (brewer.pal(n = 9, name = "RdPu")[2:9])) &
    xlim(-14,-7) & NoAxes()
dev.off()

v2r_h2mv_family_names = read.table("V2r_h2-mv_genes_families.txt", sep = "\t", header = T)
genes.to.plot = v2r_h2mv_family_names$receptor
cairo_ps("Figure 4E.eps", height = 40, width = 50, pointsize = 20)
DoHeatmap(gnao1_neurons, features = genes.to.plot) + scale_fill_gradientn(colours = rev(brewer.pal(n = 9,
    name = "RdBu"))) + theme(axis.text.y = element_text(size = 24), legend.key.size = unit(1.5,
    "in"), legend.title = element_text(size = 42), legend.text = element_text(size = 42))
dev.off()
```


Coexpression analysis
```r
# function to find co-expression counts
get_coexpression_combis <- function(expression_table, threshold, sufix) {
    temp = expression_table

    temp[temp < threshold] = 0  #thres of expression is 1.25
    combi = vector()
    no_genes = vector()
    for (i in colnames(temp)) {
        a = rownames(temp)[temp[, i] > 0]
        # name=a
        name = paste(unlist(a), collapse = ",")
        no_genes[i] = length(a)
        combi[i] = name
    }

    combi_counts = table(combi)
    total_exp = colSums(temp)

    gene_no = paste("total", sufix, "genes", sep = "_")
    combination = paste(sufix, "combinations", sep = "_")
    confi = paste(sufix, "confidence", sep = "_")
    col2 = paste(sufix, "total_expression", sep = "_")
    toreturn = data.frame(gene_no = unlist(no_genes), tot = unlist(total_exp), combination = unlist(combi))
    toreturn$confi <- ifelse(combi_counts[toreturn$combination] >= 2, "High confidence",
        ifelse(combi_counts[toreturn$combination] == 1, "Low confidence", "undefined"))
    # toreturn$confi <- ifelse( combi_counts[toreturn$combination]>= 2, 'High
    # confidence','Low confidence') ifelse(combi_counts[toreturn$combination]
    # == 2, 'Medium confidence', 'Low confidence'))

    colnames(toreturn) = c(gene_no, col2, combination, confi)
    return(toreturn)

}

#collecting names of all detected VRS
vmn1r = grep(pattern = "^Vmn1r", x = rownames(vno_neurons@assays[["RNA"]]),
    value = TRUE)  # gene names of all V1Rs
vmn2r = grep(pattern = "^Vmn2r", x = rownames(vno_neurons@assays[["RNA"]]),
    value = TRUE)  # gene names of all V2Rs

#categorizing H2-Mvs and V2R as per their family    
famc_vmn2r = c("Vmn2r1", "Vmn2r2", "Vmn2r3", "Vmn2r4", "Vmn2r5", "Vmn2r6", "Vmn2r7")
h2mv_genes = c("H2-M10.2", "H2-M10.1", "H2-M10.3", "H2-M10.4", "H2-M11", "H2-M9",
    "H2-M1", "H2-M10.5", "H2-M10.6", "H2-M5", "H2-M3", "H2-M2")
abd_vmn2r = vmn2r[!vmn2r %in% famc_vmn2r]
abd_h2mv = c(abd_vmn2r, h2mv_genes)


v2r_exp_gnao1_clus = data.frame(gnao1_neurons[["RNA"]]$data[vmn2r, ])
famc_v2r_exp_gnao1_clus = data.frame(gnao1_neurons[["RNA"]]$data[famc_vmn2r, ])
abd_v2r_exp_gnao1_clus = data.frame(gnao1_neurons[["RNA"]]$data[abd_vmn2r, ])
abd_v2r_h2mv_exp_gnao1_clus = data.frame(gnao1_neurons[["RNA"]]$data[abd_h2mv, ])
v2r_h2mv_exp_gnao1_clus = data.frame(gnao1_neurons[["RNA"]]$data[vmn2r_h2mv, ])
v1r_exp_gnai2_clus = data.frame(gnai2_neurons[["RNA"]]$data[vmn1r, ])
h2mv_exp_gnao1_clus = data.frame(gnao1_neurons[["RNA"]]$data[h2mv_genes, ])

write.table(data.frame(rowSums(h2mv_exp_gnao1_clus > 1.25)), file = "NumberOfCellsExpressingEach_ H2mv.tsv",
    sep = "\t")
write.table(data.frame(rowSums(v2r_h2mv_exp_gnao1_clus > 1.25)), file = "NumberOfCellsExpressingEach_ V2R.tsv",
    sep = "\t")
write.table(data.frame(rowSums(h2mv_exp_gnao1_clus > 1.25)), file = "NumberOfCellsExpressingEach_ H2mv.tsv",
    sep = "\t")


abd_table = get_coexpression_combis(abd_v2r_exp_gnao1_clus, 1.25, "ABD")
famc_table = get_coexpression_combis(famc_v2r_exp_gnao1_clus, 1.25, "famc")
h2mv_table = get_coexpression_combis(h2mv_exp_gnao1_clus, 1.25, "h2mv")
all_v2r_table = get_coexpression_combis(v2r_exp_gnao1_clus, 1.25, "all_v2r")
v1r_table = get_coexpression_combis(v1r_exp_gnai2_clus, 2.5, "v1r")

total_gnao1_table = cbind(all_v2r_table, famc_table, abd_table, h2mv_table)

total_gnao1_table$famc_vmn2r = ifelse(grepl("Vmn2r1", total_gnao1_table$famc_combinations),
    "Vmn2r1", ifelse(grepl("Vmn2r2", total_gnao1_table$famc_combinations), "Vmn2r2",
        "Vmn2r3-Vmn2r7"))


total_gnao1_table$famc_type = ifelse(grepl("Vmn2r1", total_gnao1_table$famc_combinations),
    "C1", "C2")
total_gnao1_table <- total_gnao1_table %>%
    replace(is.na(.), "undefined")

total_gnao1_table <- total_gnao1_table %>%
    mutate(famc_vmn2r = fct_relevel(famc_vmn2r, "Vmn2r1", "Vmn2r2", "Vmn2r3-Vmn2r7"))

total_gnao1_table <- total_gnao1_table %>%
    mutate(ABD_confidence = fct_relevel(ABD_confidence, "High confidence", "Low confidence",
        "undefined"))

total_gnao1_table <- total_gnao1_table %>%
    mutate(famc_confidence = fct_relevel(famc_confidence, "High confidence", "Low confidence",
        "undefined"))

total_gnao1_table <- total_gnao1_table %>%
    mutate(h2mv_confidence = fct_relevel(h2mv_confidence, "High confidence", "Low confidence",
        "undefined"))

total_gnao1_table <- total_gnao1_table %>%
    mutate(all_v2r_confidence = fct_relevel(all_v2r_confidence, "High confidence",
        "Low confidence", "undefined"))
v1r_table <- v1r_table %>%
    replace(is.na(.), "undefined")

v1r_table <- v1r_table %>%
    mutate(v1r_confidence = fct_relevel(v1r_confidence, "High confidence", "Low confidence",
        "undefined"))



gnao1_gene_gnao1_clus = data.frame(gnao1_neurons[["RNA"]]$data["Gnao1", ])
colnames(gnao1_gene_gnao1_clus) = "gnao1_expression"

gnai2_gene_gnai2_clus = data.frame(gnai2_neurons[["RNA"]]$data["Gnai2", ])
colnames(gnai2_gene_gnai2_clus) = "gnai2_expression"

#Figure 4F
cairo_ps("famc_coexpression_frequency.eps", height = 5, width = 3.75)
ggplot(df <- filter(total_gnao1_table, famc_confidence == "undefined" | famc_confidence ==
    "High confidence"), aes(x = factor(total_famc_genes), fill = factor(famc_type),
    group = factor(famc_type))) + geom_bar(width = 0.75) + scale_fill_manual(values = cbPalette) +
    # facet_wrap(~famc_confidence,nrow = 1)+
geom_text(aes(label = after_stat(count)), stat = "count", position = position_stack(0.5)) +
    theme_classic() + theme(text = element_text(size = 12), legend.position = "none",
    axis.title = element_text(face = "bold", size = 12)) + xlab("Number of family-C receptors\n expressed in a cell") +
    ylab("Number of cells")
# ggsave(p,file='famc.png',height=7,width=7)
dev.off()



#Figure 4G
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00",
    "#CC79A7")
cairo_ps("ABD_coexpression_frequency.eps", height = 5, width = 3.75)
ggplot(df <- filter(total_gnao1_table, ABD_confidence == "High confidence" | ABD_confidence ==
    "undefined"), aes(x = (total_ABD_genes), fill = (famc_type), group = factor(famc_type))) +
    geom_bar(width = 0.75) + coord_cartesian(xlim = c(0, 5)) + scale_fill_manual(values = cbPalette) +
    # facet_wrap(~ABD_confidence,nrow = 1) +
geom_text(aes(label = after_stat(count)), position = position_stack(0.5), stat = "count") +
    theme_classic() + theme(text = element_text(size = 12), legend.position = "none",
    axis.title = element_text(face = "bold", size = 12)) + xlab("No. of family-ABD receptors\nexpressed in a cell") +
    ylab("Number of cells")  #+xlim(-0.1,5)
dev.off()

#Figure 4H
cairo_ps("h2mv_coexpression_frequency.eps", height = 5, width = 3.75)
ggplot(df <- filter(total_gnao1_table, h2mv_confidence == "undefined" | h2mv_confidence ==
    "High confidence"), aes(x = factor(total_h2mv_genes), fill = factor(famc_type),
    group = factor(famc_type))) + geom_bar(width = 0.75) + scale_fill_manual(values = cbPalette) +
    # facet_wrap(~h2mv_confidence,nrow = 1)+
geom_text(aes(label = after_stat(count)), stat = "count", position = position_stack(0.5)) +
    theme_classic() + theme(text = element_text(size = 12), legend.position = "none",
    axis.title = element_text(face = "bold", size = 12)) + xlab("Number of H2-Mv genes\n expressed in a cell") +
    ylab("Number of cells")
# ggsave(p,file='h2mv.png',height=7,width=7)
dev.off()

#Figure 4I 

cairo_ps("famc_levels.eps", height = 5, width = 4)
t1 = data.frame(ID = rownames(gnao1_gene_gnao1_clus), total_famc_genes = total_gnao1_table$total_famc_genes,
    total_expression = total_gnao1_table$famc_total_expression, type = total_gnao1_table$famc_type)
t1$total_expression[t1$total_famc_genes >= 2 & t1$type == "C1"] = NA  #removing Vmn2r1+Vmn2r7 expressing cells from the plot
t1 = na.omit(t1)
t2 = data.frame(ID = rownames(gnao1_gene_gnao1_clus), total_famc_genes = total_gnao1_table$total_famc_genes,
    total_expression = gnao1_gene_gnao1_clus$gnao1_expression, type = rep("Gnao1",
        1018))

t = rbind(t1, t2)
t <- t %>%
    mutate(type = fct_relevel(type, "C1", "C2"))
p = ggplot(t, aes(x = factor(total_famc_genes), y = total_expression, color = factor(type))) +
    coord_cartesian(xlim = c(2, 6)) + geom_boxplot() + geom_jitter(position = position_jitterdodge(0.3),
    alpha = 0.2) + theme_classic() + scale_colour_manual(values = cbPalette) + theme(text = element_text(size = 12),
    legend.position = "none", axis.title = element_text(face = "bold", size = 12)) +
    xlab("Number of family-C V2Rs\n expressed in a cell") + ylab("Total family-C V2R expression")
ggsave(p, file = "C1_C2.png")
dev.off()



cairo_ps("abd_levels.eps", height = 5, width = 4)
t1 = data.frame(ID = rownames(gnao1_gene_gnao1_clus), total_ABD_genes = total_gnao1_table$total_ABD_genes,
    total_expression = total_gnao1_table$ABD_total_expression, type = total_gnao1_table$famc_type)

t2 = data.frame(ID = rownames(gnao1_gene_gnao1_clus), total_ABD_genes = total_gnao1_table$total_ABD_genes,
    total_expression = gnao1_gene_gnao1_clus$gnao1_expression, type = rep("Gnao1",
        1018))
t3 = data.frame(ID = rownames(creld2_gene_gnao1_clus), total_ABD_genes = total_gnao1_table$total_ABD_genes,
    total_expression = creld2_gene_gnao1_clus$creld2_expression, type = rep("Creld2",
        1018))
t = rbind(t1, t2)
p = ggplot(t, aes(x = factor(total_ABD_genes), y = total_expression, color = type)) +
    coord_cartesian(xlim = c(2, 3), ylim = c(0, 12)) + geom_boxplot() + geom_jitter(position = position_jitterdodge(0.3),
    alpha = 0.2) + theme_classic() + scale_colour_manual(values = cbPalette) + theme(text = element_text(size = 12),
    legend.position = "none", axis.title = element_text(face = "bold", size = 12)) +
    xlab("Number of family-ABD V2Rs\n expressed in a cell") + ylab("Total family-ABD V2R expression")

dev.off()


cairo_ps("h2mv_levels.eps", height = 5, width = 4)
t1 = data.frame(ID = rownames(gnao1_gene_gnao1_clus), total_h2mv_genes = total_gnao1_table$total_h2mv_genes,
    total_expression = total_gnao1_table$h2mv_total_expression, type = total_gnao1_table$famc_type)
t1$total_expression[t1$total_h2mv_genes >= 2 & t1$type == "C1"] = NA
t1 = na.omit(t1)
t2 = data.frame(ID = rownames(gnao1_gene_gnao1_clus), total_h2mv_genes = total_gnao1_table$total_h2mv_genes,
    total_expression = gnao1_gene_gnao1_clus$gnao1_expression, type = rep("Gnao1",
        1018))
t = rbind(t1, t2)
t <- t %>%
    mutate(type = fct_relevel(type, "C1", "C2", "Gnao1"))
ggplot(t, aes(x = factor(total_h2mv_genes), y = total_expression, color = factor(type))) +
    coord_cartesian(xlim = c(1, 6)) + geom_boxplot() + geom_jitter(position = position_jitterdodge(0.3),
    alpha = 0.2) + theme_classic() + scale_colour_manual(values = cbPalette) + theme(text = element_text(size = 12),
    legend.position = "none", axis.title = element_text(face = "bold", size = 12)) +
    xlab("Number of H2-Mv genes\n expressed in a cell") + ylab("Total H2-Mv expression")
# ggsave(p,file='C1_C2.png')
dev.off()
```


Figure 4-Figure supplement 1
```r
# plotting frequency of V1R, V2R and H2-Mv - Figure 4-figure supplement 1A-1C

v2r_freq = data.frame(V2R = names(sort(rowSums(v2r_exp_gnao1_clus > 1.25), decreasing = T)),
    freq = sort(rowSums(v2r_exp_gnao1_clus > 1.25), decreasing = T))
v2r_freq$V2R <- factor(v2r_freq$V2R, levels = v2r_freq$V2R)
cairo_ps("Figure 4-figure supplement 1B.eps", height = 7, width = 7)
ggplot(v2r_freq[1:20, ], aes(x = V2R, y = freq)) + geom_bar(stat = "identity") +
    theme_classic() + theme(axis.text = element_text(size = 18), axis.text.x = element_text(angle = 75,
    hjust = 1))
dev.off()

v1r_freq = data.frame(V1R = names(sort(rowSums(v1r_exp_gnai2_clus > 1.25), decreasing = T)),
    freq = sort(rowSums(v1r_exp_gnai2_clus > 1.25), decreasing = T))
v1r_freq$V1R <- factor(v1r_freq$V1R, levels = v1r_freq$V1R)
cairo_ps("Figure 4-figure supplement 1A.eps", height = 7, width = 7)
ggplot(v1r_freq[1:20, ], aes(x = V1R, y = freq)) + geom_bar(stat = "identity") +
    theme_classic() + theme(axis.text = element_text(size = 18), axis.text.x = element_text(angle = 75,
    hjust = 1))
dev.off()

h2mv_freq = data.frame(h2mv = names(sort(rowSums(h2mv_exp_gnao1_clus > 1.25), decreasing = T)),
    freq = sort(rowSums(h2mv_exp_gnao1_clus > 1.25), decreasing = T))
h2mv_freq$h2mv <- factor(h2mv_freq$h2mv, levels = h2mv_freq$h2mv)
cairo_ps("Figure 4-figure supplement 1C.eps", height = 7, width = 7)
ggplot(h2mv_freq, aes(x = h2mv, y = freq)) + geom_bar(stat = "identity") + theme_classic() +
    theme(axis.text = element_text(size = 18), axis.text.x = element_text(angle = 75,
        hjust = 1))
dev.off()


test = as.vector(unlist(v1r_exp_gnai2_clus))
nonzero_exp_values = data.frame(exp = test[which(test != 0)])
cairo_ps("v1r_expression_threshold.eps", height = 7, width = 7)
ggplot(nonzero_exp_values, aes(x = exp)) + geom_density(lwd = 1.5) + xlim(0, 10) +
    geom_vline(xintercept = 2.5, colour = "#BB0000", linetype = "dashed") + theme_clean() +
    theme(axis.text = element_text(size = 32))
dev.off()

test = as.vector(unlist(v2r_exp_gnao1_clus))
nonzero_exp_values = data.frame(exp = test[which(test != 0)])
cairo_ps("v2r_expression_threshold.eps", height = 7, width = 7)
ggplot(nonzero_exp_values, aes(x = exp)) + geom_density(lwd = 1.5) + xlim(0, 10) +
    geom_vline(xintercept = 1.25, colour = "#BB0000", linetype = "dashed") + theme_clean() +
    theme(axis.text = element_text(size = 32))
dev.off()

test = as.vector(unlist(h2mv_exp_gnao1_clus))
nonzero_exp_values = data.frame(exp = test[which(test != 0)])
cairo_ps("v2r_expression_threshold.eps", height = 7, width = 7)
ggplot(nonzero_exp_values, aes(x = exp)) + geom_density(lwd = 1.5) + xlim(0, 10) +
    geom_vline(xintercept = 1.25, colour = "#BB0000", linetype = "dashed") + theme_clean() +
    theme(axis.text = element_text(size = 32))
dev.off()



```

Figure 4-figure supplement 2
```r
cairo_ps("H2-M1_M9_M11_featureplot_new.eps", height = 4, width = 10, pointsize = 24)
FeaturePlot(gnao1_neurons, features = c("H2-M1", "H2-M9", "H2-M11"), ncol = 3) &
    scale_colour_gradientn(colours = (brewer.pal(n = 9, name = "RdPu")[3:9])) & NoAxes() &
    xlim(-14, -7)
dev.off()

cairo_ps("H2-M10_genes_featureplot_new.eps", height = 4, width = 20, pointsize = 24)
FeaturePlot(gnao1_neurons, features = c("H2-M10.1", "H2-M10.2", "H2-M10.3", "H2-M10.4",
    "H2-M10.5", "H2-M10.6"), ncol = 6) & scale_colour_gradientn(colours = (brewer.pal(n = 9,
    name = "RdPu")[3:9])) & NoAxes() & xlim(-14, -7)
dev.off()

cairo_ps("Vmn2r81_Vmn2r82_featureplot_new.eps", height = 4, width = 6.7, pointsize = 24)
FeaturePlot(gnao1_neurons, features = c("Vmn2r81", "Vmn2r82"), ncol = 2) & scale_colour_gradientn(colours = (brewer.pal(n = 9,
    name = "RdPu")[3:9])) & NoAxes() & xlim(-14, -7)
dev.off()


cairo_ps("H2-M9_Vmn2r1_scatter.eps", height = 2.2, width = 2.2)
FeatureScatter(gnao1_neurons, feature1 = "H2-M9", feature2 = "Vmn2r1", cols = rep("black",
    4), pt.size = 2) & NoLegend()
dev.off()

cairo_ps("H2-M9_Vmn2r2_scatter.eps", height = 2.2, width = 2.2)
FeatureScatter(gnao1_neurons, feature1 = "H2-M9", feature2 = "Vmn2r2", cols = rep("black",
    4), pt.size = 2) & NoLegend()
dev.off()

cairo_ps("H2-M1_Vmn2r1_scatter.eps", height = 2.2, width = 2.2)
FeatureScatter(gnao1_neurons, feature1 = "H2-M1", feature2 = "Vmn2r1", cols = rep("black",
    4), pt.size = 2) & NoLegend()
dev.off()

cairo_ps("H2-M1_Vmn2r2_scatter.eps", height = 2.2, width = 2.2)
FeatureScatter(gnao1_neurons, feature1 = "H2-M1", feature2 = "Vmn2r2", cols = rep("black",
    4), pt.size = 2) & NoLegend()
dev.off()

cairo_ps("H2-M11_Vmn2r1_scatter.eps", height = 2.2, width = 2.2)
FeatureScatter(gnao1_neurons, feature1 = "H2-M11", feature2 = "Vmn2r1", cols = rep("black",
    4), pt.size = 2) & NoLegend()
dev.off()

cairo_ps("H2-M11_Vmn2r2_scatter.eps", height = 2.2, width = 2.2)
FeatureScatter(gnao1_neurons, feature1 = "H2-M11", feature2 = "Vmn2r2", cols = rep("black",
    4), pt.size = 2) & NoLegend()
dev.off()
```

Figure 4- Figure supplement 3 - refer hills_etal_data_analysis.r

Figure 4- figure supplement 4: V1R coexpression 
```r

selected_v1rs = c("Vmn1r85", "Vmn1r86", "Vmn1r184", "Vmn1r185", "Vmn1r56", "Vmn1r57",
    "Vmn1r168", "Vmn1r177", "Vmn1r61", "Vmn1r62", "Vmn1r37", "Vmn1r38", "Vmn1r37",
    "Vmn1r236", "Vmn1r35", "Vmn1r37")
# generating heatmap for selected V1R expression in Gnai2 neurons
cairo_ps("Figure 4-figure supplement 4A .eps", height = 20, width = 45)
DoHeatmap(gnai2_neurons, features = selected_v1rs, label = FALSE, draw.lines = FALSE,
    group.bar = FALSE) + scale_fill_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")),
    limits = c(-3, 3)) + theme(axis.text.y = element_text(size = 40), legend.key.size = unit(1.5,
    "in"), legend.title = element_text(size = 42), legend.text = element_text(size = 42))
dev.off()

cairo_ps("Figure 4-figure supplement 4B.eps", height = 5, width = 3.75)
ggplot(df <- filter(v1r_table, v1r_confidence == "undefined" | v1r_confidence ==
    "High confidence"), aes(x = factor(total_v1r_genes), fill = "#CC79A7")) + geom_bar(width = 0.75,
    aes(fill = "#CC79A7")) + geom_text(aes(label = after_stat(count)), stat = "count",
    position = position_stack(0.5)) + theme_classic() + theme(text = element_text(size = 12),
    legend.position = "none", axis.title = element_text(face = "bold", size = 12)) +
    xlab("Number of V1R genes\n expressed in a cell") + ylab("Number of cells")
# ggsave(p,file='v1r.png',height=7,width=7)
dev.off()

cairo_ps("vFigure 4-figure supplement 4C.eps", height = 5, width = 4)
t1 = data.frame(ID = rownames(gnai2_gene_gnai2_clus), total_v1r_genes = v1r_table$total_v1r_genes,
    total_expression = v1r_table$v1r_total_expression, type = rep("V1R", 3647))
t2 = data.frame(ID = rownames(gnai2_gene_gnai2_clus), total_v1r_genes = v1r_table$total_v1r_genes,
    total_expression = gnai2_gene_gnai2_clus$gnai2_expression, type = rep("Gnai2",
        3348))
t = rbind(t1, t2)
t <- t %>%
    mutate(type = fct_relevel(type, "V1R", "Gnai2"))

ggplot(t, aes(x = factor(total_v1r_genes), y = total_expression, color = type)) +
    coord_cartesian(xlim = c(1, 3), ylim = c(0, 12)) + geom_boxplot() + geom_jitter(position = position_jitterdodge(0.3),
    alpha = 0.2) + theme_classic() + scale_colour_manual(values = c("#CC79A7", "#009E73")) +
    theme(text = element_text(size = 12), legend.position = "none", axis.title = element_text(face = "bold",
        size = 12)) + xlab("Number of V1R genes\n expressed in a cell") + ylab("Total V1R expression")

dev.off()
```




Figure 5A Differential expression of mature Gnao1 vs Gnai2 neurons
```r
gnao1_vs_gnai2.de= FindMarkers(vno_scrna, ident.1="2", ident.2="6")
gnao1_vs_gnai2.de=cbind(gene=row.names(gnao1_vs_gnai2.de),gnao1_vs_gnai2.de)  
library(openxlsx)
wb <- createWorkbook()
header_style <- createStyle(halign = "center", textDecoration = "bold")
addWorksheet(wb,"Gnao1_vs_Gnai2")
writeData(wb,"Gnao1_vs_Gnai2",gnao1_vs_gnai2.de,headerStyle=header_style)
freezePane(wb,"Gnao1_vs_Gnai2",firstRow = TRUE)
setColWidths(wb, "Gnao1_vs_Gnai2",cols=1:ncol(gnao1_vs_gnai2.de), widths="auto")
addWorksheet(wb,"filtered")
filtered.de=gnao1_vs_gnai2.de %>% filter((avg_log2FC > 1 | avg_log2FC < -1) & (p_val_adj < 0.0005))
writeData(wb,"filtered",filtered.de,headerStyle=header_style)
freezePane(wb,"filtered",firstRow = TRUE)
setColWidths(wb, "filtered",cols=1:ncol(filtered.de), widths="auto")
saveWorkbook(wb, file= "mature Gnao1 vs Gnai2 differential expression.xlsx", overwrite=TRUE)


library(EnhancedVolcano)
cairo_ps("volcanoplot_matureneurons_large.eps", height = 40, width = 40, pointsize = 24)
EnhancedVolcano(gnao1_vs_gnai2.de, lab = gnao1_vs_gnai2.de$gene, x = "avg_log2FC",
    y = "p_val_adj", title = "Gnao1 vs Gnai2", subtitle = "Cut-off for log2FC is >|2| and P value is 10e-6",
    xlab = bquote(~Log[2] ~ "fold change"), pointSize = 4, labSize = 6, colAlpha = 0.7,
    legendPosition = "bottom", legendLabSize = 12, legendIconSize = 6, colConnectors = "black",
    drawConnectors = TRUE, lengthConnectors = unit(0.01, "npc"))
dev.off()
selectedlabels <- gnao1_vs_gnai2.de %>% 
                filter(abs(avg_log2FC) >= 4 & p_val_adj <= 1e-50)
```


Figure 6A Gene ontology of Gnao1 enriched genes 
```r
library(clusterProfiler)
library(msigdbr)
h_genesets = msigdbr(species = "mouse", category = "H")
go_bp_genesets = msigdbr(species = "mouse", category = "C5", subcategory = "GO:BP")
go_mf_genesets = msigdbr(species = "mouse", category = "C5", subcategory = "GO:MF")
kegg_genesets = msigdbr(species = "mouse", category = "C2", subcategory = "CP:KEGG")
reactome_genesets = msigdbr(species = "mouse", category = "C2", subcategory = "CP:REACTOME")
h_geneset_t2g = h_genesets %>%
    dplyr::distinct(gs_name, gene_symbol) %>%
    as.data.frame()
go_bp_geneset_t2g = go_bp_genesets %>%
    dplyr::distinct(gs_name, gene_symbol) %>%
    as.data.frame()
go_mf_geneset_t2g = go_mf_genesets %>%
    dplyr::distinct(gs_name, gene_symbol) %>%
    as.data.frame()
kegg_geneset_t2g = kegg_genesets %>%
    dplyr::distinct(gs_name, gene_symbol) %>%
    as.data.frame()
reactome_geneset_t2g = reactome_genesets %>%
    dplyr::distinct(gs_name, gene_symbol) %>%
    as.data.frame()
gnao1vsgnai2 = dplyr::select(gnao1_vs_gnai2.de, gene, p_val_adj, avg_log2FC)
geneList = gnao1vsgnai2[gnao1vsgnai2$p_val_adj < 0.05, ]$avg_log2FC
names(geneList) = gnao1vsgnai2[gnao1vsgnai2$p_val_adj < 0.05, ]$gene
geneList = sort(geneList, decreasing = TRUE)
em <- enricher(names(geneList[geneList >= 1]), TERM2GENE = go_bp_geneset_t2g)
cairo_ps("Figure 6A.eps", height = 10, width = 10, pointsize = 24)
dotplot(em, showCategory = 20, font.size = 12, x = "GeneRatio", color = "p.adjust")
dev.off()
```

Representing RNA levels Gnao1 enriched ER genes tested by IHC using violin plot -Figure 7-figure supplement 2
```r
vno_scrna = readRDS("vno_scrna.rds")
p = VlnPlot(vno_scrna, features = "Hspa5", idents = c(2, 6))
cairo_ps("Hspa5_violin_Gnao1_Gnai2.eps", height = 4, width = 3, pointsize = 24)
q=ggplot(p$data, aes(x = ident, y = Hspa5, fill = ident)) + geom_violin(trim = FALSE,
    linewidth = 1) + geom_jitter(aes(colour = ident), position = position_jitter(0.15),
    shape = 16, alpha = 0.2) + #scale_color_brewer(palette='Dark2') + geom_boxplot(width=0.05, fill='white')+ shape
    labs(title = "Hspa5", y = "Expression level") + scale_x_discrete(labels = c("Gnai2\nneurons",
    "Gnao1\nneurons")) + theme_classic()
dev.off()
#repeat by changing gene name for other ER genes in the 

#pvalue and Log2FC is obtained by FindMarkers fuction.
FindMarkers(vno_scrna, features = c("Sec61b", "Hspa5", "Atl1", "Canx", "Rtn4", "Hsp90b1",
    "P4hb", "Gnao1", "Gnai2", "Reep5", "Ckap4"), ident.1 = "2", ident.2 = "6", min.pct = 5e-04,
    logfc.threshold = 1e-04)

```



Featureplot and heatmap of Vmn2r1 and ER genes in neurons object - Figure 9A and 9B
```r

cairo_ps("Figure 9A.eps", height = 7, width = 7, pointsize = 20)
FeaturePlot(neurons, features = c("Gnao1", "Vmn2r1", "Sdf2l1", "Manf"), ncol = 2) &
    scale_colour_gradientn(colours = (brewer.pal(n = 9, name = "RdPu")[2:9])) & NoAxes()
ggsave(p, file = "test.png", height = 3, width = 5, units = "in")
dev.off()

vsns = vno_neurons
# head(Idents(vsns))

newnames_neurons <- c("n11","n13","n8","n1","n12","n3","n2","n7","n6","n5","n10","n9","n4")
names(newnames_neurons)<- levels(vsns)
vsns <- RenameIdents(vsns, newnames_neurons)
my_levels <- c("n1", "n2", "n3", "n4", "n5", "n6", "n7", "n8", "n9", "n10", "n11",
    "n12", "n13")
vsns@active.ident <- factor(x = vsns@active.ident, levels = my_levels)
vsns[["orig.ident"]] = vsns@active.ident
temp = vsns
my_levels <- c("n5", "n6", "n4", "n3", "n1", "n2", "n7", "n8", "n9", "n10", "n11",
    "n12", "n13")
# Relevel object@ident
Idents(temp) <- factor(Idents(temp), levels = my_levels)
p = DimPlot(vsns, label = T)

cairo_ps("Figure 9B.eps", height = 4, width = 9, pointsize = 20)
cells_toplot = WhichCells(vsns, idents = c("n5", "n6", "n4", "n3", "n2", "n1"))
DoHeatmap(temp, draw.lines = FALSE, cells = cells_toplot, features = c("Neurod1",
    "Neurog1", "Gap43", "Gnao1", "Vmn2r1", "Vmn2r2", "Sdf2l1", "Creld2", "Dnajc3",
    "Hspa5", "Pdia6", "Manf")) + scale_fill_gradientn(colours = rev(brewer.pal(n = 9,
    name = "RdBu")))
dev.off()

```

### Creating Shiny app
```r
library(Seurat)
library(ShinyCell)
vno_scrna = readRDS("vno_scrna.rds")

new.cluster.IDs <- c("Solitary chemosensory cells", "Endothelial cells-1", "Gnai2 neurons",
    "Sustentacular cells", "Gnao1 neurons", "Cells from glandular tissue", "Immature neurons", "Non-sensory epithelium",
    "Macrophage (Mac3)", "T-cells", "Macrophage (Mac2)", "Tmbim+ Gnai2 neurons",
    "Macrophage (Mac1)", "Progenitor cells", "Fibroblasts", "Olfactory ensheathing cells",
    "Neutrophils", "Endothelial cells-2")

names(new.cluster.IDs) <- levels(vno_scrna)
vno_scrna = RenameIdents(vno_scrna, new.cluster.IDs)
vno_scrna[["orig.ident"]] = vno_scrna@active.ident

vno_neurons =readRDS("vno_neurons.rds")
vsns = vno_neurons
# head(Idents(vsns))

newnames_neurons <- c("n11","n13","n8","n1","n12","n3","n2","n7","n6","n5","n10","n9","n4")
names(newnames_neurons)<- levels(vsns)
vsns <- RenameIdents(vsns, newnames_neurons)
my_levels <- c("n1", "n2", "n3", "n4", "n5", "n6", "n7", "n8", "n9", "n10", "n11",
    "n12", "n13")
vsns@active.ident <- factor(x = vsns@active.ident, levels = my_levels)
vsns[["orig.ident"]] = vsns@active.ident


vno_conf = createConfig(vno_scrna)
vno_conf = delMeta(vno_conf, c("seurat_clusters", "RNA_snn_res.0.3", "RNA_snn_res.0.2",
    "nFeature_RNA", "percent.mito", "nCount_RNA"))
vno_conf = modMetaName(vno_conf, meta.to.mod = c("orig.ident"), new.name = c("Cluster"))

makeShinyFiles(vno_scrna, vno_conf, gex.assay = "RNA", gex.slot = "data", gene.mapping = TRUE,
    shiny.prefix = "vno", shiny.dir = "~/vno/scrna/scvno_explorer", default.gene1 = "Vmn2r1",
    default.gene2 = "Vmn2r2", default.multigene = c("Gnao1", "Gnai2", "Trpc2", "Omp",
        "Gap43", "Neurod1", "Ascl1", "Vmn2r1", "Vmn2r2"), default.dimred = c("UMAP_1",
        "UMAP_2"))

vsn_conf = createConfig(vsns)
vsn_conf = delMeta(vsn_conf, c("seurat_clusters", "RNA_snn_res.0.3", "RNA_snn_res.0.2",
    "RNA_snn_res.0.4", "nFeature_RNA", "percent.mito", "nCount_RNA"))
vsn_conf = modMetaName(vsn_conf, meta.to.mod = c("orig.ident"), new.name = c("Cluster"))
makeShinyFiles(vsns, vsn_conf, gex.assay = "RNA", gex.slot = "data", gene.mapping = TRUE,
    shiny.prefix = "vsns", shiny.dir = "~/vno/scrna/scvno_explorer", default.gene1 = "Gnao1",
    default.gene2 = "Vmn2r1", default.multigene = c("Gnao1", "Gnai2", "Trpc2", "Omp",
        "Gap43", "Neurod1", "Ascl1", "Vmn2r1", "Vmn2r2"), default.dimred = c("UMAP_1",
        "UMAP_2"))
citation = list(title="Single cell transcriptomics of vomeronasal neuroepithelium reveals a differential endoplasmic reticulum environment amongst neuronal subtypes", link="http://dx.doi.org/10.7554/eLife.98250.1", DOI="10.7554/elife.98250.1", publisher="eLife Sciences Publications, Ltd", author="Devakinandan, GVS and Terasaki, Mark and Dani, Adish", year="2024")



makeShinyCodesMulti(shiny.title = "scVNO explorer", shiny.footnotes = citation, shiny.prefix = c("vno",
    "vsns"), shiny.headers = c("Vomeronasal Neuroepithelium", "Vomeronasal Sensory Neurons"),
    shiny.dir = "~/vno/scrna/scvno_explorer")
```    
