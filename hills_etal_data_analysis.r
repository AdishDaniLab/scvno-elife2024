library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(SoupX)
library(RColorBrewer)
library(forcats)
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00",
    "#CC79A7")

ronyu_adult_soupX.input=load10X('../ronyu_data/run/adult_vno_combined/outs/count/')
ronyu_adult_soupX.input=autoEstCont(ronyu_adult_soupX.input)
ronyu_adult_soupX.out=adjustCounts(ronyu_adult_soupX.input)
ronyu_adult_bc=CreateSeuratObject(ronyu_adult_soupX.out)


mito.genes <- grep(pattern = "^mt-", x = rownames(ronyu_adult_bc@assays[["RNA"]]), value = TRUE)
ronyu_adult_bc[["percent.mt"]] <- PercentageFeatureSet(ronyu_adult_bc, pattern = "^mt-")
p=VlnPlot(ronyu_adult_bc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(p, file="test.png")
ronyu_adult_bc_subset <- subset(x = ronyu_adult_bc, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 5) # 5% mt-dna
ronyu_adult_bc_subset <- NormalizeData(object = ronyu_adult_bc_subset, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(ronyu_adult_bc_subset)
ronyu_adult_bc_subset=ScaleData(object = ronyu_adult_bc_subset, features=all.genes, vars.to.regress = c("nCount_RNA", "percent.mt"))

ronyu_adult_bc_subset=FindVariableFeatures(object=ronyu_adult_bc_subset, mean.function=ExpMean, dispersion.function=LogVMR, x.low.cutoff =0.0125, x.high.cutoff=3, y.cutoff=0.5, nfeatures=2000)

ronyu_adult_bc_subset <- RunPCA(object = ronyu_adult_bc_subset,  npcs = 50, verbose = FALSE)
#ronyu_adult_bc_subset <- JackStraw(object = ronyu_adult_bc_subset, reduction = "pca", dims = 50, num.replicate = 100,  prop.freq = 0.1, verbose = FALSE)  
#ronyu_adult_bc_subset <- ScoreJackStraw(object = ronyu_adult_bc_subset, dims = 1:50, reduction = "pca")
#JackStrawPlot(object = ronyu_adult_bc_subset, dims = 1:50, reduction = "pca")
p=ElbowPlot(object = ronyu_adult_bc_subset,ndims=50) # used to determine the number of principal components
ggsave(p, file="test.png")
ronyu_adult_bc_subset <- FindNeighbors(ronyu_adult_bc_subset, reduction = "pca", dims = 1:36) #dims decided from ElbowPlot
ronyu_adult_bc_subset <- FindClusters(ronyu_adult_bc_subset, resolution = 0.7) #resolution decides number of clusters
ronyu_adult_bc_subset <- RunUMAP(object = ronyu_adult_bc_subset, dims = 1:30)
p=DimPlot(ronyu_adult_bc_subset, label=T)& NoLegend()
ggsave(p, file="umap.png")
ronyu_adult_bc_subset.markers <- FindAllMarkers(ronyu_adult_bc_subset)
top10 <- ronyu_adult_bc_subset.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC)
write.table(top10,file="top10markers_log2FC.tsv", sep="\t")
p=FeaturePlot(ronyu_adult_bc_subset, features=c("Gnao1", "Gnai2", "Omp", "Trpc2", "Gap43", "Neurod1", "Neurog1", "Ascl1","Obp2a","Tmbim1","Cbr2","Ppic")) & NoLegend()
ggsave(p, file="featureplot.png", height=14, width=14, units="in")
saveRDS(ronyu_adult_bc_subset,file="ronyu_prelimdata.rds")

ronyu_adult_bc_subset <- readRDS("ronyu_prelimdata.rds")


p=VlnPlot(ronyu_adult_bc_subset, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(p, file="test.png",height=5, width=17, units="in")

hbbbs = names(which(ronyu_adult_bc_subset[["RNA"]]$data["Hbb-bs", ] > 0))
ronyu_adult_bc_subset = subset(ronyu_adult_bc_subset, cells = hbbbs, invert = TRUE)


#Gnao1 analysis
gnao1_neurons=subset(ronyu_adult_bc_subset,idents=c(6,12,17))
p=FeaturePlot(gnao1_neurons, features=c("H2-M10.1","H2-M10.3","Obp2a")) & NoLegend()
ggsave(p, file="test.png")

p=DimPlot(gnao1_neurons) & NoLegend()
q= FeaturePlot(gnao1_neurons, features=c("Vmn2r1", "Vmn2r2", "H2-M10.3")) & NoLegend()
ggsave(p+q, file="gnao1_neurons.png",height=7, width=14, units="in")



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

vmn1r = grep(pattern = "^Vmn1r", x = rownames(ronyu_adult_bc_subset@assays[["RNA"]]),
    value = TRUE)  # getting gene names of al V1Rs
vmn2r = grep(pattern = "^Vmn2r", x = rownames(ronyu_adult_bc_subset@assays[["RNA"]]),
    value = TRUE)  # getting gene names of all V2Rs
famc_vmn2r = c("Vmn2r1", "Vmn2r2", "Vmn2r3", "Vmn2r4", "Vmn2r5", "Vmn2r6", "Vmn2r7")
h2mv_genes = c("H2-M10.1", "H2-M10.3", "H2-M10.4", "H2-M11", "H2-M9","H2-M1","H2-M10.5", "H2-M10.6", "H2-M5", "H2-M3", "H2-M2")
abd_vmn2r = vmn2r[!vmn2r %in% famc_vmn2r]
abd_h2mv = c(abd_vmn2r, h2mv_genes)
vmn2r_h2mv=c(vmn2r,h2mv_genes)

#making data tables
famc_v2r_exp_gnao1_clus = data.frame(gnao1_neurons[["RNA"]]$data[famc_vmn2r, ])
abd_v2r_exp_gnao1_clus = data.frame(gnao1_neurons[["RNA"]]$data[abd_vmn2r, ])
abd_v2r_h2mv_exp_gnao1_clus = data.frame(gnao1_neurons[["RNA"]]$data[abd_h2mv, ])
v2r_h2mv_exp_gnao1_clus = data.frame(gnao1_neurons[["RNA"]]$data[vmn2r_h2mv, ])
h2mv_exp_gnao1_clus = data.frame(gnao1_neurons[["RNA"]]$data[h2mv_genes, ])
v2r_exp_gnao1_clus = data.frame(gnao1_neurons[["RNA"]]$data[vmn2r, ])

test = as.vector(unlist(v2r_exp_gnao1_clus))
nonzero_exp_values = data.frame(exp = test[which(test != 0)])
cairo_ps("v2r_expression_threshold.eps", height = 7, width = 7)
p=ggplot(nonzero_exp_values, aes(x = exp)) + geom_density(lwd = 1.5) + xlim(0, 10) +
    geom_vline(xintercept = 1.25, colour = "#BB0000", linetype = "dashed") + theme_clean() +
    theme(axis.text = element_text(size = 32))
ggsave(p, file="v2r_expression distribution.png")

abd_table = get_coexpression_combis(abd_v2r_exp_gnao1_clus, 1.25, "ABD")
famc_table = get_coexpression_combis(famc_v2r_exp_gnao1_clus, 1.25, "famc")
h2mv_table = get_coexpression_combis(h2mv_exp_gnao1_clus, 1.25, "h2mv")
all_v2r_table = get_coexpression_combis(v2r_exp_gnao1_clus, 1.25, "all_v2r")

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

cairo_ps("h2mv_coexpression_frequency.eps", height = 5, width = 3.75)
ggplot(df <- filter(total_gnao1_table, h2mv_confidence == "undefined" | h2mv_confidence ==
    "High confidence"), aes(x = factor(total_h2mv_genes), fill = factor(famc_type),
    group = factor(famc_type))) + geom_bar(width = 0.75) + scale_fill_manual(values = cbPalette) +
    # facet_wrap(~h2mv_confidence,nrow = 1)+
geom_text(aes(label = after_stat(count)), stat = "count", position = position_stack(0.5)) +
    theme_classic() + theme(text = element_text(size = 12), legend.position = "none",
    axis.title = element_text(face = "bold", size = 12)) + xlab("Number of H2-Mv genes\n expressed in a cell") +
    ylab("Number of cells")
dev.off()