metacolors <- c(rgb(119, 62, 20, maxColorValue = 255),
                rgb(236, 133, 40, maxColorValue = 255),
                rgb(59, 22, 115, maxColorValue = 255),
                rgb(52, 101, 252, maxColorValue = 255),
                rgb(242, 242, 242, maxColorValue = 255),
                rgb(52, 101, 252, maxColorValue = 255),
                rgb(225, 39, 39, maxColorValue = 255),
                rgb(72, 159,  75, maxColorValue = 255),
                rgb(20, 64, 216, maxColorValue = 255),
                rgb(226, 75, 143, maxColorValue = 255),
                rgb(158, 60, 200, maxColorValue = 255),
                rgb(241, 250, 100, maxColorValue = 255))
metalabels <- c("Ground", "Hypoxia", "EMT",
                "G1S", "UNASSIGNED",
                "G2M",  "Muscle", "Interferon", "Prolif",
                "Histone", "Apoptosis", 'UPR')
names(metacolors) <- metalabels

library(GenomicRanges)
library(JASPAR2020)
library(TFBSTools)
library(ggplot2)
library(patchwork)
set.seed(1234)
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v75)
library(EnsDb.Hsapiens.v86)
## library(EnsDb.Mmusculus.v79)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)

## obj = readRDS("MSK93202_signac.rds")
obj = readRDS("MSK93202_signac_motif_cores_vertebrates.rds")
obj$predicted.id[obj$predicted.id=='Hypoxia'] = 'Ground'
obj$predicted.id[obj$predicted.id=='Interferon'] = 'Ground'

obj_rna = readRDS('scRNA_seq/results/seurat_sara/MSK93202_hg19_seurat-object.rds')

meta = read.delim('scRNA_seq/clusters.txt',
                  row.names=1, header=T,
                  check.names=F, stringsAsFactors=F)

state = unlist(meta['MSK93202', , drop=T])
state = state[!((state=='') | (is.na(state)))]

state[state=='Hypoxia'] = 'Ground'
state[state=='Interferon'] = 'Ground'
state[state=='Unique #1'] = 'Unique'

## print(state)
levels(obj_rna$RNA_snn_res.0.8) = state
DefaultAssay(obj) <- "RNA"
transfer.anchors <- FindTransferAnchors(
  reference = obj_rna,
  query = obj,
  reduction = 'cca'
)
rna <- TransferData(
  anchorset = transfer.anchors,
  refdata = GetAssayData(obj_rna, assay = "RNA", slot = "data"),
  weight.reduction = obj[["lsi"]],
  dims = 2:30
)
# add predicted values as a new assay
obj[["predicted"]] <- rna

DefaultAssay(obj) <- "ATAC" 

## https://github.com/timoast/signac/issues/463
gr <- granges(obj)
seq_keep <- seqnames(gr) %in% seqnames(BSgenome.Hsapiens.UCSC.hg38) 
seq_keep <- as.vector(seq_keep)
feat.keep <- GRangesToString(grange = gr[seq_keep])
obj[['ATAC']] <- subset(obj[["ATAC"]], features = feat.keep)


# Get a list of motif position frequency matrices from the JASPAR database
## pfm <- getMatrixSet(
##   x = JASPAR2020,
##   opts = list(species = 9606, all_versions = FALSE)
## )

pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

## obj <- AddMotifs(
##   object = obj,
##   genome = BSgenome.Hsapiens.UCSC.hg38,
##   pfm = pfm
## )

## obj <- RunChromVAR(
##   object = obj,
##   genome = BSgenome.Hsapiens.UCSC.hg38,
## )

## saveRDS(obj, "MSK93202_signac_motif_cores_vertebrates.rds")

## obj = readRDS("MSK93202_signac_motif_cores_vertebrates.rds")

DefaultAssay(obj) <- "ATAC"
da_peaks <- FindMarkers(
    object = obj,
    group.by = 'predicted.id',
    ident.1 = 'EMT',
    ident.2 = 'Ground',
    only.pos = TRUE,
    test.use = 'LR',
    min.pct = 0.05,
    latent.vars = 'nCount_ATAC'
)
top.da.peak <- rownames(da_peaks[da_peaks$p_val < 0.005, ])

enriched.motifs <- FindMotifs(
  object = obj,
  features = top.da.peak
)

da_peaks.muscle <- FindMarkers(
    object = obj,
    group.by = 'predicted.id',
    ident.1 = 'Muscle',
    ident.2 = 'Ground',
    only.pos = TRUE,
    test.use = 'LR',
    min.pct = 0.05,
    latent.vars = 'nCount_ATAC'
)
top.da.peak.muscle <- rownames(da_peaks.muscle[da_peaks.muscle$p_val < 0.005, ])
enriched.motifs.muscle <- FindMotifs(
  object = obj,
  features = top.da.peak.muscle
)

pdf('r01_motif_analysis_singlecell_msk93202.pdf')
DefaultAssay(obj) <- "predicted"
DimPlot(
  object = obj,
  group.by = 'predicted.id',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')
DefaultAssay(obj) <- "chromvar"
FeaturePlot(
  object = obj,
  features = rownames(enriched.motifs)[which(enriched.motifs$motif.name %in% c('CEBPB', "RBPJ", "MEF2C", "LHX2"))],
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1
)
DefaultAssay(obj) <- "predicted"
FeaturePlot(
  object = obj,
  features = c("NOTCH1", "NOTCH3", "CEBPB", "MEF2C", "LHX2"),
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1
)
dev.off()


pdf('r01_motif_analysis_singlecell_msk93202_cleanup.pdf', width=18, height=12)
p01=DimPlot(
  object = obj_rna,
  group.by = 'RNA_snn_res.0.8', cols=metacolors,
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')
p02=FeaturePlot(
  object = obj_rna,
  features = c("NOTCH3"),
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1
)
p03=FeaturePlot(
  object = obj_rna,
  features = c("CEBPB"), 
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1
)
p04=FeaturePlot(
  object = obj_rna,
  features = c("MEF2C"), 
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1
)
p05=FeaturePlot(
  object = obj_rna,
  features = c("LHX2"), 
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1
)
DefaultAssay(obj) <- "predicted"
p1=DimPlot(
  object = obj,
  group.by = 'predicted.id',  cols=metacolors,
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('imputed RNA')
p2=FeaturePlot(
  object = obj,
  features = c("NOTCH3"),
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1
)
p3=FeaturePlot(
  object = obj,
  features = c("CEBPB"), 
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1
)
p4=FeaturePlot(
  object = obj,
  features = c("MEF2C"), 
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1
)
p5=FeaturePlot(
  object = obj,
  features = c("LHX2"), 
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1
)
DefaultAssay(obj) <- "chromvar"
p6=FeaturePlot(
  object = obj,
  features = rownames(enriched.motifs)[which(enriched.motifs$motif.name %in% c("RBPJ"))],
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1
)
p7=FeaturePlot(
  object = obj,
  features = rownames(enriched.motifs)[which(enriched.motifs$motif.name %in% c("CEBPB"))],
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1
)
p8=FeaturePlot(
  object = obj,
  features = rownames(enriched.motifs)[which(enriched.motifs$motif.name %in% c("MEF2C"))],
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1
)
p9=FeaturePlot(
  object = obj,
  features = rownames(enriched.motifs)[which(enriched.motifs$motif.name %in% c("LHX2"))],
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1
)
p11=DimPlot(
  object = obj,
  group.by = 'predicted.id',  cols=metacolors,
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')
print((p01+p02+p03+p04+p05+plot_layout(ncol = 5))/(p1+p2+p3+p4+p5+plot_layout(ncol = 5))/(p11+p6+p7+p8+p9+plot_layout(ncol = 5)))
dev.off()

obj2 = readRDS("MSK74711_signac.rds")

obj2$predicted.id[obj2$predicted.id=='Hypoxia'] = 'Ground'
obj2$predicted.id[obj2$predicted.id=='Interferon'] = 'Ground'

gr <- granges(obj2)
seq_keep <- seqnames(gr) %in% seqnames(BSgenome.Hsapiens.UCSC.hg38) 
seq_keep <- as.vector(seq_keep)
feat.keep <- GRangesToString(grange = gr[seq_keep])
obj2[['ATAC']] <- subset(obj2[["ATAC"]], features = feat.keep)

obj_rna2 = readRDS('/PHShome/qq06/langenau/geo_upload/geo_submission_scrms_june21/processed_data_files/MSK74711_hg19_mm10/20191031_MSK74711_seurat-object.rds')

meta = read.delim('/PHShome/qq06/langenau/projects/01_sc_rms/phaseD_comparestjde/../final_annotations/Final_clusters.txt', sep='\t', row.names=1, header=T,
                  check.names=F, stringsAsFactors=F)

state = unlist(meta['MSK74711', , drop=T])
state = state[!((state=='') | (is.na(state)))]
state[state=='Hypoxia'] = 'Ground'
state[state=='Interferon'] = 'Ground'
state[state=='Unique #1'] = 'Unique'


levels(obj_rna2$RNA_snn_res.0.8) = state

DefaultAssay(obj2) <- "RNA"
transfer.anchors <- FindTransferAnchors(
  reference = obj_rna2,
  query = obj2,
  reduction = 'cca'
)

rna2 <- TransferData(
  anchorset = transfer.anchors,
  refdata = GetAssayData(obj_rna2, assay = "RNA", slot = "data"),
  weight.reduction = obj2[["lsi"]],
  dims = 2:30
)

# add predicted values as a new assay
obj2[["predicted"]] <- rna2

## DefaultAssay(obj2) <- "ATAC"
## obj2 <- AddMotifs(
##   object = obj2,
##   genome = BSgenome.Hsapiens.UCSC.hg38,
##   pfm = pfm
## )

## obj2 <- RunChromVAR(
##   object = obj2,
##   genome = BSgenome.Hsapiens.UCSC.hg38,
## )

## saveRDS(obj2, "MSK74711_signac_motif_cores_vertebrates.rds")
## obj2 = readRDS("MSK74711_signac_motif_cores_vertebrates.rds")

DefaultAssay(obj2) <- "ATAC"
da_peaks2 <- FindMarkers(
    object = obj2,
    group.by = 'predicted.id',
    ident.1 = 'EMT',
    ident.2 = 'Ground',
    only.pos = TRUE,
    test.use = 'LR',
    min.pct = 0.05,
    latent.vars = 'nCount_ATAC'
)

top.da.peak2 <- rownames(da_peaks2[da_peaks2$p_val < 0.005, ])

enriched.motifs2 <- FindMotifs(
  object = obj2,
  features = top.da.peak2
)

pdf('r01_motif_analysis_singlecell_msk74711.pdf')
DefaultAssay(obj2) <- "predicted"
DimPlot(
  object = obj2,
  group.by = 'predicted.id',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')
DefaultAssay(obj2) <- "chromvar"
FeaturePlot(
  object = obj2,
  features = rownames(enriched.motifs2)[which(enriched.motifs2$motif.name %in% c('CEBPB', "RBPJ", "MEF2C", "LHX2"))],
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1
)
DefaultAssay(obj2) <- "predicted"
FeaturePlot(
  object = obj2,
  features = c("NOTCH1", "NOTCH3", "CEBPB", "MEF2C", "LHX2"),
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1
)
dev.off()

pdf('r01_motif_analysis_singlecell_msk74711_cleanup.pdf', width=18, height=12)
p01=DimPlot(
  object = obj_rna2,
  group.by = 'RNA_snn_res.0.8', cols=metacolors,
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')
p02=FeaturePlot(
  object = obj_rna2,
  features = c("NOTCH3"),
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1
)
p03=FeaturePlot(
  object = obj_rna2,
  features = c("CEBPB"), 
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1
)
p04=FeaturePlot(
  object = obj_rna2,
  features = c("MEF2C"), 
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1
)
p05=FeaturePlot(
  object = obj_rna2,
  features = c("LHX2"), 
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1
)
DefaultAssay(obj2) <- "predicted"
p1=DimPlot(
  object = obj2,
  group.by = 'predicted.id',  cols=metacolors,
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('imputed RNA')
p2=FeaturePlot(
  object = obj2,
  features = c("NOTCH3"),
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1
)
p3=FeaturePlot(
  object = obj2,
  features = c("CEBPB"), 
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1
)
p4=FeaturePlot(
  object = obj2,
  features = c("MEF2C"), 
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1
)
p5=FeaturePlot(
  object = obj2,
  features = c("LHX2"), 
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1
)
DefaultAssay(obj2) <- "chromvar"
p6=FeaturePlot(
  object = obj2,
  features = rownames(enriched.motifs2)[which(enriched.motifs2$motif.name %in% c("RBPJ"))],
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1
)
p7=FeaturePlot(
  object = obj2,
  features = rownames(enriched.motifs2)[which(enriched.motifs2$motif.name %in% c("CEBPB"))],
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1
)
p8=FeaturePlot(
  object = obj2,
  features = rownames(enriched.motifs2)[which(enriched.motifs2$motif.name %in% c("MEF2C"))],
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1
)
p9=FeaturePlot(
  object = obj2,
  features = rownames(enriched.motifs2)[which(enriched.motifs2$motif.name %in% c("LHX2"))],
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1
)
p11=DimPlot(
  object = obj2,
  group.by = 'predicted.id',  cols=metacolors,
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')
print((p01+p02+p03+p04+p05+plot_layout(ncol = 5))/(p1+p2+p3+p4+p5+plot_layout(ncol = 5))/(p11+p6+p7+p8+p9+plot_layout(ncol = 5)))
dev.off()


pdf('r01_singlecell_dave_genes_3d.pdf', width=16, height=9.6)
genes=c("MEF2B", "MEF2C", "MYOG", "MYOD1", "MYF5", "NOTCH1", "NOTCH3", "RBPJL", "HES1", "FZD4", "SNAI1", "SNAI2", "MYC", "MYCN", "CEBPZ", "CEBPB", "ID1", "TNNT3", "TNNI2", "TNNT3", "CD44", "CD90", "LRRN2", "DES", "SPARC", "OGN", "MGP")
motifs=c("MEF2B", "MEF2C", "MYOG", "MYOD1", "MYF5", "RBPJ", "RBPJ", "Rbpjl", "HES1", "", "SNAI1", "SNAI2", "MYC", "MYCN", "", "CEBPB", "", "", "", "", "", "", "", "", "", "", "")
p1=DimPlot(
  object = obj,
  group.by = 'predicted.id',  cols=metacolors,
  label = F,
  repel = TRUE) + NoLegend() + ggtitle('imputed RNA')
p01=DimPlot(
    object = obj_rna,
    group.by = 'RNA_snn_res.0.8', cols=metacolors,
    label = F,
    repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')
p_empty = ggplot() + theme_void()
for (g in seq_along(genes)) {
    if (length(intersect(genes[g], rownames(obj_rna)))>0) {
        p02=FeaturePlot(
            object = obj_rna,
            features = genes[g],
            min.cutoff = 'q10',
            max.cutoff = 'q90',
            pt.size = 0.1)
    } else {
        p02 = ggplot() + theme_void()
    }
    DefaultAssay(obj) <- "predicted"
    if (length(intersect(genes[g], rownames(obj)))>0) {
        p2=FeaturePlot(
            object = obj,
            features = genes[g],
            min.cutoff = 'q10',
            max.cutoff = 'q90',
            pt.size = 0.1
        )
    } else {
        p2 = ggplot() + theme_void()
    }
    DefaultAssay(obj) <- "chromvar"
    if (length(which(enriched.motifs$motif.name %in% motifs[g]))>0) {
        p6=FeaturePlot(
            object = obj,
            features = rownames(enriched.motifs)[which(enriched.motifs$motif.name %in% motifs[g])],
            min.cutoff = 'q10',
            max.cutoff = 'q90',
            pt.size = 0.1
        )
    } else {
        p6 = ggplot() + theme_void()
    }
    print((p01+p02+p_empty)/(p1+p2+p6))
}
dev.off()

pdf('r01_singlecell_dave_genes_3d_msk74711.pdf', width=16, height=9.6)
genes=c("MEF2B", "MEF2C", "MYOG", "MYOD1", "MYF5", "NOTCH1", "NOTCH3", "RBPJL", "HES1", "FZD4", "SNAI1", "SNAI2", "MYC", "MYCN", "CEBPZ", "CEBPB", "ID1", "TNNT3", "TNNI2", "TNNT3", "CD44", "CD90", "LRRN2", "DES", "SPARC", "OGN", "MGP")
motifs=c("MEF2B", "MEF2C", "MYOG", "MYOD1", "MYF5", "RBPJ", "RBPJ", "Rbpjl", "HES1", "", "SNAI1", "SNAI2", "MYC", "MYCN", "", "CEBPB", "", "", "", "", "", "", "", "", "", "", "")
p1=DimPlot(
  object = obj2,
  group.by = 'predicted.id',  cols=metacolors,
  label = F,
  repel = TRUE) + ggtitle('imputed RNA')
p01=DimPlot(
    object = obj_rna2,
    group.by = 'RNA_snn_res.0.8', cols=metacolors,
    label = F,
    repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')
p_empty = ggplot() + theme_void()
for (g in seq_along(genes)) {
    if (length(intersect(genes[g], rownames(obj_rna2)))>0) {
        p02=FeaturePlot(
            object = obj_rna2,
            features = genes[g],
            min.cutoff = 'q10',
            max.cutoff = 'q90',
            pt.size = 0.1)
    } else {
        p02 = ggplot() + theme_void()
    }
    DefaultAssay(obj2) <- "predicted"
    if (length(intersect(genes[g], rownames(obj2)))>0) {
        p2=FeaturePlot(
            object = obj2,
            features = genes[g],
            min.cutoff = 'q10',
            max.cutoff = 'q90',
            pt.size = 0.1
        )
    } else {
        p2 = ggplot() + theme_void()
    }
    DefaultAssay(obj2) <- "chromvar"
    if (length(which(enriched.motifs2$motif.name %in% motifs[g]))>0) {
        p6=FeaturePlot(
            object = obj2,
            features = rownames(enriched.motifs2)[which(enriched.motifs2$motif.name %in% motifs[g])],
            min.cutoff = 'q10',
            max.cutoff = 'q90',
            pt.size = 0.1
        )
    } else {
        p6 = ggplot() + theme_void()
    }
    print((p01+p02+p_empty)/(p1+p2+p6))
}
dev.off()


pdf('r01_motif_analysis_singlecell_msk93202_cleanup_cellstate.pdf', width=5, height=5)
p1=DimPlot(
  object = obj,
  group.by = 'predicted.id',  cols=metacolors,
  label = F,
  repel = TRUE) + NoLegend() + ggtitle('imputed RNA')
print(p1)
dev.off()
