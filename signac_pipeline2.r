## https://github.com/timoast/signac/issues/169
## https://github.com/timoast/signac/issues/173
## https://github.com/timoast/signac/issues/872
## https://github.com/timoast/signac/issues/872

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

tsne = read.csv('MSK93202_hg38/outs/analysis/tsne/2_components/projection.csv', row.names=1)

metadata <- read.csv(
  file = "MSK93202/outs/singlecell.csv",
  header = TRUE,
  row.names = 1)

metadata.hg38 <- read.csv(
  file = "MSK93202_hg38/outs/singlecell.csv",
  header = TRUE,
  row.names = 1)

metadata = subset(metadata, is_GRCh38_cell_barcode==1)

hg38.annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(hg38.annotations) <- 'UCSC'
genome(hg38.annotations) <- "hg38"

counts <- Read10X_h5(filename = "MSK93202_hg38/outs/filtered_peak_bc_matrix.h5")

chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  fragments = 'MSK93202_hg38/outs/fragments.tsv.gz',
  annotation = hg38.annotations
)

## chrom_assay@seqinfo = chrom_assay@ranges@seqinfo
## counts <- counts[grepl("GRCh38", rownames(counts)), ]
## rownames(counts) <- gsub("GRCh38_", "", rownames(counts))

obj <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "ATAC",
  meta.data = metadata.hg38
)

obj[['ATAC']]
granges(obj)

Annotation(obj) <- hg38.annotations

DefaultAssay(obj) <- "ATAC"

obj <- obj[, rownames(metadata)] # filter mouse/doublet cells

# compute nucleosome signal score per cell
obj <- NucleosomeSignal(object = obj)
# compute TSS enrichment score per cell
obj <- TSSEnrichment(object = obj, fast=F)

# add blacklist ratio and fraction of reads in peaks
obj$pct_reads_in_peaks <- obj$peak_region_fragments / obj$passed_filters * 100
obj$blacklist_ratio <- obj$blacklist_region_fragments / obj$peak_region_fragments

obj$high.tss <- ifelse(obj$TSS.enrichment > 2, 'High', 'Low')

pdf("MSK93202_QC.pdf", width=18)
VlnPlot(
  object = obj,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  ncol = 5,
  pt.size = 0
)
TSSPlot(obj, assay="ATAC", group.by = 'high.tss') + NoLegend()
obj$nucleosome_group <- ifelse(obj$nucleosome_signal > 2, 'NS > 0.75', 'NS < 0.75')
FragmentHistogram(object = obj, group.by = 'nucleosome_group')
dev.off()

obj <- subset(
  x = obj,
  subset = peak_region_fragments > 2000 &
    peak_region_fragments < 10000 &
    pct_reads_in_peaks > 40 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 2 &
    TSS.enrichment > 2
)
obj

obj <- RunTFIDF(obj)
obj <- FindTopFeatures(obj, min.cutoff = 'q0')
obj <- RunSVD(obj)

DepthCor(obj)
dev.off()

obj <- RunUMAP(object = obj, reduction = 'lsi', dims = 2:30)
obj <- FindNeighbors(object = obj, reduction = 'lsi', dims = 2:30)
obj <- FindClusters(object = obj, verbose = FALSE, algorithm = 3)

pdf("MSK93202_UMAP.pdf")
p1 = DimPlot(object = obj, label = TRUE) + NoLegend()
print(p1)
dev.off()

gene.activities <- GeneActivity(obj)

obj[['RNA']] <- CreateAssayObject(counts = gene.activities)
obj <- NormalizeData(
  object = obj,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(obj$nCount_RNA)
)

DefaultAssay(obj) <- 'RNA'

pdf("MSK93202_gene_activity.pdf", width=28, height=32)
FeaturePlot(
  object = obj,
  features = c("CD44", "EGFR", "MEOX2", "OGN", "POSTN", # EMT
               "MYOG", "MYOD1", "MYC", "MYCN", "MAX", "MYL4", "TNNT3", "TNNT2", "MEF2C", "DES",
               "BUB3", "BUB1", "E2F1", "CDK1", "MKI67", "TOP2A",
               "NOTCH1", "NOTCH3"),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 5
)
dev.off()


# Load the pre-processed scRNA-seq data for OBJs
obj_rna = readRDS('scRNA_seq/results/seurat_sara/MSK93202_hg19_seurat-object.rds')

meta = read.delim('scRNA_seq/clusters.txt',
                  row.names=1, header=T,
                  check.names=F, stringsAsFactors=F)

state = unlist(meta['MSK93202', , drop=T])
state = state[!((state=='') | (is.na(state)))]

## print(state)
levels(obj_rna$RNA_snn_res.0.8) = state

transfer.anchors <- FindTransferAnchors(
  reference = obj_rna,
  query = obj,
  reduction = 'cca'
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = obj_rna$RNA_snn_res.0.8,
  weight.reduction = obj[['lsi']],
  dims = 2:30
)

obj <- AddMetaData(object = obj, metadata = predicted.labels)

obj <- obj[, !(predicted.labels$predicted.id%in%c('UNIQUE', 'Interferon'))]

pdf("MSK93202_combined_umap.pdf", width=10, height=4.5)
plot1 <- DimPlot(
  object = obj_rna,
  group.by = 'RNA_snn_res.0.8',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')
plot2 <- DimPlot(
  object = obj,
  group.by = 'predicted.id',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')
print(plot1 + plot2)
dev.off()

library(JASPAR2020)
library(TFBSTools)

# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = 9606, all_versions = FALSE)
)

DefaultAssay(obj) <- "ATAC"

## https://github.com/timoast/signac/issues/463
gr <- granges(obj)
seq_keep <- seqnames(gr) %in% seqnames(BSgenome.Hsapiens.UCSC.hg38) 
seq_keep <- as.vector(seq_keep)
feat.keep <- GRangesToString(grange = gr[seq_keep])
obj[['ATAC']] <- subset(obj[["ATAC"]], features = feat.keep)

obj <- AddMotifs(
  object = obj,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)

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

obj = RegionStats(obj, genome = BSgenome.Hsapiens.UCSC.hg38)

# test enrichment
enriched.motifs <- FindMotifs(
  object = obj,
  features = top.da.peak
)

pdf("MSK93202_motif.pdf", width=28, height=18)
MotifPlot(
  object = obj,
  motifs = head(rownames(enriched.motifs), 50)
)
dev.off()

obj <- RunChromVAR(
  object = obj,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

DefaultAssay(obj) <- 'chromvar'

pdf("MSK93202_RBPJ.pdf", width=12, height=5)
p1 = DimPlot(object = obj, group.by='predicted.id', label = TRUE) + NoLegend()
p2 <- FeaturePlot(
  object = obj,
  features = rownames(enriched.motifs)[2],
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1
)
print(p1+p2)
dev.off()

pdf("MSK93202_MYOG.pdf", width=12, height=5)
p1 = DimPlot(object = obj, group.by='predicted.id', label = TRUE) + NoLegend()
p2 <- FeaturePlot(
  object = obj,
  features = rownames(enriched.motifs)[enriched.motifs$motif.name=='MYOG'],
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1
)
p1+p2
dev.off()

pdf("MSK93202_MYOD1.pdf", width=12, height=5)
p2 <- FeaturePlot(
  object = obj,
  features = rownames(enriched.motifs)[enriched.motifs$motif.name=='MYOD1'],
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1
)
p1+p2
dev.off()

pdf("MSK93202_MYC.pdf", width=12, height=5)
p2 <- FeaturePlot(
  object = obj,
  features = rownames(enriched.motifs)[which(enriched.motifs$motif.name == 'MYC')],
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1
)
p1+p2
dev.off()

pdf("MSK93202_CEBPB.pdf", width=12, height=5)
p2 <- FeaturePlot(
  object = obj,
  features = rownames(enriched.motifs)[which(enriched.motifs$motif.name == 'CEBPB')],
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1
)
p1+p2
dev.off()

DefaultAssay(obj) <- 'ATAC'
# gather the footprinting information for sets of motifs
obj <- Footprint(
  object = obj,
  motif.name = c("MYC", "MYOD1", "MAX", "RBPJ", "MYCN", "MEF2C", "CEBPB", "MYOG"),
  genome = BSgenome.Hsapiens.UCSC.hg38,
  group.by = 'predicted.id',
)

pdf("MSK93202_footprint.pdf", width=5, height=18)
p2 <- PlotFootprint(obj, features = c("MYC", "MYOD1", "MAX", "RBPJ", "MYCN", "MEF2C", "CEBPB"), group.by='predicted.id')
p2 + patchwork::plot_layout(ncol = 1)
dev.off()

saveRDS(obj, "MSK93202_signac.rds")
## Inference of gene expression
