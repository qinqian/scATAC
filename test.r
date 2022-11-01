library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
set.seed(1234)

# load the RNA and ATAC data

counts <- Read10X_h5("AMG510_libraries/outs/filtered_feature_bc_matrix.h5")
fragpath <- "AMG510_libraries/outs/atac_fragments.tsv.gz"

# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"
genome(annotation) <- "hg38"

# create a Seurat object containing the RNA adata
obj <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)

# create ATAC assay and add it to the object
obj[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotation
)


DefaultAssay(obj) <- "ATAC"

obj <- NucleosomeSignal(obj)
obj <- TSSEnrichment(obj)

pdf("QC.pdf", width=18)
VlnPlot(
  object = obj,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0
)
dev.off()


obj <- subset(
  x = obj,
  subset = nCount_ATAC < 100000 &
    nCount_RNA < 25000 &
    nCount_ATAC > 1000 &
    nCount_RNA > 1000 &
    nucleosome_signal < 2 &
    TSS.enrichment > 1
)
obj

peaks <- CallPeaks(obj, macs2.path = "~/miniconda3/envs/archr/bin/macs2")

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(obj),
  features = peaks,
  cells = colnames(obj)
)

# create a new assay using the MACS2 peak set and add it to the Seurat object
obj[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = fragpath,
  annotation = annotation
)

DefaultAssay(obj) <- "RNA"
obj <- SCTransform(obj)
obj <- RunPCA(obj)
obj <- FindNeighbors(obj, dims = 1:30)
obj <- FindClusters(obj, resolution = 0.5, algorithm = 3)


DefaultAssay(obj) <- "peaks"
obj <- FindTopFeatures(obj, min.cutoff = 5)
obj <- RunTFIDF(obj)
obj <- RunSVD(obj)

obj <- FindMultiModalNeighbors(
  object = obj,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:50, 2:40),
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)

# build a joint UMAP visualization
obj <- RunUMAP(
  object = obj,
  nn.name = "weighted.nn",
  assay = "RNA",
  verbose = TRUE
)

DimPlot(obj, label = TRUE, repel = TRUE, reduction = "umap") + NoLegend()
dev.off()

DefaultAssay(obj) <- "peaks"

# first compute the GC content for each peak
obj <- RegionStats(obj, genome = BSgenome.Hsapiens.UCSC.hg38)

# link peaks to genes
obj <- LinkPeaks(
  object = obj,
  peak.assay = "peaks",
  expression.assay = "SCT",
  genes.use = c("LYZ", "MS4A1")
)
