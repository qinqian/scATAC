library(ArchR)
library('Cairo')
library(GetoptLong)
library(circlize)
library(DOSE)

addArchRThreads(threads = 10)
addArchRGenome("hg19")

setwd("~/langenau/stjude/scATAC-seq/")

proj = loadArchRProject("Save-Proj/")

proj2 <- filterDoublets(proj)

## focus on one example
mast111 = proj2[proj2$Sample=='SJRHB013758_X1', ]

proj_tmp <- addIterativeLSI(
    ArchRProj = mast111,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI", 
    iterations = 4, 
    clusterParams = list( #See Seurat::FindClusters
        resolution = c(0.2), 
        sampleCells = 10000, 
        n.start = 10
    ), 
    varFeatures = 25000, 
    dimsToUse = 1:30,
    force=T,
    )
proj_tmp <- addUMAP(
    ArchRProj = proj_tmp,
    reducedDims = "IterativeLSI", 
    name = "UMAP", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine",
    force=T
)
proj_tmp <- addClusters(
    input = proj_tmp,
    reducedDims = "IterativeLSI",
    method = "Seurat",
    name = "ClustersHarmony_individualsample",
    resolution = 0.8,
    force=T
)


p4 <- plotEmbedding(ArchRProj = proj_tmp, colorBy = "cellColData", name = "ClustersHarmony_individualsample", embedding = "UMAP")
plotPDF(p4, name = "MAST111_RMS_ATAC_Plot-UMAP-Clusters.pdf", ArchRProj = proj_tmp, addDOC = FALSE, width = 5, height = 5)


mast111 = proj_tmp[proj_tmp$ClustersHarmony_individualsample!='C1', ]

proj_tmp <- addIterativeLSI(
    ArchRProj = mast111,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI", 
    iterations = 4, 
    clusterParams = list( #See Seurat::FindClusters
        resolution = c(0.2), 
        sampleCells = 10000, 
        n.start = 10
    ), 
    varFeatures = 25000, 
    dimsToUse = 1:30,
    force=T,
    )
proj_tmp <- addUMAP(
    ArchRProj = proj_tmp,
    reducedDims = "IterativeLSI", 
    name = "UMAP", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine",
    force=T
)
proj_tmp <- addClusters(
    input = proj_tmp,
    reducedDims = "IterativeLSI",
    method = "Seurat",
    name = "ClustersHarmony_individualsample",
    resolution = 0.8,
    force=T
)

p4 <- plotEmbedding(ArchRProj = proj_tmp, colorBy = "cellColData", name = "ClustersHarmony_individualsample", embedding = "UMAP")
plotPDF(p4, name = "MAST111_RMS_ATAC_Plot-UMAP-Clusters-RemoveC1.pdf", ArchRProj = proj_tmp, addDOC = FALSE, width = 5, height = 5)


markersGS <- getMarkerFeatures(
    ArchRProj = proj_tmp, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "ClustersHarmony_individualsample",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")

## refine annotation of cell state
library(clusterProfiler)
library(GetoptLong)
library(circlize)
library(DOSE)

allgenes = readRDS("allgenes_signatures.rds")

clusters = list()
for (clust in names(markerList)) {
    print(dim(markerList[[clust]]))
    genes = markerList[[clust]]$name
    y1 = enricher(genes, TERM2GENE=allgenes[,-2], maxGSSize=1000, minGSSize = 5)
    if (is.null(y1))
        clusters[[as.character(clust)]] <- y1
    else{
        y1 = y1@result
        rownames(y1) = paste0(clust, ".", rownames(y1))
        y1$clusters = clust
        clusters[[as.character(clust)]] <- y1
    }
}
clusters = do.call('rbind', clusters)
clusters = clusters %>% mutate(FoldEnrichment = parse_ratio(GeneRatio)/parse_ratio(BgRatio))
clusters = subset(clusters, p.adjust <= 0.05 & qvalue <= 0.05 & Count >= 3 & FoldEnrichment >= 1.1)
readr::write_csv(clusters, path='MAST111_RemoveC1_StJude_cluster_enrichr_filter_annotation.csv')

## leverage our MAST111 scRNA-seq
MAST111 = readRDS('/data/langenau/human_rms_pdxs/seurat_objects/20190624_seurat-object_MAST111.rds')

meta = read.delim('/PHShome/qq06/langenau/projects/01_sc_rms/phaseD_comparestjde/../final_annotations/Final_clusters.txt', sep='\t', row.names=1, header=T,
                  check.names=F, stringsAsFactors=F)

state = unlist(meta['MAST111', , drop=T])
state = state[!((state=='') | (is.na(state)))]

## print(state)
levels(MAST111$RNA_snn_res.0.8) = state

proj_tmp <- addGeneIntegrationMatrix(
    ArchRProj = proj_tmp,
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = MAST111,
    addToArrow = FALSE,
    groupRNA = "RNA_snn_res.0.8",
    nameCell = "predictedCell_Un",
    nameGroup = "predictedGroup_Un",
    nameScore = "predictedScore_Un"
)

p4 <- plotEmbedding(ArchRProj = proj_tmp, colorBy = "cellColData", name = "predictedGroup_Un", embedding = "UMAP")
plotPDF(p4, name = "MAST111_RMS_ATAC_Plot-UMAP-Clusters-RemoveC1-annotation-by-rna.pdf", ArchRProj = proj_tmp, addDOC = FALSE, width = 5, height = 5)

cM <- as.matrix(confusionMatrix(proj_tmp$ClustersHarmony_individualsample, proj_tmp$predictedGroup_Un))

preClust <- colnames(cM)[apply(cM, 1, which.max)]

preClust = cbind(preClust, rownames(cM))

clustMuscleEMTProlif <- preClust[grep("Muscle|EMT|Prolif", preClust[ ,1]), 2]

groupList <- SimpleList(
    MEP = SimpleList(
        ATAC = proj_tmp$cellNames[proj_tmp$ClustersHarmony_individualsample %in% clustMuscleEMTProlif],
        RNA = colnames(MAST111)[grep("Muscle|EMT|Prolif", MAST111$RNA_snn_res.0.8)]
    ),
    NonMEP = SimpleList(
        ATAC = proj_tmp$cellNames[!(proj_tmp$ClustersHarmony_individualsample %in% clustMuscleEMTProlif)],
        RNA = colnames(MAST111)[-grep("Muscle|EMT|Prolif", MAST111$RNA_snn_res.0.8)]
    )    
)

proj_tmp <- addGeneIntegrationMatrix(
    ArchRProj = proj_tmp,
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = MAST111,
    addToArrow = FALSE, 
    groupList = groupList,
    groupRNA = "RNA_snn_res.0.8",
    nameCell = "predictedCell_Co",
    nameGroup = "predictedGroup_Co",
    nameScore = "predictedScore_Co"
)

p4 <- plotEmbedding(ArchRProj = proj_tmp, colorBy = "cellColData", name = "predictedGroup_Co", embedding = "UMAP")
plotPDF(p4, name = "MAST111_RMS_ATAC_Plot-UMAP-CoClusters-RemoveC1-annotation-by-rna.pdf", ArchRProj = proj_tmp, addDOC = FALSE, width = 5, height = 5)


p <- plotEmbedding(
    ArchRProj = proj_tmp, 
    colorBy = "GeneScoreMatrix", 
    name = c("CD44", "EGFR", "MEOX2", "OGN", "POSTN", # EMT
             "MYOG", "MYL4", "TNNT3", "TNNT2", "MEF2C", "DES", # Muscle
             "BUB3", "BUB1", "E2F1", "CDK1", "MKI67", "TOP2A",
             "NOTCH1", "NOTCH3"),
    embedding = "UMAP",
    quantCut = c(0.01, 0.99),
    imputeWeights = NULL
)

plotPDF(plotList = p, 
        name = "MAST111-Plot-UMAP-Marker-Genes-WO-Imputation.pdf",
        ArchRProj = proj_tmp,
        addDOC = FALSE, width = 5, height = 5)

p <- plotBrowserTrack(
    ArchRProj = proj_tmp, 
    groupBy = "predictedGroup_Co",
    geneSymbol = c("CD44", "EGFR", "MEOX2", "OGN", "POSTN", # EMT
                   "MYOG", "MYL4", "TNNT3", "TNNT2", "MEF2C", "DES", # Muscle
                   "BUB3", "BUB1", "E2F1", "CDK1", "MKI67", "TOP2A", # Prolif
                   "NOTCH1", "NOTCH3"),
    upstream = 100000,
    downstream = 100000
)

plotPDF(plotList = p, 
        name = "MAST111-RMS-Plot-Tracks-Marker-Genes.pdf",
        ArchRProj = proj_tmp,
        addDOC = FALSE, width = 12, height = 6)


proj_tmp <- addGeneIntegrationMatrix(
    ArchRProj = proj_tmp, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = MAST111,
    addToArrow = TRUE,
    force= TRUE,
    groupList = groupList,
    groupRNA = "RNA_snn_res.0.8",
    nameCell = "predictedCell",
    nameGroup = "predictedGroup",
    nameScore = "predictedScore"
)

proj_tmp <- addImputeWeights(proj_tmp)

p1 <- plotEmbedding(
    ArchRProj = proj_tmp, 
    colorBy = "GeneIntegrationMatrix", 
    name = c("CD44", "EGFR", "MEOX2", "OGN", "POSTN", # EMT
             "MYOG", "MYL4", "TNNT3", "TNNT2", "MEF2C", "DES", # Muscle
             "BUB3", "BUB1", "E2F1", "CDK1", "MKI67", "TOP2A", # Prolif
             "NOTCH1", "NOTCH3"),
    continuousSet = "horizonExtra",
    embedding = "UMAP",
    imputeWeights = getImputeWeights(proj_tmp))

p2 <- plotEmbedding(
    ArchRProj = proj_tmp, 
    colorBy = "GeneScoreMatrix", 
    continuousSet = "horizonExtra",
    name = c("CD44", "EGFR", "MEOX2", "OGN", "POSTN", # EMT
             "MYOG", "MYL4", "TNNT3", "TNNT2", "MEF2C", "DES", # Muscle
             "BUB3", "BUB1", "E2F1", "CDK1", "MKI67", "TOP2A", # Prolif
             "NOTCH1", "NOTCH3"),
    embedding = "UMAP",
    imputeWeights = getImputeWeights(proj_tmp))


plotPDF(plotList = list(p1, p2), 
        name = "MAST111-Plot-UMAP-Marker-Genes-RNA-W-Imputation.pdf", 
        ArchRProj = proj_tmp, 
        addDOC = FALSE, width = 5, height = 5)

p4 <- plotEmbedding(ArchRProj = proj_tmp, colorBy = "cellColData", name = "predictedGroup", embedding = "UMAP")
plotPDF(p4, name = "MAST111_RMS_ATAC_Plot-UMAP-FinalClusters-RemoveC1-annotation-by-rna.pdf", ArchRProj = proj_tmp, addDOC = FALSE, width = 5, height = 5)


proj_tmp <- addGroupCoverages(ArchRProj = proj_tmp, groupBy = "predictedGroup")
## Use macs2
pathToMacs2 <- findMacs2()


proj_tmp <- addReproduciblePeakSet(
    ArchRProj = proj_tmp, 
    groupBy = "predictedGroup", 
    pathToMacs2 = pathToMacs2
)

getPeakSet(proj_tmp)

proj_tmp <- addPeakMatrix(proj_tmp)

getAvailableMatrices(proj_tmp)

table(proj_tmp$predictedGroup)


markersPeaks <- getMarkerFeatures(
    ArchRProj = proj_tmp,
    useMatrix = "PeakMatrix", 
    groupBy = "predictedGroup",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE)
markerList

heatmapPeaks <- markerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  transpose = TRUE
)
plotPDF(heatmapPeaks, name = "MAST111-Peak-Marker-Heatmap", width = 8, height = 6, ArchRProj = proj_tmp, addDOC = FALSE)


pma <- markerPlot(seMarker = markersPeaks, name = "EMT", cutOff = "FDR <= 0.1 & Log2FC >= 1", plotAs = "MA")
pv <- markerPlot(seMarker = markersPeaks, name = "EMT", cutOff = "FDR <= 0.1 & Log2FC >= 1", plotAs = "Volcano")
plotPDF(pma, pv, name = "MAST111-EMT-Markers-MA-Volcano", width = 5, height = 5, ArchRProj = proj_tmp, addDOC = FALSE)


p <- plotBrowserTrack(
    ArchRProj = proj_tmp, 
    groupBy = "predictedGroup", 
    geneSymbol = c("CD44", "EGFR", "MEOX2", "OGN", "POSTN", # EMT
             "MYOG", "MYL4", "TNNT3", "TNNT2", "MEF2C", "DES", # Muscle
             "BUB3", "BUB1", "E2F1", "CDK1", "MKI67", "TOP2A", # Prolif
             "NOTCH1", "NOTCH3"),
    features =  getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 1", returnGR = TRUE)["EMT"],
    upstream = 50000,
    downstream = 50000
)

plotPDF(p, name = "MAST111-EMT-Plot-Tracks-With-Features", width = 5, height = 5, ArchRProj = proj_tmp, addDOC = FALSE)


proj_tmp <- addMotifAnnotations(ArchRProj = proj_tmp,
                                motifSet = "cisbp", name = "Motif")

markerTest <- getMarkerFeatures(
  ArchRProj = proj_tmp,
  useMatrix = "PeakMatrix",
  groupBy = "predictedGroup",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "EMT",
  bgdGroups = "Ground"
)
pma <- markerPlot(seMarker = markerTest, name = "EMT", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "MA")
pv <- markerPlot(seMarker = markerTest, name = "EMT", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "Volcano")
plotPDF(pma, pv, name = "MAST111-EMT-vs-Ground-Markers-MA-Volcano", width = 5, height = 5, ArchRProj = proj_tmp, addDOC = FALSE)

## Pairwise motif test for EMT vs Ground
motifsUp <- peakAnnoEnrichment(
    seMarker = markerTest,
    ArchRProj = proj_tmp,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)
df <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))

## All comparisons
enrichMotif <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = proj_tmp,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)

enrichMotif@assays@data$mlog10Padj[grep('MEF', rownames(enrichMotif)),]
enrichMotif@assays@data$mlog10Padj[grep('RBPJ', rownames(enrichMotif)),]

## heatmapEM <- plotEnrichHeatmap(enrichMotif, n = 50, transpose = TRUE)-
heatmapEM <- plotEnrichHeatmap(enrichMotif, n=300, transpose = F) #cutOff=10,
plotPDF(heatmapEM, name = "MAST111-Motifs-Enriched-Marker-Heatmap", width = 9, height = 20, ArchRProj = proj_tmp, addDOC = FALSE)


## ENCODE enrichment

#proj_tmp <- addArchRAnnotations(ArchRProj = proj_tmp, collection = "EncodeTFBS")
#
#enrichEncode <- peakAnnoEnrichment(
#    seMarker = markersPeaks,
#    ArchRProj = proj_tmp,
#    peakAnnotation = "EncodeTFBS",
#    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
#  )
#
#heatmapEncode <- plotEnrichHeatmap(enrichEncode, n = 50, transpose = TRUE)
#plotPDF(heatmapEncode, name = "MAST111-EncodeTFBS-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = proj_tmp, addDOC = FALSE)

#proj_tmp <- addArchRAnnotations(ArchRProj = proj_tmp, collection = "ATAC")
#
#enrichATAC <- peakAnnoEnrichment(
#    seMarker = markersPeaks,
#    ArchRProj = proj_tmp,
#    peakAnnotation = "ATAC",
#    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
#)
### no enrichment at all
#heatmapATAC <- plotEnrichHeatmap(enrichATAC, n = 7, transpose = TRUE)
#plotPDF(heatmapATAC, name = "MAST111-ATAC-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = proj_tmp, addDOC = FALSE)


#proj_tmp <- addArchRAnnotations(ArchRProj = proj_tmp, collection = "Codex")
#
#enrichCodex <- peakAnnoEnrichment(
#    seMarker = markersPeaks,
#    ArchRProj = proj_tmp,
#    peakAnnotation = "Codex",
#    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
#  )
#
#heatmapCodex <- plotEnrichHeatmap(enrichCodex, n = 10, transpose = TRUE)
#plotPDF(heatmapCodex, name = "MAST111-Codex-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = proj_tmp, addDOC = FALSE)

## ChromVAR
proj_tmp <- addBgdPeaks(proj_tmp)

proj_tmp <- addDeviationsMatrix(
  ArchRProj = proj_tmp, 
  peakAnnotation = "Motif",
  force = TRUE
)

plotVarDev <- getVarDeviations(proj_tmp, name = "MotifMatrix", plot = TRUE)

plotPDF(plotVarDev, name = "MAST111-Variable-Motif-Deviation-Scores", width = 5, height = 5, ArchRProj = proj_tmp, addDOC = FALSE)


motifs = c("MYF5", "MYF6", "MYOG", "MYOD1", "MEF2C", "RBPJ", "CEBPB")

markerMotifs <- getFeatures(proj_tmp, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")

p <- plotGroups(ArchRProj = proj_tmp,
  groupBy = "predictedGroup", 
  colorBy = "MotifMatrix", 
  name = markerMotifs,
  imputeWeights = getImputeWeights(proj_tmp)
)

plotPDF(p, name = "MAST111-Plot-Groups-Deviations-w-Imputation", width = 5, height = 5, ArchRProj = proj_tmp, addDOC = FALSE)


p1 <- plotEmbedding(
    ArchRProj = proj_tmp, 
    colorBy = "MotifMatrix", 
    name = sort(markerMotifs), 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(proj_tmp)
)

markerRNA <- getFeatures(proj_tmp, select = paste(c("MYF5", "MYF6", "MYOG", "MYOD1", "MEF2C", "NOTCH1", "NOTCH3"), collapse="|"), useMatrix = "GeneScoreMatrix")
markerRNA

p2 <- plotEmbedding(
    ArchRProj = proj_tmp, 
    colorBy = "GeneScoreMatrix", 
    name = sort(markerRNA), 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(proj_tmp)
)

markerRNA <- getFeatures(proj_tmp, select = paste(c("MYF5", "MYF6", "MYOG", "MYOD1", "MEF2C", "NOTCH1", "NOTCH3"), collapse="|"), useMatrix = "GeneIntegrationMatrix")
markerRNA

p3 <- plotEmbedding(
    ArchRProj = proj_tmp,
    colorBy = "GeneIntegrationMatrix", 
    name = sort(markerRNA), 
    embedding = "UMAP",
    continuousSet = "blueYellow",
    imputeWeights = getImputeWeights(proj_tmp)
)

plotPDF(p1,p2,p3, name = "MAST111-Plot-Groups-Deviations-w-Imputation-on-UMAP", width = 5, height = 5, ArchRProj = proj_tmp, addDOC = FALSE)

## Footprinting

motifPositions <- getPositions(proj_tmp)

seFoot <- getFootprints(
  ArchRProj = proj_tmp, 
  positions = motifPositions[gsub(".*:", "", markerMotifs)],
  groupBy = "predictedGroup"
)

plot_fp1 = plotFootprints(
  seFoot = seFoot,
  ArchRProj = proj_tmp, 
  normMethod = "Subtract",
  plotName = "MAST111-Plot-Groups-Footprints-Subtract-Bias",
  addDOC = FALSE,
  smoothWindow = 5
)

plot_fp1 = plotFootprints(
  seFoot = seFoot,
  ArchRProj = proj_tmp, 
  normMethod = "Divide",
  plotName = "MAST111-Plot-Groups-Footprints-Divide-Bias",
  addDOC = FALSE,
  smoothWindow = 5
)
