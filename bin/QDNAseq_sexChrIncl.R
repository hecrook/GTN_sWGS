#!/usr/bin/env Rscript
# load required packages
    library(QDNAseq)
    #sbatch library(QDNASeq.hg19)
    library(Biobase)
    library(DNAcopy)
    library(CGHcall)
    library(future)
    library(tidyverse)

    future::plan("multisession", workers=4)

    # Get command line args
    args = commandArgs(trailingOnly = FALSE)
    print(args)

    # Get bin annotations.
    # Avialable bin sizes: 1, 5, 10, 15, 30, 50, 100, 500, and 1000 kbp. 
    binsize_str <- args[7]
    # Convert to numeric format as commandArgs creates a character, and getBinAnnnotations expects a numeric input
    binsize <- as.numeric(binsize_str)

    bins <- getBinAnnotations(binSize=binsize,genome="hg19")

    # Get bam files
    bamfiles = args[6]
    bamfiles = str_split(bamfiles, pattern = " ")[[1]]
    # bamfiles <- c(${1})
    print(bamfiles)

    # Load sequencing data
    cat("Reading counts per bin ..\n")
    readCounts <- binReadCounts(bins = bins, bamfiles = bamfiles, cache = T, pairedEnds = T, isDuplicate=F)
    save(readCounts, file="readCounts.RData")

    # Plot a raw copy number profile (read counts across the genome), and highlight bins that will be removed with default filtering
    cat("Plotting raw copy number profiles i.e. read counts across the genome..\n")
    cairo_pdf("RawCopyNumberHighlights.pdf",width = 10, height = 7,onefile = T)
    plot(readCounts, logTransform=FALSE, ylim=c(-50, 200))
    highlightFilters(readCounts, logTransform=FALSE,residual=TRUE, blacklist=TRUE)
    dev.off()

    # Apply filters and plot median read counts per bin as a function of GC content and mappability
    cat("Applying filters and plotting median read counts as a function of GC content and mappability..\n")
    readCountsFiltered <- applyFilters(readCounts,residual=TRUE, blacklist=TRUE)
    cairo_pdf("isobarPlot.pdf",width = 10, height = 7,onefile = T)
    isobarPlot(readCountsFiltered)
    dev.off()
    save(readCountsFiltered,file="readCountsFiltered.RData")

    # Estimate the correction for GC content and mappability, and make a plot for the relationship between the observed standard deviation in the data and its read depth
    cat("Plotting observed SD vs read depth..\n")
    readCountsFiltered <- estimateCorrection(readCountsFiltered)
    save(readCountsFiltered,file="readCountsFiltered_correctionEstimated.RData")
    cairo_pdf("noisePlot.pdf",width = 10, height = 7,onefile = T)
    noisePlot(readCountsFiltered)
    dev.off()

    # Re-apply the filters now that the LOESS correction has been calculated (without the sex chromosomes) to all chromosomes including the sex chromosomes.
    # Note that the GC and mappability correction was not calculated for the sex chromosomes specifically, and bins within this region may have a combination of GC content/ mappability found nowehre else in the genome, and will therefore be removed by correctBins.
    cat("Adding the sex chromsomes back in and plotting median read counts as a function of GC content and mappability. Note that the sex chromosomes were not included when calculating GC content and mappability correction..\n")
    readCountsFiltered <- applyFilters(readCountsFiltered, residual=TRUE, blacklist=TRUE, chromosomes=NA)
    cairo_pdf("isobarPlot_SEX.pdf",width = 10, height = 7,onefile = T)
    isobarPlot(readCountsFiltered)
    dev.off()
    save(readCountsFiltered,file="readCountsFiltered_correctionEstimated_sexApplied.RData")

    # apply the correction for GC content and mappability which we then normalize, smooth outliers, calculate segmentation and plot the copy number profile
    cat("GC corrections..\n")
    copyNumbers <- correctBins(readCountsFiltered)
    save(copyNumbers,file="copyNumbers.RData")
    cat("Normalisation..\n")
    copyNumbersNormalized <- normalizeBins(copyNumbers)
    save(copyNumbersNormalized,file="copyNumbersNormalized.RData")
    cat("Smoothing outliers..\n")
    copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)
    save(copyNumbersSmooth,file="copyNumbersSmooth.RData")
    cat("Plotting smooth CN profiles..\n")
    cairo_pdf("copyNumbersSmooth.pdf",width = 10, height = 7,onefile = T)
    plot(copyNumbersSmooth)
    dev.off()

    # Data is now ready to be analyzed with a downstream package of choice. For analysis with an external program or for visualizations in IGV, the data can be exported to a file
    cat("Exporting binned data in text format..\n")
    exportBins(copyNumbersSmooth, file="copyNumbersSmooth.txt")

    cat("CBS segmentation..\n")
    copyNumbersSegmented <- segmentBins(copyNumbersSmooth, transformFun="sqrt")
    save(copyNumbersSegmented,file="copyNumbersSegmented.RData")
    cat("Normalising CBS segmented data..\n")
    copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)
    save(copyNumbersSegmented,file="copyNumbersSegmented_normalised.RData")
    cat("Plotting normalised CBS segmented data..\n")
    cairo_pdf("copyNumbersSegmented.pdf",width = 10, height = 7,onefile = T)
    plot(copyNumbersSegmented)
    dev.off()
    cat("Calling CN per bins..\n")
    copyNumbersCalled <- callBins(copyNumbersSegmented)
    save(copyNumbersCalled,file="copyNumbersCalled.RData")
    cairo_pdf("copyNumbersCalled.pdf",width = 10, height = 7,onefile = T)
    plot(copyNumbersCalled)
    dev.off()


    copyNumbersCalled_segmented <-assayData(copyNumbersCalled)$segmented
    save(copyNumbersCalled_segmented,file="copyNumbersCalled_segmented.RData")