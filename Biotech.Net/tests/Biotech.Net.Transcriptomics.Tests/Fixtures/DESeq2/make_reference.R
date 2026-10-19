# Generates the DESeq2 / apeglm reference fixtures used by DESeq2ParityTests.
#
#   Rscript make_reference.R <pasilla_gene_counts.tsv> <pasilla_sample_annotation.csv> <outdir> [subset|full] [airway.RData]
#
# The pasilla input files ship with the Bioconductor data package "pasilla" (inst/extdata); airway.RData with the
# data package "airway" (data/). In "subset" mode (the
# committed fixtures, written as gzip-compressed CSV) every 7th gene of the file is kept; "full" mode runs every gene and is used for the
# one-off full-dataset comparison reported in the docs. Reference versions: R 4.6.1, DESeq2 1.52.0,
# apeglm 1.34.0 (Bioconductor 3.23).
suppressPackageStartupMessages({ library(DESeq2); library(apeglm) })
args <- commandArgs(trailingOnly = TRUE)
ctsFile <- args[1]; annoFile <- args[2]; out <- args[3]; mode <- if (length(args) > 3) args[4] else "subset"
airwayFile <- if (length(args) > 4) args[5] else NA
dir.create(out, showWarnings = FALSE, recursive = TRUE)

cts <- as.matrix(read.csv(ctsFile, sep = "\t", row.names = "gene_id"))
coldata <- read.csv(annoFile, row.names = 1)[, c("condition", "type")]
rownames(coldata) <- sub("fb", "", rownames(coldata))
coldata$condition <- factor(coldata$condition, levels = c("untreated", "treated"))
coldata$type <- factor(coldata$type)
cts <- cts[, rownames(coldata)]
if (mode == "subset") cts <- cts[seq(1, nrow(cts), by = 7), ]
storage.mode(cts) <- "integer"

fmt <- function(x) {
  if (is.logical(x)) return(ifelse(is.na(x), "NA", ifelse(x, "1", "0")))
  ifelse(is.na(x), "NA", formatC(x, digits = 13, format = "g"))
}
writeTable <- function(df, file) {
  df <- as.data.frame(df)
  outdf <- data.frame(gene = rownames(cts), stringsAsFactors = FALSE)
  for (n in names(df)) outdf[[n]] <- fmt(df[[n]])
  con <- gzfile(file.path(out, paste0(file, ".gz")), "w")
  write.table(outdf, con, sep = ",", quote = FALSE, row.names = FALSE)
  close(con)
}
scalars <- list()
addScalar <- function(scenario, key, value) {
  scalars[[length(scalars) + 1]] <<- data.frame(scenario = scenario, key = key, value = fmt(value))
}
geneCols <- c("baseMean", "baseVar", "allZero", "dispGeneEst", "dispGeneIter", "dispFit", "dispersion",
              "dispIter", "dispOutlier", "dispMAP", "betaConv", "betaIter", "deviance", "maxCooks")
dump <- function(scenario, dds, res, extra = list()) {
  mc <- mcols(dds)
  cols <- intersect(c(geneCols, "replace", resultsNames(dds), paste0("SE_", resultsNames(dds)),
                      "LRTStatistic", "LRTPvalue"), names(mc))
  df <- as.list(mc[, cols])
  df$res_baseMean <- res$baseMean; df$res_log2FoldChange <- res$log2FoldChange; df$res_lfcSE <- res$lfcSE
  df$res_stat <- res$stat; df$res_pvalue <- res$pvalue; df$res_padj <- res$padj
  for (n in names(extra)) df[[n]] <- extra[[n]]
  writeTable(df, paste0(scenario, ".csv"))
  sf <- sizeFactors(dds)
  for (j in seq_along(sf)) addScalar(scenario, paste0("sizeFactor_", j), sf[j])
  fn <- dispersionFunction(dds)
  co <- attr(fn, "coefficients")
  if (!is.null(co)) { addScalar(scenario, "asymptDisp", co[1]); addScalar(scenario, "extraPois", co[2]) }
  if (!is.null(attr(fn, "mean"))) addScalar(scenario, "meanDisp", attr(fn, "mean"))
  addScalar(scenario, "varLogDispEsts", attr(fn, "varLogDispEsts"))
  addScalar(scenario, "dispPriorVar", attr(fn, "dispPriorVar"))
  md <- metadata(res)
  if (!is.null(md$filterThreshold)) {
    addScalar(scenario, "filterThreshold", unname(md$filterThreshold))
    addScalar(scenario, "filterTheta", md$filterTheta)
  }
}

# A: ~ condition (defaults; linear-mu path for gene-wise dispersions)
ddsA <- DESeq(DESeqDataSetFromMatrix(cts, coldata, ~ condition), quiet = TRUE)
resA <- results(ddsA)
apeA <- lfcShrink(ddsA, coef = 2, type = "apeglm", quiet = TRUE, returnList = TRUE)
apeAs <- lfcShrink(ddsA, coef = 2, type = "apeglm", quiet = TRUE, svalue = TRUE)
norA <- lfcShrink(ddsA, coef = 2, type = "normal", quiet = TRUE)
dump("A_condition", ddsA, resA, list(ape_lfc = apeA$res$log2FoldChange, ape_se = apeA$res$lfcSE,
     ape_svalue = apeAs$svalue, ape_conv = apeA$fit$diag[, "conv"],
     normal_lfc = norA$log2FoldChange, normal_se = norA$lfcSE))
addScalar("A_condition", "ape_prior_scale", apeA$fit$prior.control$prior.scale)
addScalar("A_condition", "normal_betaPriorVar", priorInfo(norA)$betaPriorVar[2])
resA0 <- results(ddsA, independentFiltering = FALSE, cooksCutoff = FALSE)
writeTable(data.frame(pvalue = resA0$pvalue, padj = resA0$padj), "A_condition_nofilter.csv")

# B: ~ type + condition (GLM-mu path; contrasts)
ddsB <- DESeq(DESeqDataSetFromMatrix(cts, coldata, ~ type + condition), quiet = TRUE)
resB <- results(ddsB)
resBflip <- results(ddsB, contrast = c("condition", "untreated", "treated"))
resBnum <- results(ddsB, contrast = c(0, -1, 1))
apeB <- lfcShrink(ddsB, coef = 3, type = "apeglm", quiet = TRUE, returnList = TRUE)
norB <- lfcShrink(ddsB, coef = 3, type = "normal", quiet = TRUE)
dump("B_type_condition", ddsB, resB, list(flip_lfc = resBflip$log2FoldChange, flip_stat = resBflip$stat,
     flip_padj = resBflip$padj, num_lfc = resBnum$log2FoldChange, num_se = resBnum$lfcSE,
     num_pvalue = resBnum$pvalue, num_padj = resBnum$padj,
     ape_lfc = apeB$res$log2FoldChange, ape_se = apeB$res$lfcSE, ape_conv = apeB$fit$diag[, "conv"],
     normal_lfc = norB$log2FoldChange, normal_se = norB$lfcSE))
addScalar("B_type_condition", "ape_prior_scale", apeB$fit$prior.control$prior.scale)
addScalar("B_type_condition", "normal_betaPriorVar_type", priorInfo(norB)$betaPriorVar[2])
addScalar("B_type_condition", "normal_betaPriorVar_condition", priorInfo(norB)$betaPriorVar[3])

# C: LRT ~ type + condition vs ~ type
ddsC <- DESeq(DESeqDataSetFromMatrix(cts, coldata, ~ type + condition), test = "LRT", reduced = ~ type, quiet = TRUE)
dump("C_lrt_type", ddsC, results(ddsC))

# D: LRT ~ condition vs ~ 1 (intercept-only shortcut for the reduced model)
ddsD <- DESeq(DESeqDataSetFromMatrix(cts, coldata, ~ condition), test = "LRT", reduced = ~ 1, quiet = TRUE)
dump("D_lrt_intercept", ddsD, results(ddsD))

# E: poscounts size factors
ddsE <- DESeq(DESeqDataSetFromMatrix(cts, coldata, ~ condition), sfType = "poscounts", quiet = TRUE)
dump("E_poscounts", ddsE, results(ddsE))

# F: mean trend
ddsF <- DESeq(DESeqDataSetFromMatrix(cts, coldata, ~ condition), fitType = "mean", quiet = TRUE)
dump("F_mean", ddsF, results(ddsF))

# G: local (locfit) trend
ddsG <- DESeq(DESeqDataSetFromMatrix(cts, coldata, ~ condition), fitType = "local", quiet = TRUE)
dump("G_local", ddsG, results(ddsG))

# H: Cook's outlier replacement (minReplicatesForReplace = 3 so pasilla's 3 + 4 replicates qualify)
ddsH <- DESeq(DESeqDataSetFromMatrix(cts, coldata, ~ condition), minReplicatesForReplace = 3, quiet = TRUE)
dump("H_replace", ddsH, results(ddsH))

write.table(data.frame(gene_id = rownames(cts), cts, check.names = FALSE),
            file.path(out, if (mode == "subset") "pasilla_subset_counts.tsv" else "pasilla_full_counts.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(sample = rownames(coldata), condition = coldata$condition, type = coldata$type),
            file.path(out, "pasilla_coldata.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(c(R.version.string, paste("DESeq2", packageVersion("DESeq2")), paste("apeglm", packageVersion("apeglm"))),
           file.path(out, "versions.txt"))

# I: airway (8 samples, ~ cell + dex: 3 residual degrees of freedom, so DESeq2 estimates the dispersion prior variance by
# its seeded simulation). Subset mode keeps every 30th gene.
if (!is.na(airwayFile)) {
  suppressPackageStartupMessages(library(SummarizedExperiment))
  load(airwayFile)
  cts <- assay(airway)
  if (mode == "subset") cts <- cts[seq(1, nrow(cts), by = 30), ]
  storage.mode(cts) <- "integer"
  acd <- data.frame(cell = factor(airway$cell), dex = relevel(factor(airway$dex), "untrt"), row.names = colnames(cts))
  ddsI <- DESeq(DESeqDataSetFromMatrix(cts, acd, ~ cell + dex), quiet = TRUE)
  resI <- results(ddsI)
  apeI <- lfcShrink(ddsI, coef = "dex_trt_vs_untrt", type = "apeglm", quiet = TRUE)
  dump("I_airway", ddsI, resI, list(ape_lfc = apeI$log2FoldChange, ape_se = apeI$lfcSE))
  write.table(data.frame(gene_id = rownames(cts), cts, check.names = FALSE),
              file.path(out, if (mode == "subset") "airway_subset_counts.tsv" else "airway_full_counts.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(data.frame(sample = rownames(acd), cell = acd$cell, dex = acd$dex),
              file.path(out, "airway_coldata.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}

# Building blocks of the simulation: R's RNG stream and loess(surface = "interpolate").
set.seed(2)
rng <- c(runif(3), rnorm(3), rchisq(3, 1), rchisq(3, 2), rchisq(3, 3), rnorm(2, 0, 0), rgamma(3, 0.5), rnorm(2))
lx <- seq(0, 8, length.out = 200); ly <- sin(lx) + 0.1 * lx^2 + cos(7 * lx) / 5
lp <- predict(loess(ly ~ lx, span = 0.2), seq(0, 8, length.out = 1000))
writeLines(c(formatC(rng, digits = 17, format = "g")), file.path(out, "rng.txt"))
writeLines(c(formatC(lp, digits = 17, format = "g")), file.path(out, "loess.txt"))

write.table(do.call(rbind, scalars), file.path(out, "scalars.csv"), sep = ",", quote = FALSE, row.names = FALSE)
