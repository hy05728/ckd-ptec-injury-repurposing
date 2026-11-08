suppressMessages(library(Matrix))
suppressMessages(library(AnnotationDbi))
suppressMessages(library(org.Hs.eg.db))

raw_dir <- "E:/workbuddy/2026-11-08-00-03-12/data/raw"
out_dir <- "E:/workbuddy/2026-11-08-00-03-12/data/processed"
dir.create(out_dir, showWarnings = FALSE)

files <- sort(list.files(raw_dir, pattern = "\\.dgecounts\\.rds\\.gz$", full.names = TRUE))
cat("found", length(files), "files\n")

# Ensembl -> Symbol 映射（用第一个样本的基因集）
first <- readRDS(gzcon(gzfile(files[1])))
ens <- rownames(first$umicount$exon$all)
rm(first); gc()
sym <- select(org.Hs.eg.db, keys = ens, columns = "SYMBOL", keytype = "ENSEMBL")
sym <- sym[!is.na(sym$SYMBOL) & !duplicated(sym$ENSEMBL), ]
map <- setNames(sym$SYMBOL, sym$ENSEMBL)
symbols <- ifelse(ens %in% names(map), unname(map[ens]), ens)
cat("mapped", sum(ens %in% names(map)), "/", length(ens), "genes to symbol\n")

for (f in files) {
  sid <- sub("\\.dgecounts\\.rds\\.gz$", "", basename(f))   # GSMxxxx_control.s1
  cat("processing", sid, "...\n")
  x <- readRDS(gzcon(gzfile(f)))
  m <- x$umicount$exon$all
  rm(x); gc()
  idx <- match(rownames(m), ens)
  gsym <- symbols[idx]
  d <- file.path(out_dir, sid); dir.create(d, showWarnings = FALSE)
  writeMM(m, file.path(d, "matrix.mtx"))
  write.table(data.frame(ens = rownames(m), symbol = gsym),
              file.path(d, "genes.tsv"), sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)
  writeLines(colnames(m), file.path(d, "barcodes.tsv"))
  cat("  wrote", dim(m)[1], "genes x", dim(m)[2], "cells\n")
}
cat("DONE\n")
