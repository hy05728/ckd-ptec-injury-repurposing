suppressMessages({
  library(clusterProfiler); library(org.Hs.eg.db)
  library(enrichplot); library(ggplot2)
})
setwd("E:/workbuddy/2026-11-08-00-03-12")
out <- "results/phase3"
dir.create(out, showWarnings = FALSE)

sig <- read.csv("results/vscore_injury_normal2injured.csv", check.names = FALSE, stringsAsFactors = FALSE)
colnames(sig)[1:2] <- c("gene", "vscore")
sig$gene <- gsub("-\\d+$", "", sig$gene)               # 去掉 make_unique 后缀
sig <- sig[!duplicated(sig$gene) & sig$gene != "" & !is.na(sig$vscore), ]
cat("genes:", nrow(sig), "| up:", sum(sig$vscore > 0), "| down:", sum(sig$vscore < 0), "\n")

sig <- sig[order(sig$vscore, decreasing = TRUE), ]
up   <- head(sig$gene[sig$vscore > 0], 500)            # 损伤上调 top500
down <- tail(sig$gene[sig$vscore < 0], 500)            # 损伤下调 top500

run_go <- function(genes, tag) {
  eg <- enrichGO(genes, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP",
                 pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
  if (is.null(eg) || nrow(as.data.frame(eg)) == 0) { cat("[GO]", tag, "no hit\n"); return(invisible(NULL)) }
  df <- as.data.frame(eg)
  write.csv(df, file.path(out, paste0("enrich_GO_BP_", tag, ".csv")), row.names = FALSE)
  cat("[GO]", tag, "terms:", nrow(df), "| top:", df$Description[1], "(padj", signif(df$p.adjust[1],2), ")\n")
  p <- dotplot(eg, showCategory = 15) + ggtitle(paste("GO BP -", tag))
  ggsave(file.path(out, paste0("enrich_GO_BP_", tag, ".png")), p, width = 9, height = 6, dpi = 150)
  invisible(df)
}

run_kegg <- function(genes, tag) {
  ent <- tryCatch(bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID,
                  error = function(e) NULL)
  if (is.null(ent) || length(ent) == 0) { cat("[KEGG]", tag, "no entrez\n"); return(invisible(NULL)) }
  eg <- tryCatch(enrichKEGG(ent, organism = "hsa", keyType = "kegg",
                            pvalueCutoff = 0.05, pAdjustMethod = "BH"),
                 error = function(e) { cat("[KEGG]", tag, "failed:", conditionMessage(e), "\n"); NULL })
  if (is.null(eg) || nrow(as.data.frame(eg)) == 0) { cat("[KEGG]", tag, "no hit\n"); return(invisible(NULL)) }
  df <- setReadable(eg, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")
  df <- as.data.frame(df)
  write.csv(df, file.path(out, paste0("enrich_KEGG_", tag, ".csv")), row.names = FALSE)
  cat("[KEGG]", tag, "pathways:", nrow(df), "| top:", df$Description[1], "(padj", signif(df$p.adjust[1],2), ")\n")
  invisible(df)
}

cat("=== UP (损伤上调) ===\n");  d_up <- run_go(up, "up");  k_up <- run_kegg(up, "up")
cat("=== DOWN (损伤下调) ===\n"); d_dn <- run_go(down, "down"); k_dn <- run_kegg(down, "down")

# 写一份人读摘要
sink(file.path(out, "enrichment_signature_readable.txt"))
cat("损伤签名通路富集（GSE131882 healthy PT -> injured PT, v-score）\n\n")
if (!is.null(d_up)) { cat("## GO BP 损伤上调 top15\n"); print(head(d_up[,c("Description","GeneRatio","p.adjust")],15)) }
if (!is.null(k_up)) { cat("\n## KEGG 损伤上调 top15\n"); print(head(k_up[,c("Description","GeneRatio","p.adjust")],15)) }
if (!is.null(d_dn)) { cat("\n## GO BP 损伤下调 top15\n"); print(head(d_dn[,c("Description","GeneRatio","p.adjust")],15)) }
sink()
cat("[done] enrichment written\n")
