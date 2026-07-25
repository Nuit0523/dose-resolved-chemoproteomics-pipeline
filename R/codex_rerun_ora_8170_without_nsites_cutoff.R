suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ReactomePA)
  library(AnnotationDbi)
})

input_dir <- "E:/gradthesis/curvecurator/8170"
summary_file <- file.path(input_dir, "gene_level_pEC50_summary.csv")
mapping_file <- file.path(input_dir, "enrichment_compare_topN", "background_symbol_to_entrez.csv")
out_dir <- "C:/Users/owner/Documents/Codex/2026-07-07/ni-ha/outputs/ora_8170_no_nsites_cutoff"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

gene_summary <- read.csv(summary_file, check.names = FALSE, stringsAsFactors = FALSE)
symbol_to_entrez <- read.csv(mapping_file, check.names = FALSE, stringsAsFactors = FALSE)

required_cols <- c("gene_symbol", "best_pEC50", "n_sites")
missing_cols <- setdiff(required_cols, colnames(gene_summary))
if (length(missing_cols) > 0) {
  stop("Missing required columns in 8170 gene summary: ", paste(missing_cols, collapse = ", "))
}

symbol_to_entrez <- symbol_to_entrez[!is.na(symbol_to_entrez$SYMBOL) & !is.na(symbol_to_entrez$ENTREZID), ]
symbol_to_entrez <- symbol_to_entrez[!duplicated(symbol_to_entrez$SYMBOL), ]
background_entrez <- unique(as.character(symbol_to_entrez$ENTREZID))

# New candidate definition: keep potency cutoff, remove n_sites >= 2.
candidate <- gene_summary[!is.na(gene_summary$best_pEC50) & gene_summary$best_pEC50 >= 4, ]
candidate <- candidate[order(-candidate$best_pEC50, candidate$gene_symbol), ]

write.csv(candidate, file.path(out_dir, "candidate_genes_best_pEC50_ge4_no_nsites_cutoff.csv"), row.names = FALSE)

writeLines(c(
  "8170 ORA rerun without n_sites cutoff",
  paste0("Input gene-level summaries: ", nrow(gene_summary)),
  paste0("Candidate definition: best_pEC50 >= 4; no n_sites cutoff"),
  paste0("Candidate genes: ", nrow(candidate)),
  paste0("Ranking: best_pEC50 descending"),
  paste0("Background Entrez IDs: ", length(background_entrez)),
  paste0("TopN tested: ", paste(c(30, 50, 80), collapse = ", "))
), file.path(out_dir, "run_parameters.txt"))

write_result <- function(obj, file) {
  df <- as.data.frame(obj)
  write.csv(df, file, row.names = FALSE)
  invisible(df)
}

run_one_topn <- function(n) {
  top_dir <- file.path(out_dir, paste0("top", n))
  dir.create(top_dir, recursive = TRUE, showWarnings = FALSE)

  top_df <- head(candidate, n)
  top_map <- merge(
    data.frame(SYMBOL = top_df$gene_symbol, stringsAsFactors = FALSE),
    symbol_to_entrez,
    by = "SYMBOL",
    all.x = TRUE
  )
  top_entrez <- unique(as.character(top_map$ENTREZID[!is.na(top_map$ENTREZID)]))

  write.csv(top_df, file.path(top_dir, paste0("hit_top", n, "_table.csv")), row.names = FALSE)
  write.csv(data.frame(gene_symbol = top_df$gene_symbol), file.path(top_dir, paste0("hit_top", n, "_genes.csv")), row.names = FALSE)
  write.csv(top_map, file.path(top_dir, paste0("hit_top", n, "_symbol_to_entrez.csv")), row.names = FALSE)

  ego <- enrichGO(
    gene = top_entrez,
    universe = background_entrez,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE
  )

  ekegg <- tryCatch({
    kk <- enrichKEGG(
      gene = top_entrez,
      universe = background_entrez,
      organism = "hsa",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2
    )
    if (!is.null(kk) && nrow(as.data.frame(kk)) > 0) {
      kk <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
    }
    kk
  }, error = function(e) {
    writeLines(conditionMessage(e), file.path(top_dir, paste0("KEGG_top", n, "_error.txt")))
    NULL
  })

  ereact <- enrichPathway(
    gene = top_entrez,
    universe = background_entrez,
    organism = "human",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE
  )

  go_df <- write_result(ego, file.path(top_dir, paste0("GO_BP_top", n, ".csv")))
  kegg_df <- if (is.null(ekegg)) {
    data.frame()
  } else {
    write_result(ekegg, file.path(top_dir, paste0("KEGG_top", n, ".csv")))
  }
  react_df <- write_result(ereact, file.path(top_dir, paste0("Reactome_top", n, ".csv")))

  writeLines(c(
    paste0("TopN: ", n),
    paste0("Input hit genes: ", nrow(top_df)),
    paste0("Mapped hit genes: ", length(top_entrez)),
    paste0("Background mapped genes: ", length(background_entrez)),
    paste0("GO BP terms: ", nrow(go_df)),
    paste0("KEGG pathways: ", nrow(kegg_df)),
    paste0("Reactome pathways: ", nrow(react_df))
  ), file.path(top_dir, paste0("summary_top", n, ".txt")))

  data.frame(
    TopN = n,
    input_genes = nrow(top_df),
    mapped_entrez = length(top_entrez),
    GO_BP_terms = nrow(go_df),
    KEGG_terms = nrow(kegg_df),
    Reactome_terms = nrow(react_df),
    stringsAsFactors = FALSE
  )
}

summary_df <- do.call(rbind, lapply(c(30, 50, 80), run_one_topn))
write.csv(summary_df, file.path(out_dir, "enrichment_summary_no_nsites_cutoff.csv"), row.names = FALSE)
print(summary_df)
