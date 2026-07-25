# =========================================================
# Chemoproteomics enrichment pipeline
# compare top30 / top50 / top80
# =========================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

pkgs <- c(
  "readr", "dplyr", "stringr", "ggplot2",
  "clusterProfiler", "org.Hs.eg.db", "ReactomePA", "enrichplot"
)

for (pkg in pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

library(readr)
library(dplyr)
library(stringr)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(enrichplot)

# =========================================================
# 1. file paths
# =========================================================
input_file <- "E:/python/curvecurator/8170/gene_level_pEC50_summary.csv"
out_dir <- "E:/python/curvecurator/8170/enrichment_compare_topN"

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# =========================================================
# 2. read data
# =========================================================
df <- read_csv(input_file, show_col_types = FALSE)

cat("Input columns:\n")
print(colnames(df))

required_cols <- c("gene_symbol", "best_pEC50", "n_sites")
missing_cols <- setdiff(required_cols, colnames(df))

if (length(missing_cols) > 0) {
  stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
}

# =========================================================
# 3. clean data
# =========================================================
df_clean <- df %>%
  filter(!is.na(gene_symbol), gene_symbol != "") %>%
  mutate(gene_symbol = str_trim(gene_symbol)) %>%
  distinct(gene_symbol, .keep_all = TRUE)

cat("Total unique genes after cleaning:", nrow(df_clean), "\n")

# =========================================================
# 4. define background and candidate hits
# =========================================================
background_genes <- df_clean %>%
  pull(gene_symbol)

candidate_hits <- df_clean %>%
  filter(n_sites >= 2, best_pEC50 >= 4) %>%
  arrange(desc(best_pEC50), desc(n_sites))

cat("Candidate hits after cutoff:", nrow(candidate_hits), "\n")

# save background and candidate hits
write.csv(
  data.frame(gene_symbol = background_genes),
  file.path(out_dir, "background_genes.csv"),
  row.names = FALSE
)

write.csv(
  candidate_hits,
  file.path(out_dir, "candidate_hits_after_cutoff.csv"),
  row.names = FALSE
)

# =========================================================
# 5. helper function: convert SYMBOL to ENTREZ
# =========================================================
convert_to_entrez <- function(genes) {
  bitr(
    genes,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
}

bg_entrez <- convert_to_entrez(background_genes)

write.csv(
  bg_entrez,
  file.path(out_dir, "background_symbol_to_entrez.csv"),
  row.names = FALSE
)

cat("Background mapped to ENTREZ:", nrow(bg_entrez), "\n")

# =========================================================
# 6. helper functions for enrichment
# =========================================================
run_go <- function(hit_entrez, bg_entrez, p_cut = 0.2, q_cut = 0.5) {
  enrichGO(
    gene          = hit_entrez$ENTREZID,
    universe      = bg_entrez$ENTREZID,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENTREZID",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = p_cut,
    qvalueCutoff  = q_cut,
    readable      = TRUE
  )
}

run_kegg <- function(hit_entrez, bg_entrez, p_cut = 0.2) {
  kk <- enrichKEGG(
    gene          = hit_entrez$ENTREZID,
    universe      = bg_entrez$ENTREZID,
    organism      = "hsa",
    pvalueCutoff  = p_cut,
    pAdjustMethod = "BH"
  )
  
  kk_df <- as.data.frame(kk)
  if (nrow(kk_df) > 0) {
    kk <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  }
  kk
}

run_reactome <- function(hit_entrez, bg_entrez, p_cut = 0.2) {
  enrichPathway(
    gene          = hit_entrez$ENTREZID,
    universe      = bg_entrez$ENTREZID,
    organism      = "human",
    pvalueCutoff  = p_cut,
    pAdjustMethod = "BH",
    readable      = TRUE
  )
}

save_dotplot <- function(enrich_obj, file_prefix, out_dir, show_n = 15) {
  df <- as.data.frame(enrich_obj)
  
  if (nrow(df) == 0) {
    cat(file_prefix, ": no result, skip plot.\n")
    return(NULL)
  }
  
  pdf(file.path(out_dir, paste0(file_prefix, "_dotplot.pdf")), width = 10, height = 8)
  print(dotplot(enrich_obj, showCategory = min(show_n, nrow(df))))
  dev.off()
}

# =========================================================
# 7. compare different hit sizes
# =========================================================
top_n_list <- c(30, 50, 80)

summary_list <- list()

for (top_n in top_n_list) {
  
  cat("\n============================\n")
  cat("Running enrichment for top", top_n, "\n")
  cat("============================\n")
  
  sub_dir <- file.path(out_dir, paste0("top", top_n))
  if (!dir.exists(sub_dir)) {
    dir.create(sub_dir, recursive = TRUE)
  }
  
  # select hits
  hit_df <- candidate_hits %>%
    head(top_n)
  
  hit_genes <- hit_df$gene_symbol
  
  write.csv(
    hit_df,
    file.path(sub_dir, paste0("hit_top", top_n, "_table.csv")),
    row.names = FALSE
  )
  
  write.csv(
    data.frame(gene_symbol = hit_genes),
    file.path(sub_dir, paste0("hit_top", top_n, "_genes.csv")),
    row.names = FALSE
  )
  
  # map
  hit_entrez <- convert_to_entrez(hit_genes)
  
  write.csv(
    hit_entrez,
    file.path(sub_dir, paste0("hit_top", top_n, "_symbol_to_entrez.csv")),
    row.names = FALSE
  )
  
  cat("Input hit genes:", length(hit_genes), "\n")
  cat("Mapped hit genes:", nrow(hit_entrez), "\n")
  
  # enrichment
  ego <- run_go(hit_entrez, bg_entrez, p_cut = 0.2, q_cut = 0.5)
  ekegg <- run_kegg(hit_entrez, bg_entrez, p_cut = 0.2)
  ereact <- run_reactome(hit_entrez, bg_entrez, p_cut = 0.2)
  
  ego_df <- as.data.frame(ego)
  ekegg_df <- as.data.frame(ekegg)
  ereact_df <- as.data.frame(ereact)
  
  # save tables
  write.csv(ego_df, file.path(sub_dir, paste0("GO_BP_top", top_n, ".csv")), row.names = FALSE)
  write.csv(ekegg_df, file.path(sub_dir, paste0("KEGG_top", top_n, ".csv")), row.names = FALSE)
  write.csv(ereact_df, file.path(sub_dir, paste0("Reactome_top", top_n, ".csv")), row.names = FALSE)
  
  # save plots
  save_dotplot(ego, paste0("GO_BP_top", top_n), sub_dir)
  save_dotplot(ekegg, paste0("KEGG_top", top_n), sub_dir)
  save_dotplot(ereact, paste0("Reactome_top", top_n), sub_dir)
  
  # summary row
  top_go_desc <- if (nrow(ego_df) > 0) ego_df$Description[1] else NA
  top_go_padj <- if (nrow(ego_df) > 0) ego_df$p.adjust[1] else NA
  
  top_kegg_desc <- if (nrow(ekegg_df) > 0) ekegg_df$Description[1] else NA
  top_kegg_padj <- if (nrow(ekegg_df) > 0) ekegg_df$p.adjust[1] else NA
  
  top_react_desc <- if (nrow(ereact_df) > 0) ereact_df$Description[1] else NA
  top_react_padj <- if (nrow(ereact_df) > 0) ereact_df$p.adjust[1] else NA
  
  summary_list[[as.character(top_n)]] <- data.frame(
    top_n = top_n,
    input_hit_genes = length(hit_genes),
    mapped_hit_genes = nrow(hit_entrez),
    GO_terms = nrow(ego_df),
    KEGG_terms = nrow(ekegg_df),
    Reactome_terms = nrow(ereact_df),
    top_GO = top_go_desc,
    top_GO_padj = top_go_padj,
    top_KEGG = top_kegg_desc,
    top_KEGG_padj = top_kegg_padj,
    top_Reactome = top_react_desc,
    top_Reactome_padj = top_react_padj,
    stringsAsFactors = FALSE
  )
  
  # write text summary per group
  sink(file.path(sub_dir, paste0("summary_top", top_n, ".txt")))
  
  cat("====================================\n")
  cat("Enrichment summary for top", top_n, "\n")
  cat("====================================\n\n")
  
  cat("Input hit genes:", length(hit_genes), "\n")
  cat("Mapped hit genes:", nrow(hit_entrez), "\n")
  cat("Background mapped genes:", nrow(bg_entrez), "\n\n")
  
  cat("GO BP terms:", nrow(ego_df), "\n")
  cat("KEGG pathways:", nrow(ekegg_df), "\n")
  cat("Reactome pathways:", nrow(ereact_df), "\n\n")
  
  if (nrow(ego_df) > 0) {
    cat("Top GO terms:\n")
    print(head(ego_df[, c("Description", "Count", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "geneID")], 10))
    cat("\n")
  }
  
  if (nrow(ekegg_df) > 0) {
    cat("Top KEGG pathways:\n")
    print(head(ekegg_df[, c("Description", "Count", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "geneID")], 10))
    cat("\n")
  }
  
  if (nrow(ereact_df) > 0) {
    cat("Top Reactome pathways:\n")
    print(head(ereact_df[, c("Description", "Count", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "geneID")], 10))
    cat("\n")
  }
  
  sink()
}

# =========================================================
# 8. merge global summary
# =========================================================
summary_df <- bind_rows(summary_list)

write.csv(
  summary_df,
  file.path(out_dir, "enrichment_comparison_summary.csv"),
  row.names = FALSE
)

cat("\n====================================\n")
cat("Final comparison summary\n")
cat("====================================\n")
print(summary_df)




input_file <- "E:/python/curvecurator/8170/hit_top80_table_with_description.csv"
out_dir <- "E:/python/curvecurator/8170/protein_classification_final"

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

df <- read.csv(input_file, stringsAsFactors = FALSE)

if (!"gene_symbol" %in% colnames(df)) stop("missing gene_symbol")
if (!"prot_description" %in% colnames(df)) stop("missing prot_description")

classify_protein_function <- function(gene, desc) {
  txt <- tolower(paste(gene, desc))
  gene2 <- toupper(gene)
  
  if (grepl("splice|splicing|rna helicase|rna binding|ribonucleoprotein|mrna processing|mrna splicing|rna export|rrna processing|ribosome biogenesis|exosome", txt) ||
      grepl("^DDX|^DHX|^SRSF|^HNRNP|^PRPF|^LSM|^CPSF|^CSTF", gene2)) {
    return("RNA processing / splicing")
  }
  
  if (grepl("transcription|chromatin|histone|epigenetic|corepressor|coactivator|mediator|sin3|transcription regulator", txt) ||
      grepl("^SIN3|^KMT|^KDM|^BRD|^SMAR|^ARID|^CHD|^EP300|^CREBBP", gene2)) {
    return("Transcription / chromatin")
  }
  
  if (grepl("dna repair|repair|replication|topoisomerase|checkpoint|genome stability|dna damage", txt) ||
      grepl("^PRKDC|^TOP1|^TOP2|^RAD|^ATM|^ATR|^XRCC|^BRCA|^PARP", gene2)) {
    return("DNA repair / replication")
  }
  
  if (grepl("ribosomal|ribosome|translation|trna synthetase|initiation factor|elongation factor", txt) ||
      grepl("^RPL|^RPS|^EIF|^EEF|^AARS|^TARS|^GARS|^QARS|^MARS|^WARS|^YARS", gene2)) {
    return("Translation / ribosome")
  }
  
  if (grepl("ubiquitin|ubiquitin ligase|proteasome|chaperone|heat shock|folding|deubiquitinase", txt) ||
      grepl("^HUWE1|^UBE|^USP|^PSM|^HSP|^DNAJ|^BAG", gene2)) {
    return("Proteostasis / ubiquitin")
  }
  
  if (grepl("kinase|phosphatase|signal transduction|signaling|gtpase|receptor", txt) ||
      grepl("^MAPK|^AKT|^PRK|^PTPN|^PPP|^SRC|^JAK|^RAF|^RHO|^RAB|^GNA", gene2)) {
    return("Signaling")
  }
  
  if (grepl("metabolic|metabolism|biosynthesis|synthase|dehydrogenase|oxidase|reductase|transferase|lipid|sphingolipid|hexosamine", txt) ||
      grepl("^GFPT|^SGPL1|^ALDO|^GPI|^PKM|^LDH|^IDH|^MDH|^ACLY|^FASN", gene2)) {
    return("Metabolism")
  }
  
  if (grepl("cytoskeleton|actin|tubulin|microtubule|vesicle|trafficking|endosome|golgi|motor protein", txt) ||
      grepl("^ACT|^TUB|^KIF|^DYNC|^MYH|^MYO|^CLTC|^COP", gene2)) {
    return("Cytoskeleton / trafficking")
  }
  
  if (grepl("nuclear pore|nucleoporin|nuclear transport|importin|exportin|nucleocytoplasmic transport", txt) ||
      grepl("^NUP|^KPN|^XPO", gene2)) {
    return("Nuclear transport / pore")
  }
  
  if (grepl("nucleic acid binding|dna binding|rna binding", txt)) {
    return("Nucleic acid binding (general)")
  }
  
  return("Other / unclassified")
}

df$Functional_Class <- mapply(classify_protein_function, df$gene_symbol, df$prot_description)

write.csv(df, file.path(out_dir, "top80_functional_classification.csv"), row.names = FALSE)

tab <- as.data.frame(table(df$Functional_Class), stringsAsFactors = FALSE)
colnames(tab) <- c("Functional_Class", "Count")
tab$Percent <- round(tab$Count / sum(tab$Count) * 100, 2)

write.csv(tab, file.path(out_dir, "top80_functional_class_summary.csv"), row.names = FALSE)

print(tab)
cat("classification file exists: ", file.exists(file.path(out_dir, "top80_functional_classification.csv")), "\n")
cat("summary file exists: ", file.exists(file.path(out_dir, "top80_functional_class_summary.csv")), "\n")










# =========================================================
# Chemoproteomics enrichment pipeline
# For raw EC50 table with columns: Name, pEC50
# Example Name format: PKM_0_49 / USP9X_10095_842
# =========================================================

# =========================================================
# 0. install / load packages
# =========================================================
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

cran_pkgs <- c("readr", "dplyr", "stringr", "ggplot2")
bioc_pkgs <- c("clusterProfiler", "org.Hs.eg.db", "ReactomePA", "enrichplot")

for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

library(readr)
library(dplyr)
library(stringr)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(enrichplot)

# =========================================================
# 1. file paths
# =========================================================
input_file <- "E:/python/curvecurator/8171/8171_EC50.csv"
out_dir    <- "E:/python/curvecurator/8171/enrichment_from_raw_name_pEC50"

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# =========================================================
# 2. user-adjustable parameters
# =========================================================
# 从 site-level 汇总到 gene-level 后，再筛 candidate hits
min_sites_cutoff <- 2
min_pEC50_cutoff <- 4.0

# 对于 13615 行的原始数据，gene-level 后通常可以用这些 topN
top_n_list <- c(100, 200, 400)

# enrichment cutoff
go_p_cut    <- 0.05
go_q_cut    <- 0.20
kegg_p_cut  <- 0.05
react_p_cut <- 0.05

show_n_plot <- 15

# =========================================================
# 3. read raw data
# =========================================================
df_raw <- read_csv(input_file, show_col_types = FALSE)

cat("====================================\n")
cat("Raw input loaded\n")
cat("====================================\n")
cat("Columns in raw file:\n")
print(colnames(df_raw))
cat("Raw rows:", nrow(df_raw), "\n\n")

required_cols <- c("Name", "pEC50")
missing_cols <- setdiff(required_cols, colnames(df_raw))

if (length(missing_cols) > 0) {
  stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
}

# 保存原始表拷贝
write.csv(df_raw, file.path(out_dir, "raw_input_copy.csv"), row.names = FALSE)

# =========================================================
# 4. clean raw data
# =========================================================
df_raw_clean <- df_raw %>%
  filter(!is.na(Name), Name != "") %>%
  filter(!is.na(pEC50)) %>%
  mutate(
    Name = str_trim(Name),
    pEC50 = as.numeric(pEC50)
  ) %>%
  filter(Name != "")

cat("Rows after raw cleaning:", nrow(df_raw_clean), "\n\n")

# =========================================================
# 5. extract gene_symbol from Name
# =========================================================
# 假设 Name 格式像：
# PKM_0_49
# USP9X_10095_842
# 取第一个下划线前面的部分作为 gene_symbol

df_site <- df_raw_clean %>%
  mutate(
    gene_symbol = str_extract(Name, "^[^_]+")
  ) %>%
  filter(!is.na(gene_symbol), gene_symbol != "")

cat("Rows after gene_symbol extraction:", nrow(df_site), "\n")
cat("Unique gene symbols:", n_distinct(df_site$gene_symbol), "\n\n")

# 保存带 gene_symbol 的 site-level 表
write.csv(
  df_site,
  file.path(out_dir, "site_level_with_gene_symbol.csv"),
  row.names = FALSE
)

# =========================================================
# 6. summarize to gene-level
# =========================================================
# 这里定义：
# best_pEC50 = 每个基因下最大的 pEC50
# n_sites    = 每个基因下 site 条目数

df_gene <- df_site %>%
  group_by(gene_symbol) %>%
  summarise(
    best_pEC50 = max(pEC50, na.rm = TRUE),
    mean_pEC50 = mean(pEC50, na.rm = TRUE),
    median_pEC50 = median(pEC50, na.rm = TRUE),
    n_sites = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(best_pEC50), desc(n_sites))

cat("====================================\n")
cat("Gene-level summary generated\n")
cat("====================================\n")
cat("Gene-level rows:", nrow(df_gene), "\n\n")

write.csv(
  df_gene,
  file.path(out_dir, "gene_level_pEC50_summary_generated.csv"),
  row.names = FALSE
)

# =========================================================
# =========================================================
# 7. define background and candidate hits
# =========================================================
background_genes <- as.character(df_gene$gene_symbol)

candidate_hits <- df_gene %>%
  filter(!is.na(best_pEC50), !is.na(n_sites)) %>%
  filter(n_sites >= min_sites_cutoff, best_pEC50 >= min_pEC50_cutoff) %>%
  arrange(desc(best_pEC50), desc(n_sites)) %>%
  as.data.frame()

cat("====================================\n")
cat("Candidate hit filtering summary\n")
cat("====================================\n")
cat("Background genes:", length(background_genes), "\n")
cat("Candidate hits after cutoff:", nrow(candidate_hits), "\n")
cat("Current cutoff: n_sites >=", min_sites_cutoff,
    ", best_pEC50 >=", min_pEC50_cutoff, "\n\n")

write.csv(
  data.frame(gene_symbol = background_genes),
  file.path(out_dir, "background_genes.csv"),
  row.names = FALSE
)

write.csv(
  candidate_hits,
  file.path(out_dir, "candidate_hits_after_cutoff.csv"),
  row.names = FALSE
)

# =========================================================
# 8. helper function: convert SYMBOL to ENTREZ
# =========================================================
convert_to_entrez <- function(genes) {
  genes <- unique(genes)
  genes <- genes[!is.na(genes)]
  genes <- genes[genes != ""]
  
  res <- tryCatch(
    {
      bitr(
        genes,
        fromType = "SYMBOL",
        toType   = "ENTREZID",
        OrgDb    = org.Hs.eg.db
      )
    },
    error = function(e) {
      data.frame()
    }
  )
  
  if (nrow(res) == 0) {
    return(data.frame())
  }
  
  res <- res %>% distinct(SYMBOL, .keep_all = TRUE)
  return(res)
}

bg_entrez <- convert_to_entrez(background_genes)

write.csv(
  bg_entrez,
  file.path(out_dir, "background_symbol_to_entrez.csv"),
  row.names = FALSE
)

cat("Background mapped to ENTREZ:", nrow(bg_entrez), "\n\n")

if (nrow(bg_entrez) == 0) {
  stop("No background genes could be mapped to ENTREZID. Please check gene symbol format.")
}

# =========================================================
# 9. helper functions for enrichment
# =========================================================
run_go <- function(hit_entrez, bg_entrez, p_cut = 0.05, q_cut = 0.20) {
  if (nrow(hit_entrez) == 0) return(NULL)
  
  enrichGO(
    gene          = hit_entrez$ENTREZID,
    universe      = bg_entrez$ENTREZID,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENTREZID",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = p_cut,
    qvalueCutoff  = q_cut,
    readable      = TRUE
  )
}

run_kegg <- function(hit_entrez, bg_entrez, p_cut = 0.05) {
  if (nrow(hit_entrez) == 0) return(NULL)
  
  kk <- enrichKEGG(
    gene          = hit_entrez$ENTREZID,
    universe      = bg_entrez$ENTREZID,
    organism      = "hsa",
    pvalueCutoff  = p_cut,
    pAdjustMethod = "BH"
  )
  
  kk_df <- as.data.frame(kk)
  if (nrow(kk_df) > 0) {
    kk <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  }
  return(kk)
}

run_reactome <- function(hit_entrez, bg_entrez, p_cut = 0.05) {
  if (nrow(hit_entrez) == 0) return(NULL)
  
  enrichPathway(
    gene          = hit_entrez$ENTREZID,
    universe      = bg_entrez$ENTREZID,
    organism      = "human",
    pvalueCutoff  = p_cut,
    pAdjustMethod = "BH",
    readable      = TRUE
  )
}

safe_df <- function(enrich_obj) {
  if (is.null(enrich_obj)) {
    return(data.frame())
  }
  df <- as.data.frame(enrich_obj)
  if (is.null(df) || nrow(df) == 0) {
    return(data.frame())
  }
  return(df)
}

save_dotplot <- function(enrich_obj, file_prefix, out_dir, show_n = 15) {
  df <- safe_df(enrich_obj)
  
  if (nrow(df) == 0) {
    cat(file_prefix, ": no result, skip plot.\n")
    return(NULL)
  }
  
  pdf(file.path(out_dir, paste0(file_prefix, "_dotplot.pdf")), width = 10, height = 8)
  print(dotplot(enrich_obj, showCategory = min(show_n, nrow(df))))
  dev.off()
}

write_summary_text <- function(file, title_text, hit_genes, hit_entrez, bg_entrez,
                               ego_df, ekegg_df, ereact_df) {
  sink(file)
  
  cat("====================================\n")
  cat(title_text, "\n")
  cat("====================================\n\n")
  
  cat("Input hit genes:", length(hit_genes), "\n")
  cat("Mapped hit genes:", nrow(hit_entrez), "\n")
  cat("Background mapped genes:", nrow(bg_entrez), "\n\n")
  
  cat("GO BP terms:", nrow(ego_df), "\n")
  cat("KEGG pathways:", nrow(ekegg_df), "\n")
  cat("Reactome pathways:", nrow(ereact_df), "\n\n")
  
  if (nrow(ego_df) > 0) {
    cat("Top GO terms:\n")
    print(head(ego_df[, c("Description", "Count", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "geneID")], 10))
    cat("\n")
  }
  
  if (nrow(ekegg_df) > 0) {
    cat("Top KEGG pathways:\n")
    print(head(ekegg_df[, c("Description", "Count", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "geneID")], 10))
    cat("\n")
  }
  
  if (nrow(ereact_df) > 0) {
    cat("Top Reactome pathways:\n")
    print(head(ereact_df[, c("Description", "Count", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "geneID")], 10))
    cat("\n")
  }
  
  sink()
}

# =========================================================
# 10. determine usable topN automatically
# =========================================================
available_hit_n <- nrow(candidate_hits)

if (available_hit_n == 0) {
  stop("No candidate hits after cutoff. Try relaxing min_sites_cutoff or min_pEC50_cutoff.")
}

valid_top_n <- top_n_list[top_n_list <= available_hit_n]

if (length(valid_top_n) == 0) {
  valid_top_n <- sort(unique(c(min(30, available_hit_n), min(50, available_hit_n))))
  valid_top_n <- valid_top_n[valid_top_n > 0]
}

cat("====================================\n")
cat("TopN settings\n")
cat("====================================\n")
cat("Requested topN:", paste(top_n_list, collapse = ", "), "\n")
cat("Usable topN:", paste(valid_top_n, collapse = ", "), "\n\n")

# =========================================================
# 11. compare different hit sizes
# =========================================================
summary_list <- list()

for (top_n in valid_top_n) {
  
  cat("\n============================\n")
  cat("Running enrichment for top", top_n, "\n")
  cat("============================\n")
  
  sub_dir <- file.path(out_dir, paste0("top", top_n))
  if (!dir.exists(sub_dir)) {
    dir.create(sub_dir, recursive = TRUE)
  }
  
  hit_df <- candidate_hits %>%
    slice(1:top_n)
  
  hit_genes <- hit_df$gene_symbol
  
  write.csv(
    hit_df,
    file.path(sub_dir, paste0("hit_top", top_n, "_table.csv")),
    row.names = FALSE
  )
  
  write.csv(
    data.frame(gene_symbol = hit_genes),
    file.path(sub_dir, paste0("hit_top", top_n, "_genes.csv")),
    row.names = FALSE
  )
  
  hit_entrez <- convert_to_entrez(hit_genes)
  
  write.csv(
    hit_entrez,
    file.path(sub_dir, paste0("hit_top", top_n, "_symbol_to_entrez.csv")),
    row.names = FALSE
  )
  
  cat("Input hit genes:", length(hit_genes), "\n")
  cat("Mapped hit genes:", nrow(hit_entrez), "\n")
  
  ego    <- run_go(hit_entrez, bg_entrez, p_cut = go_p_cut, q_cut = go_q_cut)
  ekegg  <- run_kegg(hit_entrez, bg_entrez, p_cut = kegg_p_cut)
  ereact <- run_reactome(hit_entrez, bg_entrez, p_cut = react_p_cut)
  
  ego_df    <- safe_df(ego)
  ekegg_df  <- safe_df(ekegg)
  ereact_df <- safe_df(ereact)
  
  write.csv(ego_df,    file.path(sub_dir, paste0("GO_BP_top", top_n, ".csv")), row.names = FALSE)
  write.csv(ekegg_df,  file.path(sub_dir, paste0("KEGG_top", top_n, ".csv")), row.names = FALSE)
  write.csv(ereact_df, file.path(sub_dir, paste0("Reactome_top", top_n, ".csv")), row.names = FALSE)
  
  save_dotplot(ego,    paste0("GO_BP_top", top_n), sub_dir, show_n = show_n_plot)
  save_dotplot(ekegg,  paste0("KEGG_top", top_n), sub_dir, show_n = show_n_plot)
  save_dotplot(ereact, paste0("Reactome_top", top_n), sub_dir, show_n = show_n_plot)
  
  top_go_desc       <- if (nrow(ego_df) > 0) ego_df$Description[1] else NA
  top_go_padj       <- if (nrow(ego_df) > 0) ego_df$p.adjust[1] else NA
  
  top_kegg_desc     <- if (nrow(ekegg_df) > 0) ekegg_df$Description[1] else NA
  top_kegg_padj     <- if (nrow(ekegg_df) > 0) ekegg_df$p.adjust[1] else NA
  
  top_react_desc    <- if (nrow(ereact_df) > 0) ereact_df$Description[1] else NA
  top_react_padj    <- if (nrow(ereact_df) > 0) ereact_df$p.adjust[1] else NA
  
  summary_list[[as.character(top_n)]] <- data.frame(
    top_n             = top_n,
    input_hit_genes   = length(hit_genes),
    mapped_hit_genes  = nrow(hit_entrez),
    GO_terms          = nrow(ego_df),
    KEGG_terms        = nrow(ekegg_df),
    Reactome_terms    = nrow(ereact_df),
    top_GO            = top_go_desc,
    top_GO_padj       = top_go_padj,
    top_KEGG          = top_kegg_desc,
    top_KEGG_padj     = top_kegg_padj,
    top_Reactome      = top_react_desc,
    top_Reactome_padj = top_react_padj,
    stringsAsFactors  = FALSE
  )
  
  write_summary_text(
    file = file.path(sub_dir, paste0("summary_top", top_n, ".txt")),
    title_text = paste("Enrichment summary for top", top_n),
    hit_genes = hit_genes,
    hit_entrez = hit_entrez,
    bg_entrez = bg_entrez,
    ego_df = ego_df,
    ekegg_df = ekegg_df,
    ereact_df = ereact_df
  )
}

# =========================================================
# 12. merge global summary
# =========================================================
summary_df <- bind_rows(summary_list)

write.csv(
  summary_df,
  file.path(out_dir, "enrichment_comparison_summary.csv"),
  row.names = FALSE
)

cat("\n====================================\n")
cat("Final comparison summary\n")
cat("====================================\n")
print(summary_df)

# =========================================================
# 13. save parameters
# =========================================================
param_df <- data.frame(
  parameter = c(
    "input_file",
    "out_dir",
    "min_sites_cutoff",
    "min_pEC50_cutoff",
    "top_n_list",
    "go_p_cut",
    "go_q_cut",
    "kegg_p_cut",
    "react_p_cut",
    "show_n_plot"
  ),
  value = c(
    input_file,
    out_dir,
    min_sites_cutoff,
    min_pEC50_cutoff,
    paste(top_n_list, collapse = ","),
    go_p_cut,
    go_q_cut,
    kegg_p_cut,
    react_p_cut,
    show_n_plot
  ),
  stringsAsFactors = FALSE
)

write.csv(
  param_df,
  file.path(out_dir, "pipeline_parameters.csv"),
  row.names = FALSE
)

cat("\n====================================\n")
cat("Pipeline finished.\n")
cat("All results saved to:\n")
cat(out_dir, "\n")
cat("====================================\n")



# =========================================================
# Chemoproteomics GSEA pipeline
# Input table columns: Name, pEC50
# Suitable for raw EC50 table like 8171_EC50.csv
# =========================================================

# =========================================================
# 0. install / load packages
# =========================================================
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

cran_pkgs <- c("readr", "dplyr", "stringr", "ggplot2")
bioc_pkgs <- c("clusterProfiler", "org.Hs.eg.db", "ReactomePA", "enrichplot", "DOSE")

for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

library(readr)
library(dplyr)
library(stringr)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(enrichplot)
library(DOSE)

# =========================================================
# 1. file paths
# =========================================================
input_file <- "E:/python/curvecurator/8171/8171_EC50.csv"
out_dir    <- "E:/python/curvecurator/8171/GSEA_from_raw_name_pEC50"

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# =========================================================
# 2. parameters
# =========================================================
# gene-level summary generation
min_valid_pEC50 <- 0

# GSEA parameters
minGSSize <- 10
maxGSSize <- 500
pvalueCutoff <- 0.2
nPermSimple <- 10000

# plotting
show_n_plot <- 20

# =========================================================
# 3. read raw data
# =========================================================
df_raw <- read_csv(input_file, show_col_types = FALSE)

cat("====================================\n")
cat("Raw input loaded\n")
cat("====================================\n")
cat("Columns in raw file:\n")
print(colnames(df_raw))
cat("Raw rows:", nrow(df_raw), "\n\n")

required_cols <- c("Name", "pEC50")
missing_cols <- setdiff(required_cols, colnames(df_raw))

if (length(missing_cols) > 0) {
  stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
}

write.csv(df_raw, file.path(out_dir, "raw_input_copy.csv"), row.names = FALSE)

# =========================================================
# 4. clean raw data
# =========================================================
df_raw_clean <- df_raw %>%
  filter(!is.na(Name), Name != "") %>%
  filter(!is.na(pEC50)) %>%
  mutate(
    Name = str_trim(Name),
    pEC50 = as.numeric(pEC50)
  ) %>%
  filter(Name != "") %>%
  filter(is.finite(pEC50)) %>%
  filter(pEC50 >= min_valid_pEC50) %>%
  as.data.frame()

cat("Rows after raw cleaning:", nrow(df_raw_clean), "\n\n")

# =========================================================
# 5. extract gene_symbol from Name
# =========================================================
# Example:
# PKM_0_49        -> PKM
# USP9X_10095_842 -> USP9X

df_site <- df_raw_clean %>%
  mutate(
    gene_symbol = str_extract(Name, "^[^_]+")
  ) %>%
  filter(!is.na(gene_symbol), gene_symbol != "") %>%
  as.data.frame()

cat("Rows after gene_symbol extraction:", nrow(df_site), "\n")
cat("Unique gene symbols:", length(unique(df_site$gene_symbol)), "\n\n")

write.csv(
  df_site,
  file.path(out_dir, "site_level_with_gene_symbol.csv"),
  row.names = FALSE
)

# =========================================================
# 6. gene-level summary
# =========================================================
df_gene <- df_site %>%
  group_by(gene_symbol) %>%
  summarise(
    best_pEC50   = max(pEC50, na.rm = TRUE),
    mean_pEC50   = mean(pEC50, na.rm = TRUE),
    median_pEC50 = median(pEC50, na.rm = TRUE),
    n_sites      = n(),
    weighted_score = max(pEC50, na.rm = TRUE) * log2(n() + 1),
    .groups      = "drop"
  ) %>%
  arrange(desc(best_pEC50), desc(n_sites)) %>%
  as.data.frame()

cat("====================================\n")
cat("Gene-level summary generated\n")
cat("====================================\n")
cat("Gene-level rows:", nrow(df_gene), "\n\n")

write.csv(
  df_gene,
  file.path(out_dir, "gene_level_summary_for_GSEA.csv"),
  row.names = FALSE
)

# =========================================================
# 7. SYMBOL -> ENTREZ mapping
# =========================================================
symbol_to_entrez <- tryCatch(
  {
    bitr(
      unique(as.character(df_gene$gene_symbol)),
      fromType = "SYMBOL",
      toType   = "ENTREZID",
      OrgDb    = org.Hs.eg.db
    )
  },
  error = function(e) {
    data.frame()
  }
)

if (nrow(symbol_to_entrez) == 0) {
  stop("No SYMBOL could be mapped to ENTREZID. Please check gene_symbol extraction.")
}

symbol_to_entrez <- symbol_to_entrez %>%
  distinct(SYMBOL, .keep_all = TRUE) %>%
  as.data.frame()

write.csv(
  symbol_to_entrez,
  file.path(out_dir, "symbol_to_entrez_mapping.csv"),
  row.names = FALSE
)

cat("Mapped SYMBOL -> ENTREZ:", nrow(symbol_to_entrez), "\n\n")

# =========================================================
# 8. merge mapping into gene-level table
# =========================================================
df_gene_mapped <- df_gene %>%
  left_join(symbol_to_entrez, by = c("gene_symbol" = "SYMBOL")) %>%
  filter(!is.na(ENTREZID), ENTREZID != "") %>%
  as.data.frame()

cat("Gene-level rows after ENTREZ mapping:", nrow(df_gene_mapped), "\n\n")

write.csv(
  df_gene_mapped,
  file.path(out_dir, "gene_level_summary_mapped.csv"),
  row.names = FALSE
)

# =========================================================
# 9. helper: build ranked geneList
# =========================================================
# GSEA requires a named numeric vector:
# values = ranking score
# names  = ENTREZID
# sorted decreasingly
#
# If duplicated ENTREZ exists after mapping, keep the highest score.

make_geneList <- function(df, score_col) {
  tmp <- df[, c("ENTREZID", score_col)]
  colnames(tmp) <- c("ENTREZID", "score")
  
  tmp <- tmp %>%
    filter(!is.na(score), is.finite(score)) %>%
    group_by(ENTREZID) %>%
    summarise(score = max(score, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(score)) %>%
    as.data.frame()
  
  geneList <- tmp$score
  names(geneList) <- tmp$ENTREZID
  
  geneList <- sort(geneList, decreasing = TRUE)
  return(geneList)
}

geneList_best     <- make_geneList(df_gene_mapped, "best_pEC50")
geneList_mean     <- make_geneList(df_gene_mapped, "mean_pEC50")
geneList_weighted <- make_geneList(df_gene_mapped, "weighted_score")

write.csv(
  data.frame(ENTREZID = names(geneList_best), score = as.numeric(geneList_best)),
  file.path(out_dir, "rank_best_pEC50.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(ENTREZID = names(geneList_mean), score = as.numeric(geneList_mean)),
  file.path(out_dir, "rank_mean_pEC50.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(ENTREZID = names(geneList_weighted), score = as.numeric(geneList_weighted)),
  file.path(out_dir, "rank_weighted_score.csv"),
  row.names = FALSE
)

# =========================================================
# 10. helper functions for GSEA
# =========================================================
safe_gsea_df <- function(gsea_obj) {
  if (is.null(gsea_obj)) {
    return(data.frame())
  }
  df <- tryCatch(as.data.frame(gsea_obj), error = function(e) data.frame())
  if (is.null(df) || nrow(df) == 0) {
    return(data.frame())
  }
  return(as.data.frame(df))
}

run_gsea_go <- function(geneList) {
  tryCatch(
    {
      gseGO(
        geneList      = geneList,
        OrgDb         = org.Hs.eg.db,
        keyType       = "ENTREZID",
        ont           = "BP",
        minGSSize     = minGSSize,
        maxGSSize     = maxGSSize,
        pvalueCutoff  = pvalueCutoff,
        pAdjustMethod = "BH",
        verbose       = FALSE,
        by            = "fgsea"
      )
    },
    error = function(e) NULL
  )
}

run_gsea_kegg <- function(geneList) {
  tryCatch(
    {
      gseKEGG(
        geneList      = geneList,
        organism      = "hsa",
        minGSSize     = minGSSize,
        maxGSSize     = maxGSSize,
        pvalueCutoff  = pvalueCutoff,
        pAdjustMethod = "BH",
        verbose       = FALSE,
        by            = "fgsea"
      )
    },
    error = function(e) NULL
  )
}

run_gsea_reactome <- function(geneList) {
  tryCatch(
    {
      gsePathway(
        geneList      = geneList,
        organism      = "human",
        minGSSize     = minGSSize,
        maxGSSize     = maxGSSize,
        pvalueCutoff  = pvalueCutoff,
        pAdjustMethod = "BH",
        verbose       = FALSE,
        by            = "fgsea"
      )
    },
    error = function(e) NULL
  )
}

save_gsea_dotplot <- function(gsea_obj, file_prefix, out_dir, show_n = 20) {
  df <- safe_gsea_df(gsea_obj)
  if (nrow(df) == 0) {
    cat(file_prefix, ": no result, skip dotplot.\n")
    return(NULL)
  }
  
  pdf(file.path(out_dir, paste0(file_prefix, "_dotplot.pdf")), width = 10, height = 8)
  print(dotplot(gsea_obj, showCategory = min(show_n, nrow(df)), split = ".sign"))
  dev.off()
}

save_gsea_ridgeplot <- function(gsea_obj, file_prefix, out_dir, show_n = 20) {
  df <- safe_gsea_df(gsea_obj)
  if (nrow(df) == 0) {
    cat(file_prefix, ": no result, skip ridgeplot.\n")
    return(NULL)
  }
  
  pdf(file.path(out_dir, paste0(file_prefix, "_ridgeplot.pdf")), width = 11, height = 8)
  print(ridgeplot(gsea_obj, showCategory = min(show_n, nrow(df))))
  dev.off()
}

save_top_running_plots <- function(gsea_obj, geneList, file_prefix, out_dir, top_n = 5) {
  df <- safe_gsea_df(gsea_obj)
  if (nrow(df) == 0) {
    cat(file_prefix, ": no result, skip running plots.\n")
    return(NULL)
  }
  
  df <- df[order(df$p.adjust, -abs(df$NES)), , drop = FALSE]
  top_ids <- utils::head(df$ID, top_n)
  
  for (pid in top_ids) {
    safe_name <- gsub("[^A-Za-z0-9_\\-]", "_", pid)
    pdf(file.path(out_dir, paste0(file_prefix, "_", safe_name, "_runningScore.pdf")), width = 10, height = 7)
    print(gseaplot2(gsea_obj, geneSetID = pid, title = pid))
    dev.off()
  }
}

write_gsea_summary <- function(file, title_text, gsea_df) {
  sink(file)
  
  cat("====================================\n")
  cat(title_text, "\n")
  cat("====================================\n\n")
  
  cat("Total enriched terms:", nrow(gsea_df), "\n\n")
  
  if (nrow(gsea_df) > 0) {
    keep_cols <- intersect(
      c("ID", "Description", "setSize", "enrichmentScore", "NES", "pvalue", "p.adjust", "qvalues", "core_enrichment"),
      colnames(gsea_df)
    )
    print(utils::head(gsea_df[, keep_cols, drop = FALSE], 20))
  }
  
  sink()
}

# =========================================================
# 11. run GSEA for multiple ranking strategies
# =========================================================
rank_list <- list(
  best_pEC50     = geneList_best,
  mean_pEC50     = geneList_mean,
  weighted_score = geneList_weighted
)

summary_list <- list()

for (rank_name in names(rank_list)) {
  
  cat("\n====================================\n")
  cat("Running GSEA for ranking:", rank_name, "\n")
  cat("====================================\n")
  
  geneList <- rank_list[[rank_name]]
  sub_dir <- file.path(out_dir, rank_name)
  
  if (!dir.exists(sub_dir)) {
    dir.create(sub_dir, recursive = TRUE)
  }
  
  cat("Ranked genes:", length(geneList), "\n")
  
  # -------------------------
  # GO BP
  # -------------------------
  gsea_go <- run_gsea_go(geneList)
  gsea_go_df <- safe_gsea_df(gsea_go)
  
  write.csv(
    gsea_go_df,
    file.path(sub_dir, paste0("GSEA_GO_BP_", rank_name, ".csv")),
    row.names = FALSE
  )
  
  save_gsea_dotplot(gsea_go, paste0("GSEA_GO_BP_", rank_name), sub_dir, show_n = show_n_plot)
  save_gsea_ridgeplot(gsea_go, paste0("GSEA_GO_BP_", rank_name), sub_dir, show_n = show_n_plot)
  save_top_running_plots(gsea_go, geneList, paste0("GSEA_GO_BP_", rank_name), sub_dir, top_n = 5)
  
  write_gsea_summary(
    file = file.path(sub_dir, paste0("summary_GSEA_GO_BP_", rank_name, ".txt")),
    title_text = paste("GSEA GO BP summary -", rank_name),
    gsea_df = gsea_go_df
  )
  
  # -------------------------
  # KEGG
  # -------------------------
  gsea_kegg <- run_gsea_kegg(geneList)
  gsea_kegg_df <- safe_gsea_df(gsea_kegg)
  
  write.csv(
    gsea_kegg_df,
    file.path(sub_dir, paste0("GSEA_KEGG_", rank_name, ".csv")),
    row.names = FALSE
  )
  
  save_gsea_dotplot(gsea_kegg, paste0("GSEA_KEGG_", rank_name), sub_dir, show_n = show_n_plot)
  save_gsea_ridgeplot(gsea_kegg, paste0("GSEA_KEGG_", rank_name), sub_dir, show_n = show_n_plot)
  save_top_running_plots(gsea_kegg, geneList, paste0("GSEA_KEGG_", rank_name), sub_dir, top_n = 5)
  
  write_gsea_summary(
    file = file.path(sub_dir, paste0("summary_GSEA_KEGG_", rank_name, ".txt")),
    title_text = paste("GSEA KEGG summary -", rank_name),
    gsea_df = gsea_kegg_df
  )
  
  # -------------------------
  # Reactome
  # -------------------------
  gsea_react <- run_gsea_reactome(geneList)
  gsea_react_df <- safe_gsea_df(gsea_react)
  
  write.csv(
    gsea_react_df,
    file.path(sub_dir, paste0("GSEA_Reactome_", rank_name, ".csv")),
    row.names = FALSE
  )
  
  save_gsea_dotplot(gsea_react, paste0("GSEA_Reactome_", rank_name), sub_dir, show_n = show_n_plot)
  save_gsea_ridgeplot(gsea_react, paste0("GSEA_Reactome_", rank_name), sub_dir, show_n = show_n_plot)
  save_top_running_plots(gsea_react, geneList, paste0("GSEA_Reactome_", rank_name), sub_dir, top_n = 5)
  
  write_gsea_summary(
    file = file.path(sub_dir, paste0("summary_GSEA_Reactome_", rank_name, ".txt")),
    title_text = paste("GSEA Reactome summary -", rank_name),
    gsea_df = gsea_react_df
  )
  
  # -------------------------
  # Global summary row
  # -------------------------
  top_go_desc    <- if (nrow(gsea_go_df) > 0) gsea_go_df$Description[1] else NA
  top_go_nes     <- if (nrow(gsea_go_df) > 0) gsea_go_df$NES[1] else NA
  top_go_padj    <- if (nrow(gsea_go_df) > 0) gsea_go_df$p.adjust[1] else NA
  
  top_kegg_desc  <- if (nrow(gsea_kegg_df) > 0) gsea_kegg_df$Description[1] else NA
  top_kegg_nes   <- if (nrow(gsea_kegg_df) > 0) gsea_kegg_df$NES[1] else NA
  top_kegg_padj  <- if (nrow(gsea_kegg_df) > 0) gsea_kegg_df$p.adjust[1] else NA
  
  top_react_desc <- if (nrow(gsea_react_df) > 0) gsea_react_df$Description[1] else NA
  top_react_nes  <- if (nrow(gsea_react_df) > 0) gsea_react_df$NES[1] else NA
  top_react_padj <- if (nrow(gsea_react_df) > 0) gsea_react_df$p.adjust[1] else NA
  
  summary_list[[rank_name]] <- data.frame(
    ranking_method       = rank_name,
    ranked_genes         = length(geneList),
    GO_terms             = nrow(gsea_go_df),
    top_GO               = top_go_desc,
    top_GO_NES           = top_go_nes,
    top_GO_padj          = top_go_padj,
    KEGG_terms           = nrow(gsea_kegg_df),
    top_KEGG             = top_kegg_desc,
    top_KEGG_NES         = top_kegg_nes,
    top_KEGG_padj        = top_kegg_padj,
    Reactome_terms       = nrow(gsea_react_df),
    top_Reactome         = top_react_desc,
    top_Reactome_NES     = top_react_nes,
    top_Reactome_padj    = top_react_padj,
    stringsAsFactors     = FALSE
  )
}

# =========================================================
# 12. merge final summary
# =========================================================
summary_df <- bind_rows(summary_list) %>% as.data.frame()

write.csv(
  summary_df,
  file.path(out_dir, "GSEA_comparison_summary.csv"),
  row.names = FALSE
)

cat("\n====================================\n")
cat("Final GSEA comparison summary\n")
cat("====================================\n")
print(summary_df)

# =========================================================
# 13. save parameters
# =========================================================
param_df <- data.frame(
  parameter = c(
    "input_file",
    "out_dir",
    "min_valid_pEC50",
    "minGSSize",
    "maxGSSize",
    "pvalueCutoff",
    "nPermSimple",
    "show_n_plot"
  ),
  value = c(
    input_file,
    out_dir,
    min_valid_pEC50,
    minGSSize,
    maxGSSize,
    pvalueCutoff,
    nPermSimple,
    show_n_plot
  ),
  stringsAsFactors = FALSE
)

write.csv(
  param_df,
  file.path(out_dir, "GSEA_parameters.csv"),
  row.names = FALSE
)

cat("\n====================================\n")
cat("GSEA pipeline finished.\n")
cat("All results saved to:\n")
cat(out_dir, "\n")
cat("====================================\n")


library(readr)
library(dplyr)
library(ggplot2)
library(forcats)

# ==============================
# 1. 读取 GSEA 结果
# ==============================
df <- read_csv("E:/python/curvecurator/8171/GSEA_from_raw_name_pEC50/weighted_score/GSEA_GO_BP_weighted_score.csv")

# 看一下列名，确认有 Description, NES, p.adjust, setSize
colnames(df)

# ==============================
# 2. 取 top 15 term
#    这里按 p.adjust 最小排序，更适合论文图
# ==============================
df_top <- df %>%
  filter(!is.na(Description), !is.na(NES), !is.na(p.adjust)) %>%
  arrange(p.adjust, desc(abs(NES))) %>%
  head(15) %>%
  mutate(
    neglog10_padj = -log10(p.adjust),
    Description = fct_reorder(Description, NES)
  )

# ==============================
# 3. 横向 barplot
#    每个 bar 不同颜色
# ==============================
ggplot(df_top, aes(x = Description, y = NES, fill = Description)) +
  geom_col(width = 0.75) +
  coord_flip() +
  labs(
    title = "Top GO Biological Process Enrichment (GSEA)",
    x = "GO Term",
    y = "Normalized Enrichment Score (NES)"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10)
  )

# ==============================
# 4. Dotplot
#    x = NES
#    y = GO term
#    点大小 = setSize
#    点颜色 = -log10(adj.p)
# ==============================
ggplot(df_top, aes(x = NES, y = Description)) +
  geom_point(aes(size = setSize, color = neglog10_padj)) +
  labs(
    title = "Top GO Biological Process Enrichment (Dotplot)",
    x = "Normalized Enrichment Score (NES)",
    y = "GO Term",
    color = "-log10(adj.p)",
    size = "Gene Set Size"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10)
  )


library(readr)
library(dplyr)
library(stringr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

# ==============================
# 1. 读取原始 EC50 数据
# ==============================
df_raw <- read_csv("E:/python/curvecurator/8171/8171_EC50.csv", show_col_types = FALSE)

# ==============================
# 2. 提取 gene_symbol
# ==============================
df_site <- df_raw %>%
  filter(!is.na(Name), !is.na(pEC50)) %>%
  mutate(
    Name = str_trim(Name),
    gene_symbol = str_extract(Name, "^[^_]+"),
    pEC50 = as.numeric(pEC50)
  ) %>%
  filter(!is.na(gene_symbol), gene_symbol != "") %>%
  as.data.frame()

# ==============================
# 3. gene-level summary
# ==============================
df_gene <- df_site %>%
  group_by(gene_symbol) %>%
  summarise(
    best_pEC50 = max(pEC50, na.rm = TRUE),
    n_sites = n(),
    weighted_score = max(pEC50, na.rm = TRUE) * log2(n() + 1),
    .groups = "drop"
  ) %>%
  as.data.frame()

# ==============================
# 4. SYMBOL -> ENTREZ
# ==============================
map_df <- bitr(
  unique(df_gene$gene_symbol),
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
) %>%
  distinct(SYMBOL, .keep_all = TRUE)

df_gene_mapped <- df_gene %>%
  left_join(map_df, by = c("gene_symbol" = "SYMBOL")) %>%
  filter(!is.na(ENTREZID)) %>%
  as.data.frame()

# ==============================
# 5. 构建 geneList
# ==============================
geneList <- df_gene_mapped %>%
  select(ENTREZID, weighted_score) %>%
  group_by(ENTREZID) %>%
  summarise(weighted_score = max(weighted_score, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(weighted_score))

geneList_vec <- geneList$weighted_score
names(geneList_vec) <- geneList$ENTREZID
geneList_vec <- sort(geneList_vec, decreasing = TRUE)

# ==============================
# 6. 重新跑 gseGO
# ==============================
gsea_go <- gseGO(
  geneList      = geneList_vec,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  minGSSize     = 10,
  maxGSSize     = 500,
  pvalueCutoff  = 0.2,
  pAdjustMethod = "BH",
  verbose       = FALSE,
  by            = "fgsea"
)

# 看前几个 term
head(as.data.frame(gsea_go)[, c("ID", "Description", "NES", "p.adjust")])

# ==============================
# 7. 画 running score plot
#    这里用第1个term，也可以改成第2、第3个
# ==============================
gseaplot2(
  gsea_go,
  geneSetID = 1,
  title = as.data.frame(gsea_go)$Description[1]
)



library(readr)
library(dplyr)
library(ggplot2)
library(stringr)

# ==============================
# 1. 读入文件
# ==============================
df <- read_csv("E:/python/curvecurator/8171/GSEA_from_raw_name_pEC50/weighted_score/top5_pathway_driver_genes_with_scores.csv")

# 如果你现在文件就在工作目录，也可以用：
# df <- read_csv("top5_pathway_driver_genes_with_scores.csv")

# ==============================
# 2. 基础清理
# ==============================
df2 <- df %>%
  filter(!is.na(Pathway_Description),
         !is.na(SYMBOL),
         !is.na(weighted_score)) %>%
  mutate(
    Pathway_Description = str_wrap(Pathway_Description, width = 35)
  )

# 按 pathway 显著性和 score 排序
pathway_order <- df2 %>%
  distinct(Pathway_Description, p.adjust) %>%
  arrange(p.adjust) %>%
  pull(Pathway_Description)

gene_order <- df2 %>%
  arrange(weighted_score) %>%
  pull(SYMBOL) %>%
  unique()

df2$Pathway_Description <- factor(df2$Pathway_Description, levels = pathway_order)
df2$SYMBOL <- factor(df2$SYMBOL, levels = gene_order)

# ==============================
# 3. 画 heatmap
# ==============================
ggplot(df2, aes(x = Pathway_Description, y = SYMBOL, fill = weighted_score)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = round(best_pEC50, 2)), size = 3) +
  labs(
    title = "Driver Genes of Top GSEA Pathways",
    x = "Pathway",
    y = "Driver Gene",
    fill = "Weighted score"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 10, angle = 35, hjust = 1),
    axis.text.y = element_text(size = 9),
    panel.grid = element_blank()
  )