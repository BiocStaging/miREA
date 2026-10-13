## Regenerates the toy_* example datasets shipped in data/.
## Not run during R CMD build/check; re-run manually if the toy data needs
## to be rebuilt. Requires data.table, dplyr, DESeq2, and access to the raw
## TCGA BLCA data under raw_data/TCGA/ (see project root, not part of this
## package).

suppressMessages({
  library(data.table)
  library(dplyr)
  library(DESeq2)
})

raw_tcga_dir <- "/scratch/project_2011179/raw_data/TCGA"
full_data_path <- "/scratch/project_2011179/code/miREA_review_ploscb/miREA/data/input_data/hallmark/BLCA_hallmark_input_data.RData"

patient_id <- function(x) substr(x, 1, 12)

## ---- 1. Sample selection: 80 tumor (TP) + 19 normal (NT) ------------------
## NT is capped at 19: that's every TCGA-BLCA patient with both a Solid
## Tissue Normal RNA-seq AND miRNA-seq sample. 80 TP was chosen empirically:
## fewer TP samples leave too few patients for the miRNA-gene correlation
## step to produce any pair surviving BH correction across ~1e5 candidate
## pairs (tested 40 -> 0 significant pairs, 50 -> 20, 80 -> 190).
gene_samp_clean <- read.csv(file.path(raw_tcga_dir, "gene/BLCA_geneSample.csv"))
mir_samp_clean  <- read.csv(file.path(raw_tcga_dir, "miRNA/BLCA_mirSample.csv"))

gene_tp_patients <- patient_id(gene_samp_clean$col[gene_samp_clean$label == "Primary Tumor"])
gene_nt_patients <- patient_id(gene_samp_clean$col[gene_samp_clean$label == "Solid Tissue Normal" & !grepl("^SRR", gene_samp_clean$col)])
mir_tp_patients  <- patient_id(mir_samp_clean$col[mir_samp_clean$label == "Primary Tumor"])
mir_nt_patients  <- patient_id(mir_samp_clean$col[mir_samp_clean$label == "Solid Tissue Normal"])

tp_candidates <- intersect(gene_tp_patients, mir_tp_patients)
nt_candidates <- intersect(gene_nt_patients, mir_nt_patients)

## raw count files: integer gene/miRNA counts with real IDs, both tumor and
## normal; `raw_data/TCGA/gene/BLCA_geneData.csv` is a dead-end export that
## lost its gene-ID column, so the 202503/ snapshot is used instead.
gene_hgnc_path <- file.path(raw_tcga_dir, "202503/gene/BLCA_HGNC.csv")
mir_full_path  <- file.path(raw_tcga_dir, "202503/miRNA/BLCA_full.csv")
gene_samp_raw  <- read.csv(file.path(raw_tcga_dir, "202503/gene/BLCA_samples.csv"))
mir_samp_raw   <- read.csv(file.path(raw_tcga_dir, "202503/miRNA/BLCA_samples.csv"))

gene_raw_tp <- gene_samp_raw[gene_samp_raw$sample_type == "Primary Tumor", ]
gene_raw_nt <- gene_samp_raw[gene_samp_raw$sample_type == "Solid Tissue Normal" & !grepl("^NTa", gene_samp_raw$col), ]
mir_raw_tp  <- mir_samp_raw[mir_samp_raw$sample_type == "Primary Tumor", ]
mir_raw_nt  <- mir_samp_raw[mir_samp_raw$sample_type == "Solid Tissue Normal", ]

## keep only patients with exactly one column on both modalities, so no
## technical-replicate averaging is needed
unique_patients_in <- function(df, id_col) {
  tab <- table(patient_id(df[[id_col]]))
  names(tab[tab == 1])
}
tp_pool <- Reduce(intersect, list(tp_candidates,
                                   unique_patients_in(gene_raw_tp, "ID"),
                                   unique_patients_in(mir_raw_tp, "cases.submitter_id")))
nt_pool <- Reduce(intersect, list(nt_candidates,
                                   unique_patients_in(gene_raw_nt, "ID"),
                                   unique_patients_in(mir_raw_nt, "cases.submitter_id")))

N_TP <- 80
set.seed(20261012)
pick_nt <- sort(nt_pool)
tp_pool <- setdiff(tp_pool, pick_nt)  # a patient can appear in both pools (paired tumor+normal)
pick_tp <- sort(sample(tp_pool, N_TP))
pick_patients <- c(pick_tp, pick_nt)
pick_type <- c(rep("TP", length(pick_tp)), rep("NT", length(pick_nt)))

gene_col_tp <- setNames(gene_raw_tp$col, patient_id(gene_raw_tp$ID))
gene_col_nt <- setNames(gene_raw_nt$col, patient_id(gene_raw_nt$ID))
mir_col_tp  <- setNames(mir_raw_tp$sample.submitter_id, patient_id(mir_raw_tp$cases.submitter_id))
mir_col_nt  <- setNames(mir_raw_nt$sample.submitter_id, patient_id(mir_raw_nt$cases.submitter_id))
gene_cols <- c(gene_col_tp[pick_tp], gene_col_nt[pick_nt])
mir_cols  <- substr(c(mir_col_tp[pick_tp], mir_col_nt[pick_nt]), 1, 15)  # drop the plate-letter suffix BLCA_full.csv doesn't use

## ---- 2. Extract raw counts, run DESeq2 ------------------------------------
gene_id_col <- colnames(fread(gene_hgnc_path, nrows = 0))[1]
gene_dt <- fread(gene_hgnc_path, select = c(gene_id_col, unname(gene_cols)))
setnames(gene_dt, gene_id_col, "gene")
setnames(gene_dt, unname(gene_cols), names(gene_cols))

mir_id_col <- colnames(fread(mir_full_path, nrows = 0))[1]
mir_dt <- fread(mir_full_path, select = c(mir_id_col, unname(mir_cols)))
setnames(mir_dt, mir_id_col, "miRNA")
setnames(mir_dt, unname(mir_cols), names(mir_cols))

run_deseq <- function(dt, id_col) {
  mat <- as.matrix(dt[, -1])
  rownames(mat) <- dt[[id_col]]
  mat <- round(mat)  # gene counts are RSEM expected counts, not integers
  storage.mode(mat) <- "integer"
  colData <- data.frame(type = factor(pick_type, levels = c("NT", "TP")))
  rownames(colData) <- colnames(mat)
  dds <- DESeqDataSetFromMatrix(countData = mat, colData = colData, design = ~ type)
  dds <- DESeq(dds, quiet = TRUE)
  res <- as.data.frame(results(dds))
  res[!is.na(res$padj), ]
}
gene_res <- run_deseq(gene_dt, "gene")
mir_res  <- run_deseq(mir_dt, "miRNA")

toy_gene_DEdata <- data.frame(gene = rownames(gene_res), log2FoldChange = gene_res$log2FoldChange,
                               stat = gene_res$stat, padj = gene_res$padj, stringsAsFactors = FALSE)
toy_mir_DEdata <- data.frame(miRNA = rownames(mir_res), log2FoldChange = mir_res$log2FoldChange,
                              stat = mir_res$stat, padj = mir_res$padj, stringsAsFactors = FALSE)

## ---- 3. Correlation matrices: TP samples only -----------------------------
## Raw counts, no extra normalization — matches the original manuscript
## pipeline (code/TCGAdata_generate.Rmd), which computes Spearman
## correlation directly on BLCA_{gene,mir}_TP.csv without a normalization step.
gene_universe <- toy_gene_DEdata$gene
mir_universe  <- toy_mir_DEdata$miRNA

toy_gene_mat <- as.data.frame(gene_dt[gene_dt$gene %in% gene_universe, c("gene", pick_tp), with = FALSE])
rownames(toy_gene_mat) <- toy_gene_mat$gene
toy_gene_mat$gene <- NULL

toy_mir_mat <- as.data.frame(mir_dt[mir_dt$miRNA %in% mir_universe, c("miRNA", pick_tp), with = FALSE])
rownames(toy_mir_mat) <- toy_mir_mat$miRNA
toy_mir_mat$miRNA <- NULL

## ---- 4. background_MGI / background_GGI: subset of the full manuscript ----
## background to the gene/miRNA universe above
load(full_data_path)  # provides `input_data`
full_bg_mgi <- input_data$background_MGI
full_bg_ggi <- input_data$background_GGI

toy_background_MGI <- unique(full_bg_mgi[full_bg_mgi$miRNA %in% mir_universe & full_bg_mgi$gene %in% gene_universe, ])
toy_background_GGI <- unique(full_bg_ggi[full_bg_ggi$from %in% gene_universe & full_bg_ggi$to %in% gene_universe, ])
colnames(toy_background_GGI) <- c("source_genesymbol", "target_genesymbol", "interact")

## pathway: BLCA-specific gene sets (not the generic MSigDB hallmark used
## above for background_MGI/GGI) -- much smaller per-pathway gene counts,
## and more relevant to a bladder-cancer toy example
tp_blca_path <- "/scratch/project_2011179/code/miREA_review_ploscb/miREA/data/raw_data/pathway/TP/cancer_specific/TP_BLCA_gene.csv"
full_pathway_tp_blca <- read.csv(tp_blca_path)
toy_pathway <- unique(full_pathway_tp_blca[full_pathway_tp_blca$gene %in% gene_universe, ])

## ---- 5. Save ---------------------------------------------------------------
usethis::use_data(toy_gene_DEdata, overwrite = TRUE, compress = "xz")
usethis::use_data(toy_mir_DEdata, overwrite = TRUE, compress = "xz")
usethis::use_data(toy_gene_mat, overwrite = TRUE, compress = "xz")
usethis::use_data(toy_mir_mat, overwrite = TRUE, compress = "xz")
usethis::use_data(toy_background_MGI, overwrite = TRUE, compress = "xz")
usethis::use_data(toy_background_GGI, overwrite = TRUE, compress = "xz")
usethis::use_data(toy_pathway, overwrite = TRUE, compress = "xz")
