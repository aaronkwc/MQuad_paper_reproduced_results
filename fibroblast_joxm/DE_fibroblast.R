library(SingleCellExperiment)
library(SummarizedExperiment)
library(Seurat)
library(edgeR)
library(dplyr)
library(scater)
library(IHW)
library(limma)
library(org.Hs.eg.db)
library(ggrepel)
xx <- as.list(org.Hs.egENSEMBL2EG)
load(url("http://bioinf.wehi.edu.au/software/MSigDB/human_H_v5p2.rdata"))


de_res <- readRDS("filt_lenient.cell_coverage_sites.de_results_unstimulated_cells.rds")
#get count matrix to redo the de analysis
count_matrix <- de_res$dge_list[['joxm']]$counts
new_assignment <- read.csv('clone_id.csv')
de_res$sce_list_unst[["joxm"]]$assigned <- new_assignment$combined

gene_use_DE <- (!rowData(de_res$sce_list_unst[["joxm"]])$is_feature_control) # remove ERCC
gene_use_DE <- gene_use_DE & (rowMeans(counts(de_res$sce_list_unst[["joxm"]]) > 1)) # all 0s
  # 1 count in 10% of cells
  # rowMeans of counts > 0.5

dge_list_unst <- edgeR::DGEList(round(counts(de_res$sce_list_unst[["joxm"]][gene_use_DE,])))
dge_list_unst <- edgeR::calcNormFactors(dge_list_unst, method = "TMM")
de_res$sce_list_unst[["joxm"]]$cdr <- (colSums(counts(de_res$sce_list_unst[["joxm"]]) > 0) / nrow(de_res$sce_list_unst[["joxm"]]))
de_res$sce_list_unst[["joxm"]]$plate <- as.factor(de_res$sce_list_unst[["joxm"]]$plate)
design_list_unst <- model.matrix(~cdr + plate + assigned, data = colData(de_res$sce_list_unst[["joxm"]]))
dge_list_unst <- estimateDisp(dge_list_unst, design_list_unst)
fit_list_unst <- glmQLFit(dge_list_unst, design_list_unst)
num_clones <- length(unique(de_res$sce_list_unst[["joxm"]]$assigned))
qlf_list_unst <- glmQLFTest(fit_list_unst,
                                 coef = (ncol(design_list_unst) - num_clones + 2):ncol(design_list_unst))
sum(p.adjust(qlf_list_unst$table$PValue, method = "BH") <= 0.05, na.rm = TRUE)
print(summary(decideTestsDGE(qlf_list_unst)))


num_clones <- length(unique(de_res$sce_list_unst[["joxm"]]$assigned))
qlf_1st_coef_list_unst <- glmQLFTest(fit_list_unst,
                                          coef = (ncol(design_list_unst) - num_clones + 2))
print(summary(decideTestsDGE(qlf_1st_coef_list_unst)))


### Testing first clone coeeficient alone

camera_msigdb_H_1st_coef_list_unst <- list()
fry_msigdb_H_1st_coef_list_unst <- list()

  if (num_clones > 1.5) {
    qlf_1st_coef_list_unst$table$ensembl_gene_id <- strsplit2(rownames(qlf_1st_coef_list_unst$table), split = "_")[,1]
    qlf_1st_coef_list_unst$table$hgnc_symbol <- strsplit2(rownames(qlf_1st_coef_list_unst$table), split = "_")[,2]
    qlf_1st_coef_list_unst$table$entrezid <- NA
    for (j in seq_len(nrow(qlf_1st_coef_list_unst$table))) {
      if (qlf_1st_coef_list_unst$table$ensembl_gene_id[j] %in% names(xx))
        qlf_1st_coef_list_unst$table$entrezid[j] <- xx[[qlf_1st_coef_list_unst$table$ensembl_gene_id[j]]][1]
    }
    idx <- ids2indices(Hs.H, id=qlf_1st_coef_list_unst$table$entrezid)
    length(idx)
    camera_msigdb_H_1st_coef_list_unst <- camera(dge_list_unst, idx, qlf_1st_coef_list_unst$design, 
                                                      contrast = (ncol(design_list_unst) - num_clones + 2), inter.gene.cor = 0.005, use.ranks=TRUE)
    cat("                camera significant MSigDB hallmark genesets (FDR < 10%): ", 
        sum(camera_msigdb_H_1st_coef_list_unst$FDR < 0.1), "\n")
    fry_msigdb_H_1st_coef_list_unst <- fry(dge_list_unst, idx, qlf_1st_coef_list_unst$design, 
                                                contrast = (ncol(design_list_unst) - num_clones + 2))    
    cat("                fry significant MSigDB hallmark genesets (FDR < 10%): ", 
        sum(fry_msigdb_H_1st_coef_list_unst$FDR < 0.1), "\n")
  }


## visualize DE and gene set  results
  
      i="clone1.MT1_clone1.MT0"
      cam_H_pw <- camera_msigdb_H_1st_coef_list_unst
      cam_H_pw[["geneset"]] <- rownames(camera_msigdb_H_1st_coef_list_unst)
      cam_H_pw <- cam_H_pw %>% 
        dplyr::mutate(sig = FDR < 0.1) 
      print(head(cam_H_pw))
      
      
      cam_H_pw[["lab"]] <- ""
      cam_H_pw[["lab"]][1:4] <-
        cam_H_pw[["geneset"]][1:4]
      cam_H_pw[["Direction"]][cam_H_pw[["Direction"]] == "Up"] <- 
        paste("Up in", strsplit(i, "_")[[1]][1], "vs", strsplit(i, "_")[[1]][2])
      cam_H_pw[["Direction"]][cam_H_pw[["Direction"]] == "Down"] <- 
        paste("Down in", strsplit(i, "_")[[1]][1], "vs", strsplit(i, "_")[[1]][2])
      
      
      p_hallmark <- cam_H_pw %>% 
        ggplot(aes(x = Direction, y = -log10(PValue), colour = sig, 
                   label = lab)) +
        ggbeeswarm::geom_quasirandom(aes(size = NGenes)) +
        geom_label_repel(show.legend = FALSE,
                         nudge_y = 0.3, nudge_x = 0.3, fill = "gray95") +
        scale_colour_manual(values = c("gray50", "firebrick"), 
                            label = c("N.S.", "FDR < 5%"), name = "") +
        guides(alpha = FALSE,
               fill = guide_legend(override.aes = list(size = 5))) +
        xlab("Gene set enrichment direction") +
        theme_classic(20) + theme(legend.position = "right")
      print(p_hallmark)
      
    
      de_tab <- qlf_pairwise_list[[i]]$table
      de_tab[["gene"]] <- rownames(de_tab)
      de_tab <- de_tab %>% 
        dplyr::mutate(FDR = adj_pvalues(ihw(PValue ~ logCPM, alpha = 0.1)), 
                      sig = FDR < 0.1,
                      signed_F = sign(logFC) * F) 
      de_tab[["lab"]] <- ""
      int_genes_entrezid <- c(Hs.H$HALLMARK_G2M_CHECKPOINT,Hs.H$E2F_TARGETS, Hs.H$HALLMARK_REACTIVE_OXIGEN_SPECIES_PATHWAY, Hs.H$HALLMARK_MYC_TARGETS_V1)
      mm <- match(int_genes_entrezid, de_tab$entrezid)
      mm <- mm[!is.na(mm)]
      int_genes_hgnc <- de_tab$hgnc_symbol[mm]
      int_genes_hgnc <- c(int_genes_hgnc, "MYBL1")
      genes_to_label <- (de_tab[["hgnc_symbol"]] %in% int_genes_hgnc
                           & de_tab[["FDR"]] < 0.1)
      de_tab[["lab"]][genes_to_label] <-
        de_tab[["hgnc_symbol"]][genes_to_label]
      
      p_vulcan <- ggplot(de_tab, aes(x = logFC, y = -log10(PValue), fill = sig,
                                     label = lab)) +
        geom_point(aes(size = sig), pch = 21, colour = "gray40") +
        geom_label_repel(show.legend = FALSE, 
                         arrow = arrow(type = "closed", length = unit(0.25, "cm")), 
                         nudge_x = 0.2, nudge_y = 0.3, fill = "gray95") +
        geom_segment(aes(x = -1, y = 0, xend = -4, yend = 0), 
                     colour = "black", size = 1, arrow = arrow(length = unit(0.5, "cm"))) +
        annotate("text", x = -4, y = -0.5, size = 6,
                 label = paste("higher in", strsplit(i, "_")[[1]][2])) +
        geom_segment(aes(x = 1, y = 0, xend = 4, yend = 0), 
                     colour = "black", size = 1, arrow = arrow(length = unit(0.5, "cm"))) +
        annotate("text", x = 4, y = -0.5, size = 6,
                 label = paste("higher in", strsplit(i, "_")[[1]][1])) +
        scale_fill_manual(values = c("gray60", "firebrick"), 
                          label = c("N.S.", "FDR < 10%"), name = "") +
        scale_size_manual(values = c(1, 3), guide = FALSE) +
        guides(alpha = FALSE,
               fill = guide_legend(override.aes = list(size = 5))) +
        theme_classic(20) + theme(legend.position = "right")
      print(p_vulcan)
      
#write to csv
write.csv(de_tab, file='de_tab.csv')
write.csv(cam_H_pw, file='cam_H_pw.csv')
