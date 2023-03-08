####################################
############ Libraries #############
####################################

# Data import
require(MultiAssayExperiment)
# Data Wrangling
library(tidyverse)
# RNAseq related tools
library(DESeq2)
library(limma)
library(igis) 
library(rnaseqTools)
# Visualization tools
library(EnhancedVolcano)
library(pheatmap) 
# Clustering and geneset analysis
library(Seurat)
library(fgsea)
library(GSVA)
# Database
library(biomaRt)
mart.hs = useDataset("hsapiens_gene_ensembl", useMart("ENSEMBL_MART_ENSEMBL",host="www.ensembl.org"))
library(msigdbr)
msig.hs = msigdbr(species = "Homo sapiens")


###############################################################
################### Part 1: Data Pre-process ##################
###############################################################

############################################
######### QC and create gSE ################
############################################

# Helper function: QC and gSummarizedExperiment (gSE)
gseQC <- function(gse) {
  qc=rnaseqTools::gQCproc(gse) 
  qc=gQCnorm(qc, method = "tpm", qcfilterApply = FALSE) 
  qc=qc[!qc@elementMetadata$Chromosome%in%c("X","Y","MT", "KI270711.1"),]
  qc=qc[!is.na(qc@elementMetadata$Chromosome),]
  qc@assays@data@listData$log2counts1 = log2(qc@assays@data@listData$counts+1)
  qc@assays@data@listData$log2tpm1 = log2(qc@assays@data@listData$tpm+1) 
  qc = qc[!grepl("LINC", qc@elementMetadata$HGNC),]  
  # Filter
  qc.all=subset(qc, LowExpressionFlag %in% c(0,1), LowDepthFlag == 0 & TechnicalFailureFlag == 0) 
  qc.exp=subset(qc, LowExpressionFlag == 0, LowDepthFlag == 0 & TechnicalFailureFlag == 0) 
  # Output
  results = list("all" = qc.all, "exp" = qc.exp)
  return(results) 
}

# Create gSE
gseCLL14RNAqc = gseQC(gseCLL14RNA)
gseCLL14RNAqc.all = gseCLL14RNAqc$all
gseCLL14RNAqc.exp = gseCLL14RNAqc$exp


#########################################
############# gSE to CDS ################
#########################################

# Helper 1: CDS helper function
getCDSfromGSE <- function(gse,removeLowLibSize=FALSE,dropLowCountGenes=TRUE,genesToKeep=c()){
  # Create CDS
  cds = DESeqDataSet(gse, design=~1) 
  
  # Handling genenames
  igis = Igis("human")
  genes = genesById(as.character(rownames(cds)), igis=igis)
  fannot = as.data.frame(genes@elementMetadata)
  colnames(fannot) = c("HGNC_Symbol","ID", "GeneName")
  rownames(fannot) = fannot$ID
  fannot = fannot[rownames(cds),]
  rowData(cds) = fannot
  
  if(!dropLowCountGenes){
    cds <- estimateSizeFactors( cds ) 
    return(cds)}
  
  # Handling exceptions
  tmp = rowSums(counts(cds)<3)
  ind <- which(tmp > ncol(cds)*0.95)
  keep = which(rowData(cds)$HGNC_Symbol %in% genesToKeep)
  ind = ind[!(ind %in% keep)]
  
  # Handling exceptions
  if(! isEmpty(ind)) { 
    cds2 <- cds[- ind,]
    cds2 <- estimateSizeFactors( cds2 )
  } else {
    cds2 <- estimateSizeFactors( cds )
  }
  
  # Handling exceptions
  if(removeLowLibSize){
    libSizes <- sizeFactors( cds2 ) * sum(counts(cds2)/ncol(cds2))
    i = libSizes<10*10^6
    cds2 = cds2[,!i]
  }
  
  return(cds2)
}

# Helper 2: Gene name annotation
cdsAnnotGse <- function(gse, mae) {
  cds = getCDSfromGSE(gse,dropLowCountGenes=TRUE, removeLowLibSize=TRUE) 
  
  # RNAseq-clinical data association 
  mapper = as.data.frame(MultiAssayExperiment::sampleMap(mae)) %>% filter(assay=='RNAseq_rnaaccess')
  colnames(cds) = mapper$primary[match(colnames(cds),mapper$colname)]
  cds.ensembl <- cds 
  counts = assay(mae[['RNAseq_rnaaccess']])
  annot = rowData(mae[['RNAseq_rnaaccess']])
  rownames(cds) = annot$HGNC[match(rownames(cds),annot$GeneID)]
  
  # Handling excaption
  toDropDupl = duplicated(rownames(cds))
  cds.hgnc = cds[!toDropDupl,]
  
  # Output: two nomenclatures, Ensembl format or HGNC format
  results = list("cds.emsembl" = cds.ensembl, "cds.hgnc" = cds.hgnc)
  return(results) 
}

cds2 <- cdsAnnotGse(gseCLL14RNAqc.all, mae_cll14)$cds.hgnc # 25008 genes
cds3 <- cdsAnnotGse(gseCLL14RNAqc.exp, mae_cll14)$cds.hgnc # 11234 genes

# cds2 was used to create the transcriptomic SNN clustering in the manuscript.
# cds3 was subjected to Limma/Voom DE analysis in following section. 




#########################################
############# Limma/Voom ################
#########################################

# Helper 1: Design matrix for MRD, IgHV and treatment
design_MrdIgvh_trt <- function(clinData, col, matchClin) {
  mrd<- clinData[[col]][matchClin] %>% ifelse(is.na(.), "MRDmissing", .) 
  igvh <- as.character(clinData$IGVH[matchClin]) 
  igvh <- ifelse(is.na(igvh), "IGHVmissing", ifelse(igvh =="NOT EVALUABLE", "missing", igvh)) 
  group <- paste0(mrd, ".", igvh)
  group <- factor(group)
  trt <- clinData$trtArm[matchClin]
  
  design <- model.matrix(~0 + group + trt)
  colnames(design) <- c(levels(group), "OBIN_CHLO", "OBIN_GDC") 
  design
} 

# Helper 2: Limma-Voom - log2-CPM ready for linear modeling
normalizeCDS <- function(cds, design){
  libSizes <- sizeFactors( cds ) * sum(counts(cds)/ncol(cds)) 
  vm <- voom( counts( cds ), design, libSizes, plot=F )
  
  return(vm)
}

# Helper 3: Limma-Voom - LmFit and eBayes for DE analysis
linModelSigGene <- function(vm, design, contrast, lfc = 0.5, pvalue = 0.05, coef= c(1), plotMD=TRUE) {
  vfit = lmFit(vm, design) %>% 
    contrasts.fit(., contrasts = contrast) 
  # Compute moderated t-statistics of DE by empirical Bayes moderation of standard errors
  efit = eBayes(vfit) 
  # Handling exceptions
  efit = efit[!is.na(rownames(efit)),] 
  
  # topTable
  tt = list()
  for (i in coef) {
    coef = topTable(efit, coef=i, number=Inf, sort.by="none") 
    tt[[paste0("coef",i)]] = coef 
  }
  # Moderated T-test
  tfit = treat(vfit, lfc=lfc) 
  dt = decideTests(tfit, p.value=pvalue, adjust.method = "none")
  # MA plot as sanity check
  if (isTRUE(plotMD)) {
    for (i in seq(ncol(tfit$contrasts))) {
      plotMD(tfit, column=i, status=dt[,i], main=colnames(tfit)[i], xlim=c(-1,16))
    }
  }
  
  # Output
  common <- Reduce(intersect, apply(dt, 2, function(x) which(x !=0))) 
  common.genes <- rownames(tfit$t)[common] 
  
  all <- Reduce(union, apply(dt, 2, function(x) which(x !=0))) 
  all.genes <- rownames(tfit$t)[all]
  
  results = list("efit" = efit, "tt" = tt, "tfit" = tfit, "dt" = dt, "cg" = common.genes, "ag" = all.genes)
  return(results)
  
}


# Execution
design <- design_MrdIgvh_trt(clinData, "FUM3BLNGS.GS.2", matchClin) 
cont.matrix <- makeContrasts(HLvsU6_ALL = (HL_MRD.MUTATED + HL_MRD.UNMUTATED)/2- (uMRD6.MUTATED + uMRD6.UNMUTATED)/2, #contrast 1
                             HLvsU6_M = HL_MRD.MUTATED - uMRD6.MUTATED, #contrast 2
                             HLvsU6_U = HL_MRD.UNMUTATED - uMRD6.UNMUTATED, #contrast 3
                             levels=design)
vm3.eot = normalizeCDS(cds3[,matchCDS], design)


lmVm3EOT = linModelSigGene(vm3.eot, design, cont.matrix, 
                           lfc = 0.5, #tfit 
                           pvalue = 0.05, #decideTest, dt
                           coef = c(1,2,3), #contrast
                           plotMD=FALSE) 

# All tt files
tt.eot.mu <- lmVm3EOT$tt$coef1 %>% rownames_to_column("ID") %>% as_tibble # MCLL UCLL agonistic analysis

# Note: vm3 object was used to investigate patient level log2-CPM and GSVA. The final tt.eot.mu was used for volcano plot




###############################################################
######################## Part 2: Analysis #####################
###############################################################

############################################
################### fGSEA ##################
############################################

# Helper function: fGSEA and leading edge gene selection
FGSEAtable = function(pathways, stats, pathList, outputPath, name_contrast) {
  
  fgseaRes = fgsea(pathways=pathways, stats=stats, eps=1e-10,minSize=15, maxSize=500) 
  fgseaRes = fgseaRes %>% mutate(numLeadingEdge = lengths(leadingEdge)) 
  
  # Write all results in csv for downstream analysis
  write.csv(fgseaRes[, -8], file=paste0(outputPath,"/fGSEA_sig_",name_contrast,'_All.csv'), row.names=F)
  
  # Summary table (Show only pathways whose number of leading edge gene is more than 5)
  col2show=c('pathway','padj','ES','NES', 'size')
  topPathwaysUp <- fgseaRes[ES > 0 & lengths(leadingEdge) >= 5][head(order(padj), n=20), pathway] 
  df.up <- fgseaRes[fgseaRes$pathway %in% topPathwaysUp, col2show, with=FALSE] %>% arrange(padj)
  topPathwaysDown <- fgseaRes[ES < 0 & lengths(leadingEdge) >= 5][head(order(padj), n=20), pathway] 
  df.down <- fgseaRes[fgseaRes$pathway %in% topPathwaysDown, col2show, with=FALSE] %>% arrange(padj)
  
  main.title <- "fGSEA Top genesets in each direction"
  subtitle <- paste0("Input: t-statistics, msigDB genesets: ", paste(pathList,collapse=", "), 
                     " Excluded genesets whose min# of leadingEdge gene is 5, Comparison contrast is: ", name_contrast) %>%
    strwrap(width = 80) %>%
    paste(collapse = "\n")
  
  df.combined <- rbind(df.up, df.down) %>%
    mutate(across(where(is.numeric), ~ round(., 4))) #display only 3 digit
  
  t <- ggtexttable(df.combined, rows = NULL, theme = ttheme("blank"))%>%
    tab_add_hline(at.row = c(1:2, 22), row.side = "top", linewidth = 2) %>% 
    tab_add_hline(at.row = c(nrow(df.combined)+1), row.side = "bottom", linewidth = 2) %>% 
    tab_add_title(text = subtitle, face = "plain", size = 10) %>%
    tab_add_title(text = main.title, face = "bold", padding = unit(0.1, "line")) %>%
    tab_add_footnote(text = "*Displayed top 20 pathways in each direction", size = 10, face = "italic")
  
  ggsave(t, width = 15, height = 14,
         file=paste0(outputPath,"/fGSEA_sig_",name_contrast,'_top20_table.png'))
  
  
  # Write sig hit defined by FDR cut off
  FDRCutoff=.1
  fgseaResSig <- subset(fgseaRes, padj <= FDRCutoff)
  if (nrow(fgseaResSig) > 0) {
    write.csv(fgseaResSig[, -8], file=paste0(outputPath,"/fGSEA_sig_",name_contrast,'_FDRcutoff', FDRCutoff, '.csv'), row.names=F)
  } else {
    cat('No significant pathways were found  \n')
  }
  
  
  # Leading edge gene analysis
  leadingUp <- fgseaRes[fgseaRes$pathway %in% topPathwaysUp, c('pathway', "leadingEdge", "padj"), with=FALSE] %>% arrange(padj)
  leadingDown <- fgseaRes[fgseaRes$pathway %in% topPathwaysDown, c('pathway', "leadingEdge", "padj"), with=FALSE] %>% arrange(padj)
  
  # Helper function to create leading df
  leadingDF = function(df) {
    l <- df[["leadingEdge"]]
    l <- setNames(l, df[["pathway"]])
    leading_all <- unique(unlist(l,use.names=FALSE))
    
    ls = list()
    for (i in seq(length(names(l)))) {
      ls[[i]]<-ifelse(is.na(match(leading_all, l[[i]])), 0,1)
    }
    names(ls) <- names(l)
    
    df <- data.frame(Reduce(rbind, ls))
    colnames(df) <- leading_all
    rownames(df) <- names(l)
    
    df.ordered <- df[, order(colSums(-df, na.rm=T))]#align low based on high colsum
    df.ordered
  }
  
  upregulated <- leadingDF(leadingUp)
  downregulated <- leadingDF(leadingDown)
  
  
  # Helper function to create heatmap from leading df (fGSEA)
  leadingHeatMap = function(df, outputPath) {
    n = ifelse (dim(t(df))[1]>200, 200, dim(t(df))[1])
    ph = pheatmap(t(df)[1:n,], 
                  cutree_rows = n%/%20, #if 100 genes, 5 cluster 
                  main = paste("Top commonly enriched leading edge genes\n from the top genesets of", deparse(substitute(df)), "pathways in",paste(msig_list,collapse="_")))
    ggsave(ph, width = 10, height = n%/%6.5, #integer devision
           file=paste0(outputPath,"/fGSEA_leading_edge_genes_top_",deparse(substitute(df)),"_",name_contrast,'.png'))
    write.csv(t(df), file=paste0(outputPath,"/fGSEA_leading_edge_genes_top_",deparse(substitute(df)),"_",name_contrast,'.csv'), row.names=T)
  }
  
  
  leadingHeatMap(upregulated, outputPath)
  leadingHeatMap(downregulated, outputPath)
}


# Hallmark geneset
msig_list = c("H") 
genesets_select <-msig.hs[msig.hs$gs_cat %in% msig_list, ] 
msigdbr_gs_list = split(x = genesets_select$gene_symbol, f = genesets_select$gs_name)

# Contrast selection and stats
select_contr = 1 # All patients
name_contrast = colnames(lmVm3EOT$dt[,select_contr]) 
ranks = sort(lmVm3EOT$efit$t[,select_contr], decreasing = FALSE) 

# Output path
outputPath = paste0(file.path(homeDir, "exports", "fGSEA_EOT_MU"), "_", paste(msig_list,collapse="_")) 
dir.create(outputPath, showWarnings = FALSE) #will be ignored if the folder already exist

# Run fGSEA and create table, export sig and plots in output path
FGSEAtable(pathways=msigdbr_gs_list, stats=ranks, pathList=msig_list, 
           outputPath, name_contrast)

# Note: The exported leading edge genes are subjected to downstream analysis including heatmap visualization



############################################
################### GSVA ###################
############################################

# Example from baseline vs. Progression comparison
# Selecting gene sets to investigate (Hallmark)
gs_ven_resistance <- c("HALLMARK_IL2_STAT5_SIGNALING", "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_INFLAMMATORY_RESPONSE")
annot_top_gs <- msigdbr_gs_list[gs_ven_resistance]
annot_top_gs<- annot_top_gs[order(names(annot_top_gs))]
names(annot_top_gs)

# Running gsva using vm3
gsva.es.allgenes.prog <- gsva(vm3.prog[["E"]], annot_top_gs, verbose = FALSE)

# Clinical data comparing baseline vs progression sample 
clinData.Prog.gsva <-  clinData.prog44 %>% 
  filter(trtArm != "NOTASSGN") %>% 
  filter(!is.na(seq_TIME_POINT)) %>% 
  mutate(trtArm = factor(trtArm, 
                         levels = c("OBIN_CHLO", "OBIN_GDC"),
                         labels = c("Clb-Obi", "Ven-Obi"))) %>% 
  mutate(seq_TIME_POINT = factor(seq_TIME_POINT, 
                                 levels = c("Baseline", "Progression")))

# Visualization
fig.GSVA.prog.parentgs <- left_join(clinData.Prog.gsva %>% dplyr::select(uID_time, seq_TIME_POINT, trtArm), t(gsva.es.allgenes.prog) %>% as_tibble(rownames = "uID_time"), by="uID_time") %>% 
  gather(pathway, gsva_score, 4:6) %>% 
  filter(!is.na(seq_TIME_POINT)) %>% 
  mutate(pathway = factor(pathway, 
                          levels = c("HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE", 
                                     "HALLMARK_IL2_STAT5_SIGNALING"),
                          labels = c("HALLMARK\nINFLAMMATORY\nRESPONSE", "HALLMARK\nINTERFERON\nGAMMA\nRESPONSE", 
                                     "HALLMARK\nIL2_STAT5\nSIGNALING"))) %>% 
  ggplot(., aes(x=seq_TIME_POINT, y=gsva_score)) +
  geom_boxplot(aes(fill=seq_TIME_POINT)) +
  geom_signif(test="wilcox.test", comparisons = list(c("Baseline", "Progression"))) +
  scale_fill_manual(values=c("purple", "#69b3a2")) +
  ylab("GSVA score distribution\n from Hallmark gene sets") +
  ylim(NA, 0.8)+
  facet_grid(trtArm~pathway, labeller = label_value) +
  theme_bw() +
  theme(legend.position = "none",
        strip.background = element_rect(colour="black",fill="white"), 
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=12, color="black")) 

fig.GSVA.prog.parentgs



############################################
################# Heatmap ##################
############################################

# Leading edge (LE) gene heatmap analysis example from baseline vs. prog comparison
# Run fGSEA 
msig_list = c("H") 
genesets_select <-msig.hs[msig.hs$gs_cat %in% msig_list, ] 
msigdbr_gs_list = split(x = genesets_select$gene_symbol, f = genesets_select$gs_name)
select_contr = 1 
name_contrast = colnames(lmVm3PROG_Alltrt$dt[,select_contr]) 
ranks = sort(lmVm3PROG_Alltrt$efit$t[,select_contr], decreasing = FALSE) 
fgseaRes = fgsea(pathways=msigdbr_gs_list, stats=ranks, eps=1e-10,minSize=15, maxSize=500) 

# Extract the LE genes from GSEA results
le <- fgseaRes %>% unnest(leadingEdge) %>% dplyr::select(pathway, leadingEdge)
g2m <-le[le$pathway == "HALLMARK_G2M_CHECKPOINT",]%>% dplyr::select(leadingEdge) %>% pull
tnfa <- le[le$pathway == "HALLMARK_TNFA_SIGNALING_VIA_NFKB",] %>% dplyr::select(leadingEdge) %>% pull
myc_v1 <-le[le$pathway == "HALLMARK_MYC_TARGETS_V1",]%>% dplyr::select(leadingEdge) %>% pull

# Assign category and ordering table
tt.eot.mu.alltrt.prog.lab <- tt.eot.mu.alltrt.prog %>% mutate(cat = ifelse(ID %in% g2m, "HALLMARK_G2M_CHECKPOINT", 
                                                                           ifelse(  ID %in% tnfa, "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
                                                                                    ifelse(ID %in% myc_v1,"HALLMARK_MYC_TARGETS_V1", "Others")))) %>% 
  mutate(cat = factor(cat, levels = c("HALLMARK_G2M_CHECKPOINT","HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_MYC_TARGETS_V1", "Others"))) 

# Helper function for heatmap
HeatmapBySigPROG_ordered <- function(vm, df_top, clinData_subset, col_order) {
  mat = vm[["E"]][rownames(vm[["E"]]) %in% unique(df_top[["ID"]]),]
  my_gene_col = df_top %>% dplyr::select(ID, cat) %>% group_by(ID) %>%  slice(1) %>%  column_to_rownames("ID")
  my_sample_col = clinData_subset %>% 
    dplyr::select(uID_time, trtArm, seq_TIME_POINT, IGVH) %>% 
    remove_rownames() %>%  
    column_to_rownames("uID_time")
  # Color setting
  my_colour = list(
    cat = c(`HALLMARK_TNFA_SIGNALING_VIA_NFKB` = "#d1495b", `HALLMARK_G2M_CHECKPOINT`= "#00798c", `HALLMARK_MYC_TARGETS_V1`= "#edae49"), 
    trtArm = c(OBIN_CHLO = "#CCCCCC", OBIN_GDC = "#333333", NOTASSGN = "#FFFFFF"), 
    seq_TIME_POINT = c(Baseline = "#F4EDCA", Progression = "#C4961A"), 
    IGVH = c(UNMUTATED = "#B95FBB", MUTATED = "#D4D915", `NOT EVALUABLE` = "#FFFFFF")
  )
  # Pheatmap using z-score from log2-cpm
  cal_z_score <- function(x){(x - mean(x)) / sd(x)}
  mat_z <- t(apply(mat, 1, cal_z_score))
  pheatmap(mat_z[,col_order], 
           annotation_colors = my_colour,
           annotation_row = my_gene_col, 
           annotation_col = my_sample_col,
           cluster_cols = FALSE, 
           cutree_rows = 1,
           cutree_cols = 2,
           main = "z-score derived from Log2-CPM")
}

# Setup
clinSubset = clinData.prog44[matchClinPROG,]
vm3Subset = vm3.prog 
df_top = tt.eot.mu.alltrt.prog.lab %>% 
  filter(!cat %in% c("Others")) %>% 
  filter(logFC >0.5) %>% 
  filter(P.Value <0.01) # Top 45 genes

# Order patients to separate baseline vs prog
order <- clinSubset %>% arrange(seq_TIME_POINT, trtArm, IGVH) %>% dplyr::select(uID_time) %>% pull

# Heatmap result
h.order <- HeatmapBySigPROG_ordered(vm3Subset,
                                    df_top,
                                    clinSubset,
                                    order)
h.order



############################################
######## Transcriptomic Clustering #########
############################################

# Helper function
SeuratUMAP <- function(cds) {
  set.seed("0306")
  data_umap <- CreateSeuratObject(counts = as.data.frame(assays(cds)[["tpm"]])) 
  
  all.genes <- rownames(data_umap)
  data_umap <- ScaleData(data_umap, features = all.genes)
  data_umap <- FindVariableFeatures(object = data_umap)
  seurat_umap <- RunPCA(data_umap, features = VariableFeatures(object = data_umap))
  
  seurat_umap <- JackStraw(seurat_umap, num.replicate = 100)
  seurat_umap <- ScoreJackStraw(seurat_umap, dims = 1:15)
  
  seurat_umap <- FindNeighbors(seurat_umap, dims = 1:10) 
  seurat_umap <- FindClusters(seurat_umap, resolution = 0.6) 
  seurat_umap = RunUMAP(seurat_umap, dims = 1:10)
  return(seurat_umap)
}

# Run clustering using cds2 object 
umap.cds2 <- SeuratUMAP(cds2) 

# Reorder cds2 to match clinical data 
matchToCDS2 <- match(rownames(umap.cds2@meta.data), clinData$uID) 

cat.1 <- c("IGVH", "TP53", "NOTCH1", "SF3B1", "ATM", "XPO1", "RPS15", "POT1", "NFKBIE", "BIRC3", "EGR2", "MYD88") 
cat.2 <- c("Q11", "T12", "Q13C", "P17", "CYTOHIER") 
cat.3 <- c("Sex", "AgeGR", "trtArm", "EOTBLPC3", "EOTBMPC3", "FUM3BLNGS", "FUM3BLNGS.2", "FUM3BLNGS.GS", "FUM3BLNGS.GS.2") 
cat.4 <- c("PFSINV_AVAL_M_5YR", "OS_AVAL_M_5YR") 
cat.all <- c(cat.1, cat.2, cat.3, cat.4)

# Assign clindata into umap meta data for visualization
for (i in cat.all) {
  umap.cds2@meta.data[[i]] = clinData[matchToCDS2,i] 
}

# Update labels for visualization
umap.cds2@meta.data<- umap.cds2@meta.data %>% 
  mutate(trtArm = factor(trtArm, 
                         levels = c("OBIN_CHLO", "OBIN_GDC"),
                         labels = c("Clb-Obi", "Ven-Obi"))) %>% 
  mutate(FUM3BLNGS.GS.2 = factor(FUM3BLNGS.GS.2, # Replced with German Study Group's MRD 
                                 levels = c("HL_MRD", "uMRD45", "uMRD6"),
                                 labels = c("MRD+", "MRD4/5", "MRD6"))) %>% 
  rename(IgHV = IGVH, `Del 17p`=P17, `Trisomy 12`=T12, `Del 13q`=Q13C,
         `NGS MRD at EOT`=FUM3BLNGS.GS.2, `Treatment Arm` = trtArm)

var <- c("Treatment Arm", "NGS MRD at EOT", "IgHV", "TP53", "Del 17p", "NOTCH1", "SF3B1", "Trisomy 12", "Del 13q")  # Var of interests
var2 <- c("Sex" ,"AgeGR", "ATM", "XPO1", "RPS15", "POT1", "NFKBIE", "BIRC3", "EGR2", "MYD88", "Q11") # Var category 2
var.all <- c(var, var2)

for (val in var.all) {
  if (length(unique(umap.cds2@meta.data[[val]])) == 3) {
    results <- DimPlot(umap.cds2, reduction = "umap",group.by = val,
                       cols = c("#B95FBB", "#D4D915", "#aeadb3")) +
      theme(plot.title = element_text(size = 12, face = "bold"), legend.position = "bottom", legend.text=element_text(size=8)) 
    assign(val, results)
  } else if (length(unique(umap.cds2@meta.data[[val]])) == 4) {
    results <- DimPlot(umap.cds2, reduction = "umap",group.by = val,
                       cols = c("#B95FBB", "#D4D915", "#333333", "#aeadb3")) +
      theme(plot.title = element_text(size = 12, face = "bold"), legend.position = "bottom", legend.text=element_text(size=8)) 
    assign(val, results)
  } else if (length(unique(umap.cds2@meta.data[[val]])) == 2) {
    results <- DimPlot(umap.cds2, reduction = "umap",group.by = val,
                       cols = c("#B95FBB", "#D4D915")) +
      theme(plot.title = element_text(size = 12, face = "bold"), legend.position = "bottom", legend.text=element_text(size=8)) 
    assign(val, results)
  } else {print("The numbers of variables are more than 4")}
} 


# Visualize selected biomarkers
fig = annotate_figure(ggarrange(`Treatment Arm`, `NGS MRD at EOT`, IgHV, TP53, `Del 17p`, NOTCH1, SF3B1, `Trisomy 12`, `Del 13q`, ncol = 3, nrow = 3), 
                      top = text_grob("Transcriptomic SNN clustering and UMAP visualization of key markers", color = "Black", face = "bold", size = 14),
                      bottom = text_grob(paste0("Gene filtering criteria: 3 counts min of 5% patients, X,Y, mitocondrial and LINC removed", 
                                                     "\nTotal of", nrow(umap.cds2) -1, " genes were used"), 
                                              color = "grey27",
                                              hjust = 1, x = 1, face = "italic", size = 10))
fig

