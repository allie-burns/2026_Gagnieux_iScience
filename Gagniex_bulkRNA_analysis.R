#################################################################################
## Set up environment
#################################################################################
library(tidyverse) ## version 2.0.0
library(DESeq2) ## version 1.40.2
library(gprofiler2) ## version 0.2.2
library(ggVennDiagram) ## version 1.5.4
library(eulerr) ## version 7.0.2
library(writexl) ## version 1.5.4

## Get Sample Information
sampInfo <- readRDS("./data/SampleInfo.rds")

## Define Paths
alignPath <- "./data/1_alignment/"
dePath <- "./data/2_DifferentialExpression/"
pkdSig <- "./data/3_PKDsignature/"

## Define parameters
fc <- 0
fdr <-  0.05
org <- "mmusculus"  

#################################################################################
## Alignment Quality Control
#################################################################################
## Load STAR log data
log.fls <- list.files(paste0(alignPath,"logFinalout"), full.names = TRUE)
names(log.fls) <- gsub("_Log.final.out","",basename(log.fls))  ## name files
names(log.fls) <- sub("_[^_]+$", "", names(log.fls))

## Read and reformat stats information
getStats <- function(fls) {
    ## Read file 
    stats <- read.delim(fls,
                        header = FALSE,
                        col.names = c("variable","value"))
    ## Organize statistics of interest into a dataframe
    data <- data.frame(total_reads = stats[5,2],
                       avg_read_length = stats[6,2],
                       uniq_mapped = stats[8,2],
                       uniq_mapped_perc = stats[9,2],
                       avg_mapped_length = stats[10,2],
                       multi_mapped_loci = stats[23,2],
                       multi_mapped_loci_perc = stats[24,2],
                       multi_mapped_too_many_loci = stats[25,2],
                       multi_mapped_too_many_loci_perc = stats[26,2],
                       unmapped_mismatch= stats[28,2],
                       unmapped_mismatch_perc= stats[29,2],
                       unmapped_short = stats[30,2],
                       unmapped_short_perc= stats[31,2],
                       unmapped_other = stats[32,2],
                       unmapped_other_perc= stats[33,2]
                       )
}

stats <- lapply(log.fls, getStats) ## Run function on each file
stats <- do.call(rbind,stats)  ## create matrix of stats for all files

## Reformat read counts for plotting
raw <- stats |>
    select(!grep("perc|avg|unmapped", colnames(stats))) |>
    rownames_to_column( var = "sample") |>
    gather(key = "stat", value = "value", -sample) |>
    mutate(value = as.numeric(value),
           stat  = factor(stat, levels = unique(stat))) |>
    mutate(treat = sampInfo$treat[match(sample, sampInfo$sample)],
           sex   = sampInfo$sex[match(sample, sampInfo$sample)],
           a1    = sampInfo$animal_n1[match(sample, sampInfo$sample)],
           a2    = sampInfo$animal_n2[match(sample, sampInfo$sample)]) |>
    mutate(label = paste(treat,sex,a1,a2, sep = "_"))

## Plot read counts
alnStats.raw <-
    ggplot(raw, aes(x = label, y = value, fill = treat, pattern = sex)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap(~stat, scales = "free_y",nrow = 1) +
    theme_bw() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

## Plot percent of counts 
perc <- stats |>
    select(grep("perc", colnames(stats))) |>
    rownames_to_column( var = "sample") |>
    gather(key = "stat", value = "value", -sample) |>
    mutate(perc  = as.numeric(gsub ("%","",value)),
           stat  = factor(stat, levels = unique(stat))) |>
    mutate(treat = sampInfo$treat[match(sample, sampInfo$sample)],
           sex   = sampInfo$sex[match(sample, sampInfo$sample)],
           a1    = sampInfo$animal_n1[match(sample, sampInfo$sample)],
           a2    = sampInfo$animal_n2[match(sample, sampInfo$sample)]) |>
    mutate(label = paste(treat,sex,a1,a2, sep = "_"))

alnStats.perc <-
    ggplot(perc, aes(x = label, y = perc, fill = stat)) +
    geom_bar(stat = "identity", position = "stack") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.title.x = element_blank()) +
    geom_text(data = subset(perc, perc > 1),
              aes(label = perc),
              position = "stack", vjust=1, size=3) +
    labs(y = "% total reads")

#################################################################################
## Get counts matrix 
#################################################################################
## Get files
counts.fls <- list.files(paste0(alignPath,"ReadsPerGene"), full.names = TRUE)
names(counts.fls) <- gsub("_ReadsPerGene.out.tab","",basename(counts.fls))
names(counts.fls) <- sub("_[^_]+$", "", names(counts.fls))

## Read counts files
counts <- lapply(counts.fls, function(fls) {
    x <- read.table(fls) 
    gns <- x$V1
    x <- x$V4
    names(x) <- gns
    x
})
counts <- do.call(cbind,counts)
counts <- counts[-c(1:4),]

#################################################################################
## Differential Expression Analysis
#################################################################################
## Get gene names and descriptions for all genes in counts table
gnInfo <- gconvert(rownames(counts), organism = "mmusculus")

## Order counts to match sample info
counts <- counts[,match(sampInfo$sample,colnames(counts))]

## Create and filter DEseq object
deseq <- DESeqDataSetFromMatrix(countData = counts,
                                colData = sampInfo,
                                design = ~group)

keep <- rowSums(counts(deseq) > 10 ) >= 4
deseq <- deseq[keep,]

## Run PCA analysis
vsd <- vst(deseq, blind = FALSE) ## Variance stabilizing transformation 
pca.plot <-
    plotPCA(vsd, intgroup=c("treat", "sex"), returnData = FALSE ) +
    theme_bw() +
    geom_text(aes(label = colnames(vsd)), position = position_nudge(y = 1), size = 2.5) + 
    labs(title = "Principal Component Analysis") +
    xlim(-25,25) +
    ylim(-20,20)

## Run DEanalysis
deseq <- DESeq(deseq)
contrasts <- list(cKO_ctrl_male = c("group", "cKO.KSP_MALE","CTRL_MALE"),
                  cKO_ctrl_fem  = c("group", "cKO.KSP_FEM","CTRL_FEM"))
desc <- gconvert(rownames(deseq), organism = "mmusculus")
results <- lapply(contrasts, function(con) {
    ## Get contrast of interest
    r <- results(deseq, contrast = c(con[1],con[2],con[3]))
    ## Format results
    r <- r[order(r$padj),] ## order by p.value
    r <- rownames_to_column(data.frame(r), var = "ensembl_id") |> ## add ens id
        ## Add gene name and gene description
        mutate(gene_name = desc$name[match(gsub("\\..*","",ensembl_id),desc$input)],
               description = desc$description[match(gsub("\\..*","",ensembl_id),desc$input)])  
})

## Volcano Plots
getVolcanos <- function(i, results) {
    comp <- names(results[i])
    res <- data.frame(results[[i]])
    ## Assign Directionality
    res$diffexpressed <- "NO"
    res$diffexpressed[res$log2FoldChange >= fc & res$padj <= fdr] <- "up"
    res$diffexpressed[res$log2FoldChange <= -fc & res$padj <= fdr] <- "dw"
    ## Colors
    cols <- c(up = "#FF0000", dw = "#63B2B3", NO = "#A6A6A6")
    ## Build Plot
    ggplot(data=res,
           aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label = gene_name)) +
        geom_point() +
        scale_color_manual(values = cols) +
        theme_classic() +
        theme(panel.border = element_rect(colour = "black", fill=NA, size=1.5),
              text = element_text(size=14, family = "Helvetica"),
              axis.ticks.y = element_line(size = 0.8, color = "black"),
              axis.text.y = element_text(size = 12,color = "black"),
              axis.text.x = element_text(size = 12,color = "black"),
              legend.position = "none",
              legend.background = element_blank()) +
        labs(title = comp,
             x = expression("log"[2]*'FC') ,
             y = expression("-log"[10]*italic("(P")*" value)")) +
        geom_vline(xintercept=c(-fc, fc), col="grey") +
        guides(color = guide_legend(override.aes = list(size=3))) +
        xlim(-10,10) +
        geom_text(aes(label=ifelse(abs(log2FoldChange) >= 2,
                                   as.character(gene_name),'')), vjust = -1) + 
        annotate(geom = "text",
                 x = -9,##min(res$log2FoldChange),
                 y = max(-log10(na.omit(res$pvalue))), 
                 label = paste("n =",sum(na.omit(res)$padj <= fdr &
                                                     na.omit(res)$log2FoldChange < -fc)),
                 hjust = 0.2, vjust = 3, size = 4) +
        annotate(geom = "text",
                 x = 7, ##max(res$log2FoldChange) ,
                 y = max(-log10(na.omit(res$pvalue))), 
                 label = paste("n =",sum(na.omit(res)$padj <= fdr &
                                                     na.omit(res)$log2FoldChange > fc)),
                 hjust = 0.2, vjust = 3, size = 4)
}

## Plot Volcano Plots
volPlots <- lapply(seq(1:length(results)), getVolcanos, results)
volcano.plots <- cowplot::plot_grid(plotlist = volPlots, nrow = 1)

#################################################################################
## Define genes that are up or down regulated in both males and females (overlap)
#################################################################################
## Get Gene Lists
up.list <- lapply(results, function(res) {
    r <- na.omit(res$gene_name[res$log2FoldChange > fc & res$padj <= fdr])
    r[r != "None"]##  r[-grepl("None", r)]
})

dw.list <- lapply(results, function(res) {
    r <- na.omit(res$gene_name[res$log2FoldChange < -fc & res$padj <= fdr])
    r[r != "None"]##  r[-grepl("None", r)]
})

## Make Venn Diagrams
param <- theme(legend.position = "none")
sexOverlap.venn <- 
    cowplot::plot_grid(
                 plotlist = list(
                     ggVennDiagram(c(up.list[1], up.list[2])) + labs (title = "up_male - up_female") + param,
                     ggVennDiagram(c(dw.list[1], dw.list[2])) + labs (title = "dw_male - dw_female") + param),
                 nrow = 1)

#################################################################################
## Ontology analysis for up/dw in both males and females
#################################################################################
## Get overlaps of interest
up.both <- intersect(up.list[[1]], up.list[[2]])
dw.both <- intersect(dw.list[[1]], dw.list[[2]])

## Create overlap tables
formatOverlaps <- function(gn.list) {
    ol <- data.frame(gene_name = gn.list)
    ol$ensembl_id <- results[[1]]$ensembl_id[match(ol$gene_name, results[[1]]$gene_name)]
    ol$log2FC_male <- results[[1]]$log2FoldChange[match(ol$gene_name, results[[1]]$gene_name)]
    ol$padj_male <- results[[1]]$padj[match(ol$gene_name, results[[1]]$gene_name)]
    ol$log2FC_female <- results[[2]]$log2FoldChange[match(ol$gene_name, results[[2]]$gene_name)]
    ol$padj_female <- results[[2]]$padj[match(ol$gene_name, results[[2]]$gene_name)]
    ol
}

up.in.both <- formatOverlaps(up.both)
dw.in.both <- formatOverlaps(dw.both)

## Setup functions
getGost <- function(x) {
    ## Define list of genes for over-representation analysis
    query <- na.omit(c(x$gene_name))
    if(length(query) == 0 ){ query <- "none"}
    ## Define list of genes for comparison
    universe <- as.character(na.omit(unique(c(results[[1]]$gene_name, results[[2]]$gene_name))))
    ## Run Gene Enrichment Analysis
    myGost <- gost(query = query,             
                   organism = org,
                   ordered_query = FALSE,      
                   significant = TRUE,        
                   user_threshold = 0.05,     
                   correction_method = "fdr", 
                   sources = c("GO","KEGG"),  
                   evcodes = TRUE,            
                   domain_scope = "custom",   
                   custom_bg = universe)      
    x <- myGost
    if(is.null(x$result)) {
        data.frame(source = as.character(),
                   p_value = as.numeric(),
                   term_name = as.character(),
                   term_size = as.numeric(),
                   query_size = as.numeric(),
                   intersection_size = as.numeric(),
                   intersection = as.character())
    }else{
        x$result |>
            select(source, p_value, term_name, term_size,
                   query_size, intersection_size, intersection) |>
            filter(p_value <= 0.05) |>
            mutate(p_value = round(p_value, 3)) 
    }
}

## Run ontology analysis
up.gost <- getGost(up.in.both)
dw.gost <- getGost(dw.in.both)

#################################################################################
## Comparison to Malas data - overlap Venn Diagram
#################################################################################
## Add entrez id to DE results tables for matching to PKD signautre
entrez <- gconvert(query = results[[1]]$ensembl_id,
                   organism = "mmusculus",
                   target = "ENTREZGENE_ACC")
results[[1]]$entrez <- entrez$target[match(results[[1]]$ensembl_id, entrez$input)]
results[[2]]$entrez <- entrez$target[match(results[[2]]$ensembl_id, entrez$input)]
results[[1]] <- results[[1]][match(results[[2]]$ensembl_id, results[[1]]$ensembl_id),]

## PKD DEG information (Malas = pval < 0.05)
malasPath <- "./data/2017_malas/Data1_DEGs.xlsx"
malasSheet <- readxl::excel_sheets(malasPath)
malasSheet <- malasSheet[!grepl("Key",malasSheet)]
malasDEG <- readxl::read_xlsx(malasPath, sheet = malasSheet[[1]], skip = 5)
malasDEG <- janitor::clean_names(malasDEG)

## load Meta PKD analysis (Malas = pval < 0.05)
metaPath <- "./data/2017_malas/Data2_meta.xlsx"
metaSheets <- readxl::excel_sheets(metaPath)
metaSheets <- metaSheets[!grepl("Key",metaSheets)]
pkdMeta <- lapply(metaSheets, function(x) { readxl::read_xlsx(metaPath, sheet = x, skip = 3) })
names(pkdMeta) <- metaSheets
pkdMeta <- lapply(pkdMeta, function(x) { janitor::clean_names(x) })
pkdMeta <- lapply(pkdMeta, function(x) { x[-grep("STX17", x$gene_symbol),] }) ## Remove STX17 (on both up and down lists)

## Load Injury Repear Programs
injuryPath <- "./data/2017_malas/Data4_InjuryRepairProfiles.xlsx"
injurySheet <- readxl::excel_sheets(injuryPath)
injurySheet <- injurySheet[!grepl("Key",injurySheet)]
injuryDEG <- list(
    exp_early = readxl::read_xlsx(injuryPath, sheet = injurySheet[[1]], skip = 6),
    exp_late = readxl::read_xlsx(injuryPath, sheet = injurySheet[[2]], skip = 5),
    lit = readxl::read_xlsx(injuryPath, sheet = injurySheet[[3]], skip = 3)
)
injuryDEG <- lapply(injuryDEG, function(x) { janitor::clean_names(x) })
injuryDEG <- injuryDEG[1:2]


## Get up and down lists separate - merge male /female - must be sig in both
allCKO <- tolower(unique(na.omit(
    ifelse(results[[1]]$padj <= fdr & results[[2]]$padj <= fdr,
           results[[1]]$gene_name,
           NA)
)))

## Up in both male and female
upCKO <- tolower(unique(
    na.omit(
        ifelse(results[[1]]$padj <= fdr & results[[1]]$log2FoldChange > fc &
               results[[2]]$padj <= fdr & results[[2]]$log2FoldChange > fc ,
               results[[1]]$gene_name,
               NA)
    )))

## Down in both male and female
dwCKO <- tolower(unique(
    na.omit(
        ifelse(results[[1]]$padj <= fdr & results[[1]]$log2FoldChange < fc &
               results[[2]]$padj <= fdr & results[[2]]$log2FoldChange < fc ,
               results[[1]]$gene_name,
               NA)
    )))

## PKD signature (up and down lists separate)
allPKD <- tolower(unique(c(pkdMeta[[1]]$gene_symbol,pkdMeta[[2]]$gene_symbol )))
upPKD <-  tolower(unique(c(pkdMeta[[1]]$gene_symbol)))
dwPKD <-  tolower(unique(c(pkdMeta[[2]]$gene_symbol)))

## Injury Responese signature (early and late separate, up and dw separate)
allINJ <- tolower(unique(c(injuryDEG[[1]]$associated_gene_name,
                           injuryDEG[[2]]$associated_gene_name)))

## Venn Diagram -  setup gene lists for venns
all <- list("Experimental\nInjury Repair Signature" = unique(allINJ),
            "PKD signature" = unique(allPKD),
            "Significant\ncKO genes" = unique(allCKO))
up <- list("Experimental\nInjury Repair Signature" = unique(allINJ),
           "PKD signature\n(up-regulated)" = unique(upPKD),
           "Up-regulated\ncKO genes" = unique(upCKO))
dw <- list("Experimental\nInjury Repair Signature" = unique(allINJ),
           "PKD signature\n(down-regulated)" = unique(dwPKD),
           "Down-regulated\ncKO genes" = unique(dwCKO))

## Make Venn Diagram
pkdsig.venn <-
    plot(euler(all),
         fill = c("grey","salmon","cyan"),
         alpha = rep(0.75, 3), lwd = 0,
         main = "All Genes (as in Malas Figure 3D)",
         quantities = TRUE)

#################################################################################
## Comparison to Malas data - Meta analysis table
#################################################################################
## Filter results for significant values (p<= 0.05) 
sigRes <- lapply(results, function(x) { x <- x[x$padj <= fdr,] })

## add data to meta data table
pkdMeta <- lapply(pkdMeta, function(x) {
    x$constam_male <- sigRes[[1]]$log2FoldChange[match(tolower(x$gene_symbol),
                                                       tolower(sigRes[[1]]$gene))]
    x$constam_female <- sigRes[[2]]$log2FoldChange[match(tolower(x$gene_symbol),
                                                         tolower(sigRes[[2]]$gene))]
    x$bicc1_cKOmf_up <- ifelse(!is.na(x$constam_male) & !is.na(x$constam_female),
                               "X",
                               "-")
    x$cKO_FC_avgMF <- ifelse(!is.na(x$constam_male) & !is.na(x$constam_female)
                            ,(x$constam_male + x$constam_female) / 2 ,
                             NA)
    x$baseMean_cKO <- ifelse(!is.na(x$constam_male) & !is.na(x$constam_female),
                             sigRes[[2]]$baseMean[match(tolower(x$gene_symbol),
                                                        tolower(sigRes[[2]]$gene))],
                             NA)
    x$padj_M <- ifelse(!is.na(x$constam_male) & !is.na(x$constam_female),
                       sigRes[[1]]$padj[match(tolower(x$gene_symbol),
                                              tolower(sigRes[[1]]$gene))],
                       NA)
    x$padj_F <- ifelse(!is.na(x$constam_male) & !is.na(x$constam_female),
                       sigRes[[2]]$padj[match(tolower(x$gene_symbol),
                                              tolower(sigRes[[2]]$gene))],
                       NA)
    x
})


injuryDEG <- lapply(injuryDEG, function(x) {
    x$cKO_M_log2FC <- sigRes[[1]]$log2FoldChange[match(tolower(x$associated_gene_name),
                                                       tolower(sigRes[[1]]$gene))]
    x$cKO_F_log2FC <- sigRes[[2]]$log2FoldChange[match(tolower(x$associated_gene_name),
                                                       tolower(sigRes[[2]]$gene))]
    x$bicc1_cKO_mf <- ifelse(!is.na(x$cKO_M_log2FC) & !is.na(x$cKO_F_log2FC),
                             "X",
                             "-")
    x$cKO_FC_avgMF <- ifelse(!is.na(x$cKO_M_log2FC) & !is.na(x$cKO_F_log2FC)
                            ,(x$cKO_M_log2FC + x$cKO_F_log2FC)/2,
                             NA)
    x$baseMean_cKO <- ifelse(!is.na(x$cKO_M_log2FC) & !is.na(x$cKO_F_log2FC),
                             sigRes[[2]]$baseMean[match(tolower(x$associated_gene_name),
                                                        tolower(sigRes[[2]]$gene))],
                             NA)
    x$cKO_M_padj <- ifelse(!is.na(x$cKO_M_log2FC) & !is.na(x$cKO_F_log2FC),
                           sigRes[[1]]$padj[match(tolower(x$associated_gene_name),
                                                  tolower(sigRes[[1]]$gene))],
                           NA)
    x$cKO_F_padj <- ifelse(!is.na(x$cKO_M_log2FC) & !is.na(x$cKO_F_log2FC),
                           sigRes[[2]]$padj[match(tolower(x$associated_gene_name),
                                                  tolower(sigRes[[2]]$gene))],
                           NA)
    x
})


## Genes that are significant in both male and female
genes <- lapply(results, function(x) {na.omit(x$gene_name[x$padj <= fdr])})
sigGenes <- intersect(as.character(genes[[1]]),as.character(genes[[2]]))
mgi <- gconvert(sigGenes, organism = "mmusculus", target = "MGI_ACC")
mgi <- mgi[!duplicated(mgi$input),]

dat <- tibble(data.frame(
    ##Gene Information
    gene        = sigGenes,
    mgi_symbol  = mgi$target[match(sigGenes,mgi$input)],
    ensembl_id  = as.character(na.omit(results[[1]]$ensembl_id[match(sigGenes,results[[1]]$gene_name)])),
    entrez_id   = as.character(results[[1]]$entrez[match(sigGenes, results[[1]]$gene_name)]),
    description = as.character(results[[1]]$description[match(sigGenes,results[[1]]$gene_name)]),
    baseMean = na.omit(results[[1]]$baseMean[match(sigGenes,results[[1]]$gene_name)]),
    ## Male information
    FC_M     = 2^(na.omit(results[[1]]$log2FoldChange[match(sigGenes,results[[1]]$gene_name)])),
    log2FC_M = na.omit(results[[1]]$log2FoldChange[match(sigGenes,results[[1]]$gene_name)]),
    lfcSE_M  = na.omit(results[[1]]$lfcSE[match(sigGenes, results[[1]]$gene_name)]),
    stat_M   = na.omit(results[[1]]$stat[match(sigGenes, results[[1]]$gene_name)]),
    pval_M   = na.omit(results[[1]]$pvalue[match(sigGenes,results[[1]]$gene_name)]),
    padj_M   = na.omit(results[[1]]$padj[match(sigGenes,results[[1]]$gene_name)]),
    ## Female information
    FC_F     = 2^(na.omit(results[[2]]$log2FoldChange[match(sigGenes,results[[2]]$gene_name)])),
    log2FC_F = na.omit(results[[2]]$log2FoldChange[match(sigGenes, results[[2]]$gene_name)]),
    lfcSE_F  = na.omit(results[[2]]$lfcSE[match(sigGenes, results[[2]]$gene_name)]),
    stat_F   = na.omit(results[[2]]$stat[match(sigGenes, results[[2]]$gene_name)]),
    pval_F   = na.omit(results[[2]]$pvalue[match(sigGenes, results[[2]]$gene_name)]),
    padj_F   = na.omit(results[[2]]$padj[match(sigGenes, results[[2]]$gene_name)])
))

## Add Merged Data to data frame
dat <- dat |>
    mutate(log2FC_avgMF = (log2FC_M + log2FC_F) / 2,.before = FC_M) |>
    mutate(cko_change = ifelse(log2FC_avgMF > 0,"up", "dw"))

## Add PKD meta data
pkd.up <- cbind(pkdMeta[[1]][1:4],
                apply(pkdMeta[[1]][5:11], 2, function(x) {gsub("X","up",x)}))
pkd.dw <- cbind(pkdMeta[[2]][1:4],
                apply(pkdMeta[[2]][5:11], 2, function(x) {gsub("X","dw",x)}))
pkd <- rbind(pkd.up,pkd.dw)
pkd <- pkd[match(tolower(dat$gene),tolower(pkd$gene_symbol)),] ## Reorder table

dat <- bind_cols(dat,pkd[5:11]) ## Add pkd data to data table

## Add injury data
inj <- lapply(injuryDEG, function(x) {
    x <- x$associated_gene_name[match(tolower(dat$gene),tolower(x$associated_gene_name))]
    x[!is.na(x)] <- "X"
    x
})
inj <- do.call(cbind,inj)
colnames(inj) <- c("inj. early", "inj. late")

res <- bind_cols(dat,inj) ## Add injury data to data table

################################################################################
## Comparison to Malas data - Dotplot
################################################################################
## Order DEGs by malas list
res <- res[res$ensembl_id %in% malasDEG$ensembl_gene_id,] ## filter res for in malas data
## order and filter malas data to res
malasDEG <- malasDEG[na.omit(match(res$ensembl_id, malasDEG$ensembl_gene_id)),]

## Create merged table
mergeData <- data.frame(ensembl_id = res$ensembl_id,
                        gene_name_malas = malasDEG$associated_gene_name,
                        gene_name_constam = res$gene,
                        description_malas = malasDEG$description,
                        description_constam = res$description,
                        logFC_malas = malasDEG$log_fc,
                        fdr_malas = malasDEG$fdr,
                        logFC_MFavg_constam = (res$log2FC_F + res$log2FC_M)/2,
                        padj_male_constam = res$padj_M,
                        padj_female_constam = res$padj_F
                        )
mergeData <- na.omit(mergeData)

## Define significance - for point colors
mergeData$sig_male <- "no"
mergeData$sig_male[mergeData$padj_male_constam <= fdr] <-  "yes"
mergeData$sig_female <- "no"
mergeData$sig_female[mergeData$padj_female_constam <= fdr] <-  "yes"

## Only plot points that are significant in both males and females
mergeData <- mergeData[mergeData$sig_male == "yes" & mergeData$sig_female == "yes",]

## Define genes to label - same diredtion (Gene lists from October 14, 2025)
dw.goi <- c("Miox","Slc23a1","Fam151a","Slc4a1","Osbpl6","Egf","Slc12a3","Aqp3",
            "Klk1","Slc43a2","Pde8a","Slc2a9","Mettl8","Gcgr","Atp4a")
up.goi <- c("Col8a1","Nnmt","Anxa1","Gstt1","Ankrd1","Lcn2","Ccl2","Tacstd2",
            "Gdf15","Cd14","Socs3","Havcr1","Lgals3","Thbs1","Glipr2","Klf6",
            "Klf5","Pcsk5","Pkd2")

## Malas opposite genes
mal.opp <- c("Aif1","Tpd52l1","Stmn1","Snrpd1","Slc35g1","Cldn8","Nudt4",
             "Plekho2","Pcsk5","Slit3","Slc23a2","Pkd2","Abcc3","Plbd2",
             "Sec14l1","Col4a2","Taok3","Nrp1","Slc6a18","Mep1a","Bcl9l",
             "Cdhr2")

## labels that are sig in both
mergeData$sig_label <- 
    mergeData$gene_name_malas %in% c(dw.goi,up.goi) &
    mergeData$sig_male == "yes" &
    mergeData$sig_female == "yes"

mergeData$sig_opp_label <- 
    mergeData$gene_name_malas %in% c(mal.opp) & ## 12 genes
    mergeData$sig_male == "yes" &
    mergeData$sig_female == "yes"

## Create plots - Use male statistics (logFC, etc)
pkd_bicc_overlap <-
    ggplot(mergeData, aes(x=logFC_malas, y=logFC_MFavg_constam,
                          label = gene_name_constam)) +
    geom_point(color = "#E69F00") + 
    geom_text( 
      data=mergeData %>% filter(sig_label == TRUE),
      aes(label=gene_name_constam),
      nudge_y = 0.04,
      check_overlap = F,
      color = "black",
      size = 3) +
    geom_text( 
      data=mergeData %>% filter(sig_opp_label == TRUE),
      aes(label=gene_name_constam),
      nudge_y = 0.04, 
      check_overlap = F,
      color = "red",
      size = 3) +
    scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
    scale_x_continuous(minor_breaks = seq(from = round(min(mergeData$logFC_malas),0),
                                          to =round(max(mergeData$logFC_malas),0),
                                          by = 0.5)) + 
    scale_y_continuous(minor_breaks = seq(from = round(min(mergeData$logFC_MFavg_constam),0),
                                          to =round(max(mergeData$logFC_MFavg_constam),0),
                                          by = 0.5)) + 
    labs(title = "Meta analysis with Malas (full DEG) - MF average") +
    theme_bw()  


#################################################################################
## Print Plots
#################################################################################
## Fig S7B - perc mapped plot
alnStats.perc

## Fig S7C - raw mapped plot
alnStats.raw

## Fig 5A - PCA plot
pca.plot

## Fig 5B - sex overlap venns 
sexOverlap.venn 

## Fig 5C - Volcano Plot (male and female)
volcano.plots

## Fig 5D - PKD signature and Injury Repair venn diagram
pkdsig.venn

## Fig 5E - Malas PKD signature (DEGs) and Bicc1 KO dotplot (male only)
pkd_bicc_overlap

#################################################################################
## Print Tables
#################################################################################
## Table S3 - DEGs
writexl::write_xlsx(
             lapply(results, function(x) { x[x$padj <= fdr,] }),
             "./data/2_DifferentialExpression/DESeq2_DEGs.xlsx"
         )

## Table S4 - Gene Ontologies
writexl::write_xlsx(
             list(
                 upGenes_both = up.in.both,
                 dwGenes_both = dw.in.both,
                 upGO = up.gost,
                 dwGO = dw.gost
             ),
             "./data/2_DifferentialExpression/GO_terms.xlsx"
         )

## Table S5 - Meta analysis overlap
writexl::write_xlsx(
             res,
             "./data/3_PKDsignature/PKDsig_BICC1ko_overlaps.xlsx"
         )

