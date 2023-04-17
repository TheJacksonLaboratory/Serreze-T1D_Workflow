---
title: "QTL Analysis - binary [Hetero: ICI Myocarditis Yes vs ICI Myocarditis No]"
author: "Belinda Cornes"
date: "`r Sys.Date()`"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = FALSE,warning=FALSE)

library("kableExtra")
library("knitr")
options(stringsAsFactors = FALSE)
library(data.table) 
library(tidyr)     
library(mclust)     
library(rhdf5)      
library(optparse)
library(dplyr)
library(cluster)
library(readxl)
library(psych)
library(tibble)
library(ggplot2)
library(reshape2)
library(qtl2)
#library(qtl)
library(abind)

setwd("/projects/serreze-lab/USERS/corneb/qtl2/workflowr/2023.04.13/Serreze-T1D_Workflow")

```

## Data Information

+ Myocarditis Status: Yes, No
+ Murine MHC KO Status: Het
+ Drug Treatment(s): ICI
+ Clinical Phenotype(s): Sick, EOI
+ Covariates (Myocarditis Score, sex)

## Loading Data

We will load the data and subset indivials out that are in the groups of interest.

```{r data, echo=TRUE}

load("data/gm_allqc_7.batches_myo_mis.RData")

#gm_allqc
gm=gm_allqc
gm

#pr <- readRDS("data/serreze_probs_allqc.rds")
#pr <- readRDS("data/serreze_probs.rds")

##extracting animals with ici and pbs group status
#miceinfo <- gm$covar[gm$covar$group == "PBS" | gm$covar$group == "ICI",]
#table(miceinfo$group)
#mice.ids <- rownames(miceinfo)

#gm <- gm[mice.ids]
#gm
#table(gm$covar$group)

#gm$covar$het.ici.myo.yes_vs_het.ici.myo.no <- ifelse(gm$covar$group == "PBS", 0, 1)
#gm.full <- gm

covars <- read_csv("data/covar_corrected_het-ici-myo-yes.vs.het-ici-myo-no_7.batches_myo_mis.csv")
#removing any missing info
#covars <- subset(covars, covars$het.ici.myo.yes_vs_het.ici.myo.no!='')
nrow(covars)
table(covars$"Myocarditis Status")
table(covars$"Murine MHC KO Status")
table(covars$"Drug Treatment")
table(covars$"clinical pheno")

#keeping only informative mice
gm <- gm[covars$Mouse.ID]
gm
table(gm$covar$"Myocarditis Status")
table(gm$covar$"Murine MHC KO Status")
table(gm$covar$"Drug Treatment")
table(gm$covar$"clinical pheno")

#pr.qc.ids <- pr
#for (i in 1:20){pr.qc.ids[[i]] = pr.qc.ids[[i]][covars$Mouse.ID,,]}

##removing problmetic marker

#gm <- drop_markers(gm, "UNCHS013106")

##dropping monomorphic markers within the dataset

g <- do.call("cbind", gm$geno)

gf_mar <- t(apply(g, 2, function(a) table(factor(a, 1:2))/sum(a != 0)))
#gn_mar <- t(apply(g, 2, function(a) table(factor(a, 1:2))))

gf_mar <- gf_mar[gf_mar[,2] != "NaN",]

count <- rowSums(gf_mar <=0.05)
low_freq_df <- merge(as.data.frame(gf_mar),as.data.frame(count), by="row.names",all=T)
low_freq_df[is.na(low_freq_df)] <- ''
low_freq_df <- low_freq_df[low_freq_df$count == 1,]
rownames(low_freq_df) <- low_freq_df$Row.names

low_freq <- find_markerpos(gm, rownames(low_freq_df))
low_freq$id <- rownames(low_freq)

nrow(low_freq)

low_freq_bad <- merge(low_freq,low_freq_df, by="row.names",all=T)
names(low_freq_bad)[1] <- c("marker")

gf_mar <- gf_mar[gf_mar[,2] != "NaN",]
MAF <- apply(gf_mar, 1, function(x) min(x))
MAF <- as.data.frame(MAF)
MAF$index <- 1:nrow(gf_mar)
gf_mar_maf <- merge(gf_mar,as.data.frame(MAF), by="row.names")
gf_mar_maf <- gf_mar_maf[order(gf_mar_maf$index),]

gfmar <- NULL
gfmar$gfmar_mar_0 <- sum(gf_mar_maf$MAF==0)
gfmar$gfmar_mar_1 <- sum(gf_mar_maf$MAF< 0.01)
gfmar$gfmar_mar_5 <- sum(gf_mar_maf$MAF< 0.05)
gfmar$gfmar_mar_10 <- sum(gf_mar_maf$MAF< 0.10)
gfmar$gfmar_mar_15 <- sum(gf_mar_maf$MAF< 0.15)
gfmar$gfmar_mar_25 <- sum(gf_mar_maf$MAF< 0.25)
gfmar$gfmar_mar_50 <- sum(gf_mar_maf$MAF< 0.50)
gfmar$total_snps <- nrow(as.data.frame(gf_mar_maf))

gfmar <- t(as.data.frame(gfmar))
gfmar <- as.data.frame(gfmar)
gfmar$count <- gfmar$V1

gfmar[c(2)] %>%
  kable(escape = F,align = c("ccccccccc"),linesep ="\\hline") %>%
  kable_styling(full_width = F) %>%
  kable_styling("striped", full_width = F)  %>%
  row_spec(8 ,bold=T,color= "white",background = "black")


gm_qc <- drop_markers(gm, low_freq_bad$marker)
gm_qc <- drop_nullmarkers(gm_qc)

gm = gm_qc
gm 

### dropping disproportionate markers
#dismark <- read.csv("data/het-ici-myo-yes.vs.het-ici-myo-no_marker.freq_low.geno.freq.removed_geno.ratio_7.batches_myo_mis.csv")
#nrow(dismark)
#names(dismark)[1] <- c("marker")
#ismark <- dismark[!dismark$Include,]
#nrow(dismark)

#gm_qc_dis <- drop_markers(gm_qc, dismark$marker)
#gm_qc_dis <- drop_nullmarkers(gm_qc_dis)

#gm = gm_qc_dis
#gm

markers <- marker_names(gm)
gmapdf <- read.csv("data/genetic_map_7.batches_myo.csv")
pmapdf <- read.csv("data/physical_map_7.batches_myo.csv")
#mapdf <- merge(gmapdf,pmapdf, by=c("marker","chr"), all=T)
#rownames(mapdf) <- mapdf$marker
#mapdf <- mapdf[markers,]
#names(mapdf) <- c('marker','chr','gmapdf','pmapdf')
#mapdfnd <- mapdf[!duplicated(mapdf[c(2:3)]),]

pr.qc <- calc_genoprob(gm)

colnames(covars) <- gsub(" ", ".", colnames(covars))

```


## Genome-wide scan
For each of the phenotype analyzed, permutations were used for each model to obtain genome-wide LOD significance threshold for p < 0.01, p < 0.05, p < 0.10,  respectively, separately for X and automsomes (A).  

The table shows the estimated significance thresholds from permutation test. 

We also looked at the kinship to see how correlated each sample is. Kinship values between pairs of samples range between 0 (no relationship) and 1.0 (completely identical). The darker the colour the more indentical the pairs are.

```{r permutation, fig.height = 6, fig.width = 9.5, fig.align = "center", echo=TRUE}

#Xcovar <- get_x_covar(gm)
addcovar = model.matrix(~Histology.Score, data = covars)[,-1]
covars$het.ici.myo.yes_vs_het.ici.myo.no= as.numeric(covars$het.ici.myo.yes_vs_het.ici.myo.no)

kinship <- calc_kinship(pr.qc)
heatmap(kinship)

operm <- scan1perm(pr.qc, covars["het.ici.myo.yes_vs_het.ici.myo.no"], model="binary", addcovar=addcovar, n_perm=1000, perm_Xsp=TRUE, chr_lengths=chr_lengths(gm$gmap))

summary_table<-data.frame(unclass(summary(operm, alpha=c(0.01,  0.05, 0.1))))
names(summary_table) <- c("autosomes","X")
summary_table$significance.level <- rownames(summary_table)

rownames(summary_table) <- NULL

summary_table[c(3,1:2)] %>%
  kable(escape = F,align = c("ccc")) %>%
  kable_styling("striped", full_width = T) %>%
  column_spec(1, bold=TRUE)

```
The figures below show QTL maps for each phenotype

```{r Genome-wide scan, fig.height = 6, fig.width = 9.5, fig.align = "center", echo=TRUE}

#out <- scan1(pr.qc, covars["het.ici.myo.yes_vs_het.ici.myo.no"], Xcovar=Xcovar, model="binary")
out <- scan1(pr.qc, covars["het.ici.myo.yes_vs_het.ici.myo.no"], model="binary",addcovar=addcovar)

summary_table<-data.frame(unclass(summary(operm, alpha=c(0.01,  0.05, 0.1))))


plot_lod<-function(out,map){
  for (i in 1:dim(out)[2]){
    #png(filename=paste0("/Users/chenm/Documents/qtl/Jai/",colnames(out)[i],  "_lod.png"))
    
    ymx <- maxlod(out) # overall maximum LOD score
    plot(out, map, lodcolumn=i, col="slateblue", ylim=c(0, ymx+0.5))
    #legend("topright", lwd=2, colnames(out)[i], bg="gray90")
    title(main = paste0(colnames(out)[i], " [positions in cM]"))
    add_threshold(map,  summary(operm,alpha=0.1), col = 'purple')
    add_threshold(map,  summary(operm, alpha=0.05), col = 'red')
    add_threshold(map,  summary(operm, alpha=0.01), col = 'blue')

    ##par(mar=c(5.1, 6.1, 1.1, 1.1))
    #ymx <- 11 # overall maximum LOD score
    #plot(out, map, lodcolumn=i, col="slateblue", ylim=c(0, ymx+0.5))
    ##legend("topright", lwd=2, colnames(out)[i], bg="gray90")
    #title(main = paste0(colnames(out)[i], " [positions in cM] \n(using same scale as eoi vs ici for easier comparison)"))
    #add_threshold(map,  summary(operm, alpha=0.1), col = 'purple')
    #add_threshold(map,  summary(operm, alpha=0.05), col = 'red')
    #add_threshold(map,  summary(operm, alpha=0.01), col = 'blue')
    ##for (j in 1: dim(summary_table)[1]){
    ##  abline(h=summary_table[j, i],col="red")
    ##  text(x=400, y =summary_table[j, i]+0.12, labels = paste("p=", row.names(summary_table)[j]))
    ##}
    ##dev.off()
  }
}

plot_lod(out,gm$gmap)

```

## LOD peaks

The table below shows QTL peaks associated with the phenotype. We use the 95% threshold from the permutations to find peaks. 

### Centimorgan (cM)

```{r LOD peaks, fig.height = 6, fig.width = 9.5, fig.align = "center", echo=TRUE, results="asis"}

peaks <- find_peaks(out, gm$gmap, threshold=summary(operm,alpha=0.05)$A, thresholdX = summary(operm,alpha=0.05)$X, peakdrop=3, drop=1.5)

if(nrow(peaks) >0){
peaks$marker <- find_marker(gm$gmap, chr=peaks$chr,pos=peaks$pos)
names(peaks)[2] <- c("phenotype")
peaks <- peaks[-1]

rownames(peaks) <- NULL
print(kable(peaks, escape = F, align = c("cccccccc"), "html") 
  %>% kable_styling("striped", full_width = T)%>%
  column_spec(1, bold=TRUE)
  )

#plot only peak chromosomes

plot_lod_chr<-function(out,map,chrom){
  for (i in 1:dim(out)[2]){
    #png(filename=paste0("/Users/chenm/Documents/qtl/Jai/",colnames(out)[i],  "_lod.png"))
    
    #par(mar=c(5.1, 6.1, 1.1, 1.1))
    ymx <- maxlod(out) # overall maximum LOD score
    plot(out, map, chr = chrom, lodcolumn=i, col="slateblue", ylim=c(0, ymx+0.5))
    #legend("topright", lwd=2, colnames(out)[i], bg="gray90")
    title(main = paste0(colnames(out)[i], " - chr", chrom, " [positions in cM]"))
    add_threshold(map,  summary(operm,alpha=0.1), col = 'purple')
    add_threshold(map,  summary(operm, alpha=0.05), col = 'red')
    add_threshold(map,  summary(operm, alpha=0.01), col = 'blue')
    #for (j in 1: dim(summary_table)[1]){
    #  abline(h=summary_table[j, i],col="red")
    #  text(x=400, y =summary_table[j, i]+0.12, labels = paste("p=", row.names(summary_table)[j]))
    #}
    #dev.off()

    
    #ymx <- 11
    #plot(out, map, chr = chrom, lodcolumn=i, col="slateblue", ylim=c(0, ymx+0.5))
    ##legend("topright", lwd=2, colnames(out)[i], bg="gray90")
    #title(main = paste0(colnames(out)[i], " - chr", chrom, " [positions in cM]\n(using same scale as eoi vs. ici for easier comparison)"))
    #add_threshold(map,  summary(operm,alpha=0.1), col = 'purple')
    #add_threshold(map,  summary(operm, alpha=0.05), col = 'red')
    #add_threshold(map,  summary(operm, alpha=0.01), col = 'blue')

  }
}


for(i in unique(peaks$chr)){
#for (i in 1:nrow(peaks)){
  #plot_lod_chr(out,gm$gmap, peaks$chr[i])
  plot_lod_chr(out,gm$gmap, i)
}

} else {
  print(paste0("There are no peaks that have a LOD that reaches suggestive (p<0.05) level of ",summary(operm,alpha=0.05)$A, " [autosomes]/",summary(operm,alpha=0.05)$X, " [x-chromosome]"))
}
```

### Megabase (MB)

```{r LOD peaks_mba, fig.height = 6, fig.width = 9.5, fig.align = "center", echo=TRUE, results="asis"}

print("peaks in MB positions")

peaks_mba <- find_peaks(out, gm$pmap, threshold=summary(operm,alpha=0.05)$A, thresholdX = summary(operm,alpha=0.05)$X, peakdrop=3, drop=1.5)

if(nrow(peaks) >0){
peaks_mba$marker <- find_marker(gm$pmap, chr=peaks_mba$chr,pos=peaks_mba$pos)
names(peaks_mba)[2] <- c("phenotype")
peaks_mba <- peaks_mba[-1]


rownames(peaks_mba) <- NULL
print(kable(peaks_mba, escape = F, align = c("cccccccc"), "html") 
  %>% kable_styling("striped", full_width = T)%>%
  column_spec(1, bold=TRUE)
  )

plot_lod_chr_mb<-function(out,map,chrom){
  for (i in 1:dim(out)[2]){
    #png(filename=paste0("/Users/chenm/Documents/qtl/Jai/",colnames(out)[i],  "_lod.png"))
    
    #par(mar=c(5.1, 6.1, 1.1, 1.1))
    ymx <- maxlod(out) # overall maximum LOD score
    plot(out, map, chr = chrom, lodcolumn=i, col="slateblue", ylim=c(0, ymx+0.5))
    #legend("topright", lwd=2, colnames(out)[i], bg="gray90")
    title(main = paste0(colnames(out)[i], " - chr", chrom, " [positions in MB]"))
    add_threshold(map,  summary(operm,alpha=0.1), col = 'purple')
    add_threshold(map,  summary(operm, alpha=0.05), col = 'red')
    add_threshold(map,  summary(operm, alpha=0.01), col = 'blue')
    #for (j in 1: dim(summary_table)[1]){
    #  abline(h=summary_table[j, i],col="red")
    #  text(x=400, y =summary_table[j, i]+0.12, labels = paste("p=", row.names(summary_table)[j]))
    #}
    #dev.off()

    #ymx <- 11
    #plot(out, map, chr = chrom, lodcolumn=i, col="slateblue", ylim=c(0, ymx+0.5))
    ##legend("topright", lwd=2, colnames(out)[i], bg="gray90")
    #title(main = paste0(colnames(out)[i], " - chr", chrom, " [positions in MB]\n(using same scale as eoi vs. ici for easier comparison)"))
    #add_threshold(map,  summary(operm,alpha=0.1), col = 'purple')
    #add_threshold(map,  summary(operm, alpha=0.05), col = 'red')
    #add_threshold(map,  summary(operm, alpha=0.01), col = 'blue')


  }
}

for(i in unique(peaks_mba$chr)){
#for (i in 1:nrow(peaks_mba)){
  #plot_lod_chr_mb(out,gm$pmap, peaks_mba$chr[i])
  plot_lod_chr_mb(out,gm$pmap,i)
}

} else {
  print(paste0("There are no peaks that have a LOD that reaches suggestive (p<0.05) level of ",summary(operm,alpha=0.05)$A, " [autosomes]/",summary(operm,alpha=0.05)$X, " [x-chromosome]"))
}

```

## QTL effects

For each peak LOD location we give a list of gene

```{r QTL effects, fig.height = 6, fig.width = 9.5, fig.align = "center", echo=TRUE, results="asis"}

query_variants <- create_variant_query_func("code/cc_variants.sqlite")
query_genes <- create_gene_query_func("code/mouse_genes_mgi.sqlite")

if(nrow(peaks) >0){
for (i in 1:nrow(peaks)){


  g <- maxmarg(pr.qc, gm$gmap, chr=peaks$chr[i], pos=peaks$pos[i], return_char=TRUE)
  #png(filename=paste0("/Users/chenm/Documents/qtl/Jai/","qtl_effect_", i, ".png"))
  #par(mar=c(4.1, 4.1, 1.5, 0.6))
  plot_pxg(g, covars[,peaks$phenotype[i]], ylab=peaks$phenotype[i], sort=FALSE)
  title(main = paste0("chr: ", chr=peaks$chr[i],"; pos: ", peaks$pos[i], "cM /",peaks_mba$pos[i],"MB\n(",peaks$phenotype[i]," )"), line=0.2)
  ##dev.off()

  chr = peaks$chr[i]

# Plot 2
  pr_sub <- pull_genoprobint(pr.qc, gm$gmap, chr, c(peaks$ci_lo[i], peaks$ci_hi[i]))
  blup <- scan1blup(pr.qc[,chr], covars[peaks$phenotype[i]],addcovar = addcovar)
  blup_sub <- scan1blup(pr_sub[,chr], covars[peaks$phenotype[i]], addcovar = addcovar)

  write.csv(as.data.frame(blup_sub), paste0("data/het-ici-myo-yes.vs.het-ici-myo-no_blup_sub_chr-",chr,"_peak.marker-",peaks$marker[i],"_lod.drop-1.5_snpsqc_7.batches_myo_mis.csv"), quote=F)

  plot_coef(blup, 
       gm$gmap, columns=1:2,
       bgcolor="gray95", legend="bottomleft", 
       main = paste0("chr: ", chr=peaks$chr[i],"; pos: ", peaks$pos[i], "cM /",peaks_mba$pos[i],"MB\n(",peaks$phenotype[i]," [scan1blup; positions in cM])")
       )

  plot_coef(blup_sub, 
       gm$gmap, columns=1:2,
       bgcolor="gray95", legend="bottomleft", 
       main = paste0("chr: ", chr=peaks$chr[i],"; pos: ", peaks$pos[i], "cM /",peaks_mba$pos[i],"MB\n(",peaks$phenotype[i],"; 1.5 LOD drop interval [scan1blup; positions in cM])")
       )

  #Table 1
  chr = peaks_mba$chr[i]
  start=as.numeric(peaks_mba$ci_lo[i])
  end=as.numeric(peaks_mba$ci_hi[i])

  genesgss = query_genes(chr, start, end)

  write.csv(genesgss, file=paste0("data/het-ici-myo-yes.vs.het-ici-myo-no_genes_chr-",chr,"_peak.marker-",peaks$marker[i],"_lod.drop-1.5_snpsqc_7.batches_myo_mis.csv"), quote=F)

  rownames(genesgss) <- NULL
  genesgss$strand_old = genesgss$strand
  genesgss$strand[genesgss$strand=="+"] <- "positive"
  genesgss$strand[genesgss$strand=="-"] <- "negative"

  print(kable(genesgss[,c("chr","type","start","stop","strand","ID","Name","Dbxref","gene_id","mgi_type","description")], "html") %>% kable_styling("striped", full_width = T))


}

} else {
  print(paste0("There are no peaks that have a LOD that reaches suggestive (p<0.05) level of ",summary(operm,alpha=0.05)$A, " [autosomes]/",summary(operm,alpha=0.05)$X, " [x-chromosome]"))
}

```

```{r QTL effects blup all chrs, echo=TRUE, include=FALSE}

#for (i in names(pr.qc)){
#  chr = i
#  blup <- scan1blup(pr.qc[,chr], covars["het.ici.myo.yes_vs_het.ici.myo.no"])
#  write.csv(as.data.frame(blup), paste0("data/het-ici-myo-yes.vs.het-ici-myo-no_blup.qc_chr-",chr,"_snpsqc_7.batches_myo_mis.csv"), quote=F)
#}

########################


#gm.full

#pr.full <- calc_genoprob(gm.full)

#for (i in names(pr.full)){
#  chr = i
#  blup <- scan1blup(pr.full[,chr], gm.full$covar["het.ici.myo.yes_vs_het.ici.myo.no"])
#  write.csv(as.data.frame(blup), paste0("data/het-ici-myo-yes.vs.het-ici-myo-no_blup.full_chr-",chr,"_snpsqc_7.batches_myo_mis.csv"), quote=F)
#}

#save(out, operm, gm.full, gm, pr.qc, pr.full, file="data/het-ici-myo-yes.vs.het-ici-myo-no_scanone_snpsqc_7.batches_myo_mis.Rdata")

```

## R/qtl

```{r converting files, echo=TRUE, include=FALSE}

#converting qtl2 data into a qtl data object
#genotypes
#genos <- as.data.frame(do.call("cbind", gm$geno))
genos <- read.csv("data/sample_geno_bc_7.batches_myo.csv")
rownames(genos) <- genos$marker
genos <- genos[markers,]
colnames(genos) <- gsub("\\.","-", colnames(genos))
colnames(genos) <- gsub("^X","",colnames(genos))
genos <- genos[,qtl2::ind_ids(gm)]
genos[genos=="AA"] <- "A"
genos[genos=="AB"] <- "H"
#genos[genos=="-"] <- NA
genos[genos=="-"] <- "-"

#phenotypes
phenos <- covars
phenos$id = phenos$Mouse.ID


#genetic map
gen.map <- as.data.frame(do.call("abind", gm$gmap))
names(gen.map) <- c("pos")
gen.map$marker = rownames(gen.map)
## getting chromsomes
gmapdf$index <- 1:nrow(gmapdf)
gen.map.chr <- merge(gmapdf,gen.map, by=c("marker","pos"), all.y=T)
gen.map.chr <- gen.map.chr[order(gen.map.chr$index),]
#gen.map.chr <- gen.map.chr[,-gen.map.chr$index]
gen.map.chr <- gen.map.chr[c("marker","chr","pos")]
rownames(gen.map.chr) <- gen.map.chr$marker


#combining all
qtl.d1 <- as.data.frame(rbind(t(gen.map.chr[,-1]),t(genos)))
qtl.d1$ID <- rownames(qtl.d1) 
qtl.d1$index <- 1:nrow(qtl.d1)
qtl.d2 <- merge(qtl.d1, phenos[,c("id","het.ici.myo.yes_vs_het.ici.myo.no","sex","Histology.Score")], by.x=c("ID"), by.y=c("id"), all=T)
qtl.d2 <- qtl.d2[order(qtl.d2$index),]
qtl.d3 <- qtl.d2[c(1,(ncol(qtl.d2)-2),(ncol(qtl.d2)-1),ncol(qtl.d2),2:(ncol(qtl.d2)-4))]
#qtl.d3[!is.na(qtl.d3$sex),]$sex <- 0

#putting missing in csv files
qtl.d3[1:2,1:4] <- ""

qtl.d3[qtl.d3$sex=="F",]$sex = 0
qtl.d3[qtl.d3$sex=="M",]$sex = 1
write.csv(qtl.d3,"data/het-ici-myo-yes.vs.het-ici-myo-no_gm_qtl_snpsqc_7.batches_myo_mis.csv", quote=F, row.names=F)

```

### scanone

```{r scanone, fig.height = 6, fig.width = 9.5, fig.align = "center", echo=FALSE, eval=FALSE, include=FALSE}



gm

#detach("package:qtl2", unload=TRUE)
#library(qtl)

cross <- qtl::read.cross("csv", file = "data/het-ici-myo-yes.vs.het-ici-myo-no_gm_qtl_snpsqc_7.batches_myo_mis.csv",alleles=c("A","B"))
cross <- qtl::jittermap(cross)

summary(cross)

cross.probs <- qtl::calc.genoprob(cross)
addcovar <- pull.pheno(cross, c("Histology.Score"))

print("method == hk")

scanone.hk <-qtl::scanone(cross.probs, pheno.col="het.ici.myo.yes_vs_het.ici.myo.no" , addcovar = addcovar, model="binary", method="hk")
operm.hk <- qtl::scanone(cross.probs, method = "hk", pheno.col="het.ici.myo.yes_vs_het.ici.myo.no", addcovar = addcovar, n.perm = 1000, perm.Xsp = TRUE, model="binary", verbose=FALSE)
#plot(operm.hk)
print(summary(operm.hk, alpha=c(0.01,  0.05, 0.1)))

#plot(scanone.hk, bandcol = "grey90",lty=1, cex=1, col = "steelblue")  
#qtl::add.threshold(scanone.hk,  perms= operm.hk, alpha=0.01, col = 'blue')
#qtl::add.threshold(scanone.hk,  perms= operm.hk, alpha=0.05, col = 'red')
#qtl::add.threshold(scanone.hk,  perms= operm.hk, alpha=0.1, col = 'purple')

ymx <- maxlod(out) # overall maximum LOD score
plot(scanone.hk, bandcol = "grey90",lty=1, cex=1, col = "slateblue", ylim=c(0, ymx+0.5))
title(main = paste0(colnames(out), " [positions in cM]"))  
qtl::add.threshold(scanone.hk,  perms= operm.hk, alpha=0.01, col = 'blue')
qtl::add.threshold(scanone.hk,  perms= operm.hk, alpha=0.05, col = 'red')
qtl::add.threshold(scanone.hk,  perms= operm.hk, alpha=0.1, col = 'purple')

ymx <- 11
plot(scanone.hk, bandcol = "grey90",lty=1, cex=1, col = "slateblue", ylim=c(0, ymx+0.5))
title(main = paste0(colnames(out), " [positions in cM]\n(using same scale as eoi vs ici for easier comparison)"))
qtl::add.threshold(scanone.hk,  perms= operm.hk, alpha=0.01, col = 'blue')
qtl::add.threshold(scanone.hk,  perms= operm.hk, alpha=0.05, col = 'red')
qtl::add.threshold(scanone.hk,  perms= operm.hk, alpha=0.1, col = 'purple')

print(as.data.frame(summary(scanone.hk, perms=operm.hk, pvalues=TRUE, format="allpeaks")))

print("all peaks with a p-value less or equal to 0.05 (suggestive)")
print(as.data.frame(summary(scanone.hk, perms=operm.hk, alpha=0.05, pvalues=TRUE, format="allpeaks")))

#print("method == ehk")

#scanone.ehk <-qtl::scanone(cross.probs, pheno.col="het.ici.myo.yes_vs_het.ici.myo.no" , model="binary", method="ehk")
#operm.ehk <- qtl::scanone(cross.probs, method = "ehk", pheno.col="het.ici.myo.yes_vs_het.ici.myo.no", n.perm = 1000, perm.Xsp = TRUE, model="binary", verbose=FALSE)
#plot(operm.ehk)
#print(summary(operm.ehk, alpha=c(0.01,  0.05, 0.1)))

#plot(scanone.ehk, bandcol = "grey90",lty=1, cex=1, col = "steelblue")  
#qtl::add.threshold(scanone.ehk,  perms= operm.ehk, alpha=0.01, col = 'blue')
#qtl::add.threshold(scanone.ehk,  perms= operm.ehk, alpha=0.05, col = 'red')
#qtl::add.threshold(scanone.ehk,  perms= operm.ehk, alpha=0.1, col = 'purple')

#print(as.data.frame(summary(scanone.ehk)))
#print(as.data.frame(summary(scanone.ehk, perms=operm.ehk, alpha=0.05, pvalues=TRUE, format="allpeaks")))

```
