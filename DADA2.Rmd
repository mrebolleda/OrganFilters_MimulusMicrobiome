---
title: "Dada_diversity_bigdata"
author: "Maria Rebolleda-Gomez"
date: "19 de marzo de 2019"
output: html_document
---

```{r}
library("vegan")
library("phyloseq")
library("ggplot2")
library("tidyr")
library(plyr)
library(VennDiagram)
library(picante)
library(ape)
library(DECIPHER)
library(phangorn)
library(dada2)
```


### ASV table and analyses
Check quality and trim:
```{r}
path <- "/media/mariargz/Extra Drive 1/MimulusMicrobiome/raw/Demultiplexed_Rebolleda_Gomez_fastqs_R2/"
fnFs <- sort(list.files(path, pattern="L001_R1", full.names = T))
fnRs <- sort(list.files(path, pattern="L001_R2", full.names = T))

#sample.names=substr(basename(fnFs),12,14)
sample.names=sapply(strsplit(basename(fnFs),"_"),'[',1)

plotQualityProfile(fnFs[(sample(1:length(fnFs),3))])
plotQualityProfile(fnRs[(sample(1:length(fnRs),3))])

filtFs=file.path(path,"filtered",paste0(sample.names,"_F_filt.fastq.gz"))
filtRs=file.path(path,"filtered",paste0(sample.names,"_R_filt.fastq.gz"))


out1=filterAndTrim(fnFs,filtFs,fnRs,filtRs, truncLen=c(150,140),
              maxN=0, maxEE=c(1,3), truncQ=2, rm.phix=TRUE,
              compress=TRUE) 
#For read 2 there were many files getting lost and when I checked back those files took up very small memory so I removed them from this analysis (SEE TABLE X). 

hist(out1[,2]/out1[,1]) #We are mantaining most of reads. Specially in second round. 

#Round 1: Total of 4058243 reads in, and 
#Round 2: Total of 14028683 reads in, and 13190164 out
```

Infer sequence variants:
```{r}
#Learn error rates
set.seed(100)
errF=learnErrors(filtFs,multithread = T,randomize = T)
errR=learnErrors(filtRs,multithread = T,randomize = T)
plotErrors(errF,nominalQ=F)

#Dereplicate
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

names(derepFs) <- sample.names
names(derepRs) <- sample.names

#Infer variants
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR,  multithread=TRUE)


```

Merge sequences and make sequence table
```{r}
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
head(mergers[[1]])
seqtab1=makeSequenceTable(mergers)
saveRDS(seqtab1,"~/Dropbox/Rebolleda_Gomez_AshmanLab/Mimulus/Analyses/DADA/seqtabR2.rds")
```

Run for both sequencing runs separately and then merge
```{r}
st1=readRDS("~/Dropbox/Rebolleda_Gomez_AshmanLab/Mimulus/Analyses/DADA/seqtabR1.rds")
st2=readRDS("~/Dropbox/Rebolleda_Gomez_AshmanLab/Mimulus/Analyses/DADA/seqtabR2.rds")
st.all=mergeSequenceTables(st1,st2)
saveRDS(st.all,"~/Dropbox/Rebolleda_Gomez_AshmanLab/Mimulus/Analyses/DADA/seqtabALL.rds")
```

Remove chimeras and small sequences
```{r}

#Check sequence lengths (must be around 250bp)
hist(nchar(getSequences(st.all)))

#Percentage of sequences longer than 256 --> 1.7%
sum(nchar(getSequences(st.all))>256)/dim(st.all)[2]*100
##Percentage of sequences shorter than 248 --> 6.5% !!!!
sum(nchar(getSequences(st.all))<248)/dim(st.all)[2]*100

#Remove sequences much shorter or longer than expected
seqtab2 <-st.all[,nchar(colnames(st.all)) %in% seq(248,256)]

seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
```

Assign OTUs
```{r}
taxa <- assignTaxonomy(seqtab.nochim, "/media/mariargz/Extra Drive 1/MimulusMicrobiome/raw/gg_13_8_train_set_97.fa.gz", multithread=TRUE)
#taxa <- assignTaxonomy(seqtab.nochim, "/media/mariargz/Extra Drive 1/MimulusMicrobiome/raw/silva_nr_v132_train_set.fa.gz", multithread=TRUE)
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
taxa.print.df=as.data.frame(taxa.print)
write.table(taxa,"taxa.txt")
```


Load as a phylloseq object to remove chloroplast
```{r}
metadata=read.table("~/Dropbox/Rebolleda_Gomez_AshmanLab/Mimulus/RAW/MIGU_metadata.txt", header = T, sep="\t")
metadata$SampleID=as.character(metadata$SampleID)
metadata=metadata[match(rownames(seqtab.nochim),(metadata$SampleID)),]
metadata$SampleID==rownames(seqtab.nochim)
rownames(metadata)=rownames(seqtab.nochim)

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(metadata), 
               tax_table(taxa))
no_chloroplast <- subset_taxa(ps,Class!="c__Chloroplast")
no_chloroplast <-subset_taxa(no_chloroplast,Family!="f__mitochondria")
```

Plot histogram of sample reads with or without mitochondria
```{r}
noeuk=rowSums(otu_table(no_chloroplast))
all=rowSums(otu_table(ps))

#Create data.frame with reads for plot
read.samp=data.frame(counts=c(noeuk,all),group=rep(c("no_euk","all"),each=length(all)))

library(ggplot2)
ggplot(read.samp, aes(x=log10(counts),fill=group))+
  geom_histogram(alpha=0.7,position = "identity", bins=50)+
  scale_fill_manual(values = c("olivedrab3","gray40"))+
  geom_vline(xintercept =log10(700), linetype=3)+
  theme_bw()

quantile(noeuk,probs=0.05) #Many with mostly eukaryotes
sum(noeuk<=1200) #11 samples!!! under 1100 reads

toolow=metadata[match(names(noeuk[noeuk<=1200]),metadata$SampleID),]
table(toolow$Organ)
#Anthers Community    Leaves    Petals      Seed      Soil  Style 
# 2         6         0         2         0         0        1
```

Remove rare ASVs
```{r}
GPr  = transform_sample_counts(no_chloroplast, function(x) x / sum(x) )
GPfr = filter_taxa(GPr, function(x) mean(x) > 1e-5, TRUE)
keeptaxa = taxa_names(GPfr)
alltaxa = taxa_names(ps)
myTaxa = alltaxa %in% keeptaxa
filter_ps=prune_taxa(myTaxa,ps)
otu_tab=otu_table(filter_ps)
```

Construct tree for dada2 data
```{r}
seqs <- colnames(otu_tab)
names(seqs) <- seqs # This propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)
phangAlign <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phangAlign)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phangAlign)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
        rearrangement = "stochastic", control = pml.control(trace = 0))

#loglikelihood: -247345.4 

#unconstrained loglikelihood: -1729.875 
#Proportion of invariant sites: 0.1679502 
#Discrete gamma model
#Number of rate categories: 4 
#Shape parameter: 0.530229 

#Rate matrix:
#          a         c         g        t
#a 0.0000000 0.9720363 3.2790317 1.507019
#c 0.9720363 0.0000000 0.7651643 4.504563
#g 3.2790317 0.7651643 0.0000000 1.000000
#t 1.5070187 4.5045627 1.0000000 0.000000

#Base frequencies:  
#0.2264066 0.1898124 0.3042643 0.2795168 
```

Save files
```{r}
write.tree(fitGTR$tree, file="~/Dropbox/Rebolleda_Gomez_AshmanLab/Mimulus/Analyses/DADA/Tree.tre")
write.table(otu_tab,file="~/Dropbox/Rebolleda_Gomez_AshmanLab/Mimulus/Analyses/DADA/OTU_noeuk.txt")
```
