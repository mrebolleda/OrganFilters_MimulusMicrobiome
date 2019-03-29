---
title: Sequence processing
---


##Paper: Floral organs act as environmental filters and interact with pollinators to structure the yellow monkeyflower (Mimulus guttatus) floral microbiome. 
###Code by: María Rebolleda Gómez
####Last updated: March 2018

##Data processing protocol
1. Download data and check quality
```{bash eval=FALSE}
parent_dir=/home/mariargz/Desktop/MimulusMicrobiome

cd $parent_dir
cd $parent_dir/raw

#Get the data from server (First round data)
wget ftp://ftp.igsb.anl.gov/jobs/d365c17ebcb41211aec43a9e4d354bb2/170816_RebolledaGomez_fastq/Undetermined_S0_L001_I1_001.fastq.gz
wget ftp://ftp.igsb.anl.gov/jobs/d365c17ebcb41211aec43a9e4d354bb2/170816_RebolledaGomez_fastq/Undetermined_S0_L001_R1_001.fastq.gz
wget ftp://ftp.igsb.anl.gov/jobs/d365c17ebcb41211aec43a9e4d354bb2/170816_RebolledaGomez_fastq/Undetermined_S0_L001_R2_001.fastq.gz



#Uncompress files 
gunzip *.gz

fastaF=$parent_dir/raw/Undetermined_S0_L001_R1_001.fastq
fastaQ=$parent_dir/raw/Undetermined_S0_L001_R2_001.fastq
map=$parent_dir/raw/MIGU_metadata.txt

#Check quality of reads
 ../code/FastQC/fastqc $input_dir/*.fastq

```


2. Use PEAR to merge pair end reads (Zhang et al., 2014, Bioinformatics. 30:614-620) and re-assign barcodes. 
```{bash eval=FALSE}
pear -f $fastaF -r $fastaR -o reads.join.PEAR | tee pear_assembly_log.txt
#Get code by Daniel Smith to re-assign the proper barcodes
wget https://www.dropbox.com/s/hk33ovypzmev938/fastq-barcode.pl?dl=1 

#With fasta-barcode remove barcodes that no longer have a corresponding read (code by Smith, Daniel)
$parent_dir/code/fastq-barcode.pl $parent_dir/raw/Undetermined_S0_L001_I1_001.fastq reads.join.PEAR.assembled.fastq > barcodes.join.PEAR.fastq

```

3. Demultiplex and remove reads with low score. High quality filtering Q20 and better
```{bash eval=FALSE}
source activate qiime1 
#Validate map
validate_mapping_file.py -m $map  -o validate_map

split_libraries_fastq.py -i reads.join.PEAR.assembled.fastq -m $map -b barcodes.join.PEAR.fastq -o demult_flowers_pear -q 19 --store_qual_scores --store_demultiplexed_fastq --barcode_type 12 
```

The sequencing facility had not sequenced the other half. They sent new results on different files. Re-ran the first steps in the same way as before but with the new files. To avoid overwriting data I changed the name of the old files as: "firstRun*.fastq"

EXTRA STEP: Concatenate our files and count how many sequences we were able to assemble.
```{bash eval=FALSE }
source activate qiime1
cat demult_flowers_pear/seqs.fna und_demult_flowers_pear/seqs.fna > demult
_flowers_pear_all.fna

mothur
#In mothur
summary.seqs(fasta=demult_flowers_pear_all.fna)


Output File Names: 
demult_flowers_pear_all.summary
# of Seqs:	15948447
#Output file: demult_flowers_pear_all.summary
summary.seqs(fasta=demult_flowers_pear/seqs.fna)
# of Seqs:	12416802
#Output file: demult_flowers_pear/seqs.summary
summary.seqs(fasta=und_demult_flowers_pear/seqs.fna)
# of Seqs:	3531645
#Output file: und_demult_flowers_pear/seqs.summary

```
```{r}
12416802+3531645==15948447 #All sequences are in the file now. 
```

4. Align sequences and classify into OTUs
```{bash  eval=FALSE}
cd $parent_dir/clean
source activate qiime1
pick_open_reference_otus.py -i demult_flowers_pear_all.fna -o uclust_greengenes/ -a -O 4 
```

Check number of failed to align sequences 
```{bash  eval=FALSE}
cd uclust_greengenes/pynast_aligned_seqs
count_seqs.py -i rep_set_failures.fasta
#2798  : rep_set_failures.fasta (Sequence lengths (mean +/- std): 235.7173 +/- 21.7937)
#2798  : Total
```
```{r}
2798/15948447*100
```

5. Remove low abundance counts in conservative manner (according to recommendation from Bokulich et al., 2013)/ We did all analyses without this filter first and obtained the same results. 
```{bash eval=FALSE}
filter_otus_from_otu_table.py -i otu_table_mc2_w_tax_no_pynast_failures.biom -o  otu_table_mc2_nofailures_filtered.biom --min_count_fraction 0.00005

biom summarize-table -i otu_table_mc2_nofailures_filtered.biom > summary_table_mc2_nofailures_filtered.txt
```

6. Remove chloroplast and mitochondria
```{bash eval=FALSE}
filter_taxa_from_otu_table.py -i otu_table_mc2_nofailures_filtered.biom -o otu_table_filtered_no_euks.biom -n c__Chloroplast,f__mitochondria

biom summarize-table -i otu_table_filtered_no_euks.biom > summary_table_filtered_no_euks.txt
```

Evaluate numbers of sequences before and after chloroplast removal
```{r}
#Read summary tables and plot histograms before and after removing euk
noeuk=read.table("/media/mariargz/Extra Drive 1/MimulusMicrobiome/clean/uclust_greengenes/summary_table_filtered_no_euks.csv",header = T,sep="\t",stringsAsFactors = F)
all.samp=read.table("/media/mariargz/Extra Drive 1/MimulusMicrobiome/clean/uclust_greengenes/summary_table_mc2_nofailures_filtered.csv",header = T,sep="\t",stringsAsFactors = F)

#Create data.frame with reads for plot
read.samp=data.frame(counts=c(noeuk[,2],all.samp[,2]),group=rep(c("no_euk","all"),each=length(all.samp[,2])))

library(ggplot2)
ggplot(read.samp, aes(x=log10(counts),fill=group))+
  geom_histogram(alpha=0.7,position = "identity", bins=50)+
  scale_fill_manual(values = c("olivedrab3","gray40"))+
  geom_vline(xintercept =log10(1100), linetype=3)+
  theme_bw()

quantile(noeuk[,2],probs=0.1)

```

