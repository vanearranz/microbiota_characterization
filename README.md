# PIPELINE MICROBIOTA ANALYSIS
Bioinformatic pipeline to analyze microbiota from Arbacia lixula and Paracentrotus lividus using different body compartments.

## connect to VENUS
ssh intern@161.116.67.159

### BUILD ENVIRONMENT AND INSTALL ###
```
module load anaconda3/2020.11

wget https://data.qiime2.org/distro/core/qiime2-2022.11-py38-linux-conda.yml
conda env create -n qiime2-2022.11 --file qiime2-2022.11-py38-linux-conda.yml
rm qiime2-2022.11-py38-linux-conda.yml

conda activate qiime2-2022.11
```


### STEP 1: IMPORT RAW FILES INTO QIIME2 ###
```
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path Microbiota_adults/ \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path demux-seqs.qza
```


### STEP 2: TRIM PRIMERS OFF RAW READS #
```
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences demux-seqs.qza
  --p-front-f CCTACGGGNGGCWGCAG \
  --p-front-r GACTACHVGGGTATCTAATCC \
  --p-error-rate 0 \
  --o-trimmed-sequences trimmed-seqs.qza
```

### STEP 3: JOIN/MERGE PAIRED-ENDS ###
```
qiime vsearch join-pairs \
  --i-demultiplexed-seqs trimmed-seqs.qza\
  --o-joined-sequences joined-seqs.qza\
  --p-minovlen 20 \
  --p-maxdiffs 10 \
  --p-minmergelen 350 \
  --p-maxmergelen 550 \
  --p-allowmergestagger \
  --p-truncqual 10 \
  --p-minlen 100 \
  --p-qmax 41 \
  --p-qmaxout 41
```


### STEP 4: QUALITY CONTROL PAIRED-END READS ###
```
qiime quality-filter q-score-joined \
  --i-demux joined-seq.qza\
  --o-filtered-sequences filtered-joined-seqs.qza \
  --o-filter-stats filtered-joined-stats.qza \
  --p-quality-window 5 \
  --p-min-quality 25 
```


### STEP 5: REIMPORT FILTERED FILES WITH MODIFIED NAMES ###
```
#NOTE: THIS IS OPTIONAL
qiime tools export \
  --input-path filtered-joined-seqs.qza
  --output-path Microbiota_adults/

qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path Microbiota_adults/ \
  --output-path filtered-joined-seqs2.qza \
  --input-format SingleEndFastqManifestPhred33
```


### STEP 6: DENOISE AND CLUSTERING INTO ASV ###
```
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs trimmed-seqs.qza \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f 260 \
  --p-trunc-len-r 260 \
  --p-n-threads 4
  --o-representative-sequences rep-seqs-dada2.qza \
  --o-table table-dada2.qza \
  --o-denoising-stats denoising-stats-dada2.qza
```
####  VISUALIZATION AND SUMMARY	
```
qiime feature-table summarize \
--i-table table-dada2.qza \
--o-visualization table-dada2.qzv 

qiime feature-table tabulate-seqs \
--i-data rep-seqs-dada2.qza \
--o-visualization rep-seqs-dada2.qzv
```
#### CREATE FEATURE TABLE WITH METADATA 
```
qiime feature-table summarize \
  --i-table table-dada2.qza \
  --o-visualization table-dada2-summary.qzv \
  --m-sample-metadata-file METADATA_microbiota.tsv
```
### STEP 7: CREATE PHYLOGENY ###
```
#ALIGNMENT OF REPRESENTATIVE SEQUENCES
qiime alignment mafft \
  --i-sequences rep-seqs-deblur.qza \
  --o-alignment alignment-rep-seqs.qza \

#MASK HIGHLY VARIABLE NOISY POSITIONS IN ALIGNMENT
qiime alignment mask \
  --i-alignment alignment-rep-seqs.qza \
  --o-masked-alignment masked-alignment-rep-seqs.qza \

#CREATE PHYLOGENY WITH FASTTREE
qiime phylogeny fasttree \
  --i-alignment alignment-rep-seqs.qza \
  --o-tree fasttree-alignment-rep-seqs.qza \

#ROOT PHYLOGENY AT MIDPOINT
qiime phylogeny midpoint-root \
  --i-tree fasttree-alignment-rep-seqs.qza \
  --o-rooted-tree rooted-tree-rep-seqs.qza \
```

### STEP 8: ASSIGN TAXONOMY WITH TRAINED CLASSIFIER ###
```
wget \
-O "silva-138-99-nb-classifier.qza" \
"https://data.qiime2.org/2022.11/common/silva-138-99-nb-classifier.qza"

qiime feature-classifier classify-sklearn \
  --i-classifier silva-138-99-nb-classifier.qza \
  --i-reads rep-seqs-dada2.qza \
  --o-classification taxonomy.qza \
  --p-n-jobs 1
  --p-confidence 0.7 \
  --p-read-orientation same \
```


### STEP 9: FILTER ARCHAEA, CHLOROPLAST, SINGLETONS ###
```
qiime taxa filter-table \
  --i-table table-dada2.qza \
  --i-taxonomy taxonomy.qza \
  --p-exclude Archaea \
  --p-mode contains \
  --o-filtered-table table-ddada2-noAr.qza 

qiime taxa filter-table \
  --i-table table-clean.qza \
  --i-taxonomy taxonomy.qza \
  --p-exclude Eukaryota \
  --p-mode contains \
  --o-filtered-table table-dada2-clean.qza 

qiime feature-table filter-features \
  --i-table table-dada2-clean.qza \
  --p-min-frequency 2 \
  --o-filtered-table table-clean.qza 
```

### STEP 11: RAREIFY TABLE ###
```
qiime diversity alpha-rarefaction \
  --i-table table-dada2.qza \
  --i-phylogeny rooted-tree.qza \
  --p-max-depth 7000 \
  --m-metadata-file metadata_microbiota.tsv \
  --o-visualization alpha-rarefaction.qzv
```

### STEP 12: CREATE TAXA BARPLOT ### 
```
qiime taxa barplot \
  --i-table table-clean.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file METADATA_microbiota.tsv \
  --o-visualization taxa-bar-plots.qzv
``` 

### STEP 13: DIVERSITY ANALYSIS ###
```
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table ARB-table.qza \
  --p-sampling-depth 1103 \
  --m-metadata-file metadata_microbiota.tsv \
  --output-dir core-metrics-results
```
 #### ALPHA DIVERSITY ####
```
qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/faith_pd_vector.qza \
  --m-metadata-file metadata_microbiota.tsv \
  --o-visualization core-metrics-results/faith-pd-group-significanceboth.qzv
```
  #### BETA DIVERSITY ####
```  
qiime emperor plot \
  --i-pcoa core-metrics-results/bray_curtis_pcoa_results.qza \
  --m-metadata-file METADATA_microbiota.tsv \
  --o-visualization core-metrics-results/bray_curtis_pcoa_results.qzv
```

### STEP 14: EXPORT FROM QIIME2 TO R ###
Based on this [tutorial](https://github.com/pjtorres/Bioinformatics-BC/blob/master/Phyloseq/QIIME2_Phyloseq/Exporting%20QIIME2%20data%20for%20Phyloseq%20Analysis.ipynb)

```
# Export tree
qiime tools export --input-path unrooted-tree.qza\          --output-path Phyloseq
# Export taxonomy
qiime tools export --input-path taxonomy.qza\ 
--output-path Phyloseq
# Export table
qiime tools export --input-path table-clean.qza \
    --output-path Phyloseq
```
The first line in our taxonomy file must be changes to #OTUID

```
sed 's/Feature ID/#OTUID/' Phyloseq/taxonomy.tsv | sed 's/Taxon/taxonomy/' | sed 's/Confidence/confidence/' > Phyloseq/biom-taxonomy.tsv
```

Add the Taxonomy data to the Biom file 

```
biom add-metadata \
    -i Phyloseq/feature-table.biom \
    -o Phyloseq/table-with-taxonomyv2.biom \
    --observation-metadata-fp Phyloseq/biom-taxonomy.tsv \
    --sc-separated taxonomy 
```
Add the taxonomy data to your biom file
```
biom add-metadata \
    -i Phyloseq/feature-table.biom \
    -o Phyloseq/table-with-taxonomyv2.biom \
    --observation-metadata-fp Phyloseq/biom-taxonomy.tsv \
    --sc-separated taxonomy 
```

Now that we have the necessary files we will hop onto Phyloseq in R 
```
# Move the PhyloSeq folder to my local computer
scp -r intern@161.116.67.159:~/data/data/Microbiota_adults/Phyloseq /Users/vanessaarranz/Desktop/Microbiota adultos Lea/

# Add reference sequences for functional analisis
scp -r intern@161.116.67.159:~/data/data/Microbiota_adults/Phyloseq/dna-sequences.fasta /Users/vanessaarranz/Desktop/Microbiota_adultos_Lea/Phyloseq
```

### STEP 15 : FAPROTAX for functional annotations
```
Download the pacakge 
http://www.loucalab.com/archive/FAPROTAX/lib/php/index.php?section=Download

Directory /Users/vanessaarranz/Desktop/Microbiota_adultos_Lea/Phyloseq/FAPROTAX_1.2.10

# Convert biom table to tsv, because with biom didn't work
biom convert -i table-with-taxonomyv2.biom -o table-with-taxonomyv2.txt --to-tsv --header-key taxonomy

# Run command
python3	./collapse_table.py -i table-with-taxonomyv2.txt -g FAPROTAX.txt -f -o functional_otu_table.tsv -r report.txt --column_names_are_in last_comment_line  --keep_header_comments --non_numeric consolidate -v --row_names_are_in_column "taxonomy" --omit_columns 0 --normalize_collapsed columns_before_collapsing --group_leftovers_as 'other'
```

### STATISTICAL ANALYSIS AND VISUALIZATIONS IN R 

Statistical analysis and plots 
**[Statiscal analysis and visualizations](https://github.com/vanearranz/microbiota_characterization/tree/main#:~:text=12%20Commits-,Functional_predictions2.Rmd,-Create%20Functional_predictions2.Rmd)**: Use the R Script

Functional analysis including statistical comparisons with FAPROTAX database
**[Fuctional analysis](https://github.com/vanearranz/microbiota_characterization/tree/main#:~:text=now-,microbiome_Adult_seaurchin2.Rmd,-Create%20microbiome_Adult_seaurchin2.Rmd)** 

All the results are included in the manuscript "XXXXX" by Vanessa Arranz et al. 

