# PIPELINE MICROBIOTA ANALYSIS
Bioinformatic pipeline to analyze microbiota from Arbacia lixula and Paracentrotus lividus using different tissues.
## we do analysis in VENUS 

## connect to VENUS
ssh intern@161.116.67.159
password Marthasterias

### BUILD ENVIRONMENT AND INSTALL ###
module load anaconda3/2020.11

wget https://data.qiime2.org/distro/core/qiime2-2022.11-py38-linux-conda.yml
conda env create -n qiime2-2022.11 --file qiime2-2022.11-py38-linux-conda.yml
rm qiime2-2022.11-py38-linux-conda.yml
conda activate qiime2-2022.11



### STEP 1: IMPORT RAW FILES INTO QIIME2 ###
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path Microbiota_adults/ \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path demux-seqs.qza



### STEP 2: TRIM PRIMERS OFF RAW READS #
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences demux-seqs.qza
  --p-front-f CCTACGGGNGGCWGCAG \
  --p-front-r GACTACHVGGGTATCTAATCC \
  --p-error-rate 0 \
  --o-trimmed-sequences trimmed-seqs.qza



### STEP 3: JOIN/MERGE PAIRED-ENDS ###
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



### STEP 4: QUALITY CONTROL PAIRED-END READS ###
qiime quality-filter q-score-joined \
  --i-demux joined-seq.qza\
  --o-filtered-sequences filtered-joined-seqs.qza \
  --o-filter-stats filtered-joined-stats.qza \
  --p-quality-window 5 \
  --p-min-quality 25 



### STEP 5: REIMPORT FILTERED FILES WITH MODIFIED NAMES ###
	#NOTE: THIS IS OPTIONAL
qiime tools export \
  --input-path filtered-joined-seqs.qza
  --output-path Microbiota_adults/

qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path Microbiota_adults/ \
  --output-path filtered-joined-seqs2.qza \
  --input-format SingleEndFastqManifestPhred33



### STEP 6: DENOISE PROCESSED READS ###
qiime deblur denoise-16S \
  --i-demultiplexed-seqs filtered-joined-seqs2.qza \
  --p-trim-length 400 \
  --o-representative-sequences rep-seqs-deblur.qza \
  --o-table table-deblur.qza \
  --p-sample-stats \
  --o-stats stats-deblur.qza \

###  VISUALIZATION AND SUMMARY	
    #NOTE: THIS IS OPTIONAL
 
qiime feature-table summarize \
  --i-table table-deblur.qza \
  --o-visualization table-deblur.qzv \
  --m-sample-metadata-file sample-metadata.tsv
qiime feature-table tabulate-seqs \
  --i-data rep-seqs-deblur.qza \
  --o-visualization rep-seqs-deblur.qzv


### STEP 7: CREATE PHYLOGENY ###
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



### STEP 8: ASSIGN TAXONOMY WITH TRAINED CLASSIFIER ###
wget \
-O "silva-138-99-nb-classifier.qza" \
"https://data.qiime2.org/2022.11/common/silva-138-99-nb-classifier.qza"

qiime feature-classifier classify-sklearn \
  --i-classifier silva-138-99-nb-classifier.qza \
  --i-reads rep-seqs-deblur.qza \
  --o-classification taxonomy.qza \
  --p-confidence 0 \
  --p-read-orientation same \



### STEP 9: FILTER ARCHAEA, CHLOROPLAST, SINGLETONS ###
qiime taxa filter-table \
  --i-table table-deblur.qza \
  --i-taxonomy taxonomy.qza \
  --p-exclude Archaea \
  --p-mode contains \
  --o-filtered-table table-deblur-noAr.qza \

qiime taxa filter-table \
  --i-table table-deblur-noAr.qza
  --i-taxonomy taxonomy.qza
  --p-exclude Chloroplast \
  --p-mode contains \
  --o-filtered-table table-deblur-noArnoCh.qza \

qiime feature-table filter-features \
  --i-table table-deblur-noArnoCh.qza
  --p-min-frequency 2 \
  --o-filtered-table table-deblur-clean.qza \



### STEP 10: FILTER FEATURES FROM DNA SAMPLES ###
#IDENTIFY FEATURES IN DNA KIT BLANKS
qiime tools export \
  --input-path 
  --output-path 

biom convert \
  -i 
  -o 
  --to-tsv

#FILTER FEATURES FOUND IN DNA KIT BLANKS
qiime feature-table filter-features \
  --i-table 
  --m-metadata-file 
  --o-filtered-table 



### STEP 11: RAREIFY TABLE ###
qiime feature-table rarefy \
  --i-table 
  --p-sampling-depth 
  --p-with-replacement \
  --o-rarefied-table 
