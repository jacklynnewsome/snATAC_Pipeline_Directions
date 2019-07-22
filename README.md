# Gaulton Lab snATAC Pipeline Directions
### Required software for entire pipeline:
	Python 3
	Anaconda, probably
## 1. Run Cell Ranger on the raw FASTQ snATAC files
(https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger)  

### Commands: 
	``` I don't actually know what this pipeline entails. I have never run this ```
### Outputs:
Bam file
others? 


## 2. Run the 10xGenomics processing script that Josh wrote

### Commands:
```python snATAC_pipeline_10X.py -b [cell ranger output bam file] -o [output directory ? ] -n [project name prefix] -t 24 -m 2 ```  


### Example:
``` /home/data/anaconda3/bin/python /nfs/lab/projects/pbmc_snATAC/scripts/snATAC_pipeline_10X.py -b cellranger_v1.1/possorted_bam.bam -o lab_pipeline -n pbmc1 -t 24 -m 2 ```  

### Required files: 
Cell ranger output bam file: `cellranger_v1.1/outs/possorted_bam.bam `     
Josh's 10x pipeline script ('snATAC_pipeline_10X.py'): `/nfs/lab/projects/pbmc_snATAC/scripts/snATAC_pipeline_10X.py`   

### Required software:   
##### Python Packages:
	pysam
	numpy
	pandas
	scipy
	multiprocessing


## 3. Filter low quality and doublet cells out and cluster cells by type
### Required files: 
Jupyter notebook: `http://localhost:8888/notebooks/Documents/R4/pbmc_snATAC_islet.ipynb`   
#### A working directory containing 'snATAC_pipeline_10X.py' outputs:
```
sample_name.barcodes
sample_name.long_fmt_mtx.txt.gz
sample_name.mtx.gz
sample_name.narrowPeak
sample_name.qc_metrics.txt
sample_name.regions
```

The name of the sample, eg "sample_name" or "SAMN1120500120-hi-cyto"   
A reference file of promoter names and locations: `/home/jacklyn/Documents/R4/gencode.v19.5kb_promoter_names.txt`   
The names and abbreviations of all cell type markers relevant to the sample   
(optional) An mp3 to serve as an alert for when a jupyter notebook cell is finished running: `/home/jacklyn/Documents/R4/sweetalertsound5.wav`   
### Required software:   
	Jupyter notebooks
##### Python Packages:
	numpy 
	pandas
	scanpy
	matplotlib
	seaborn
	statsmodels.api
	sklearn
	scipy
	IPython.display
	anndata
	rpy2.ipython

### Outputs: 
Figures showing clustering before, during, and after quality filtration
A file with cells labelled with cell types, based on leiden cluster: ` (example)  Islet.cluster_labels.txt`


### Directions: 
Follow the directions in the Jupyter notebook to perform this step







## 4. Call Peaks
### A. Call Cluster Peaks
#### Required software:   
`bedGraphToBigWig`
`macs2`

#### Required files: 
Python script: `/home/jacklyn/PycharmProjects/r4/peakCallScripts_orig/1_call_cluster_peaks_ORIG.py`  
Tag align file: `(example) pbmc1.tagalign.gz`  
barcode cluster file (output by step 3): ` (example)  Islet.cluster_labels.txt`  
Reference chromosome sizes file: `/home/joshchiou/references/hg19.chrom.sizes`  
#### Other: 
Label prefix: `(example) pbmc1`  
#### Output:
bdg file : `(example) pbmc1_treat_pileup.bdg`  
bdg file: `(example) pbmc1_.norm.bdg`  
bed file : `(example) pbmc1_summits.bed`  
xls file: `(example) _peaks.xls`  

#### Directions: 
command: `python 1_call_cluster_peaks_ORIG.py [tag align file name] [barcode cluster file] [sample prefix]`  
example: `python 1_call_cluster_peaks_ORIG.py pbmc1.tagalign.gz pbmc1.cluster_labels.txt pbmc1`  

### B. Merge the bed files
#### Required software
`bedtools`
#### Required files: 
Bash script: `/home/jacklyn/PycharmProjects/r4/peakCallScripts_orig/2_make_merged_bed.sh`  
ENCODE blacklist file: `/home/jnewsome/references/ENCODE.hg19.blacklist.bed`
#### Other:
A list of the cell types that were used to label the clusters / barcodes in step 3
narrow peak files from step 3: `(example) /home/jnewsome/pipeline_islet/peakCalls/narrowPeak/islet.Alpha_cell_peaks.narrowPeak`  
#### Output: 
bed file: `(example) islet.Alpha_cell_peaks.bed`  
bed file: `islet.bed`  
bed file: `islet.sorted.bed`
bed file: `islet.sorted.filtered.bed`
bed file: `islet.sorted.merged.bare.bed`  

### C. Create the Market matrix file
#### Required software
`bedtools`
#### Required files: 
Deduplicated tag align file: `(example) pbmc1.filt.rmdup.tagAlign.gz`    
sorted merged bed file: `(example) pbmc1.sorted.merged.bed`   
Bash script: `/home/jacklyn/PycharmProjects/r4/peakCallScripts_orig/3_create_matrix_ORIG.sh`   
#### Output:
market matrix file for the merged peaks: `(example) pbmc1.merged_peaks.long_fmt.mtx.gz`  
### D. build csr int
#### Required software
`bedtools`
##### Python Packages:
`pandas`  
`scipy`  
#### Required files: 
Python script: `/home/jacklyn/PycharmProjects/r4/peakCallScripts_orig/4_build_csr_int_ORIG.py`  
Merged Bed file: `(example) Islet_123.combined.merged.bed`  
Deduplicated tag align file: `(example) pbmc1.filt.rmdup.tagAlign.gz`   
barcode cluster file (output by step 3): ` (example)  Islet.cluster_labels.txt`   
#### Other:
a prefix for the output files  
#### Output:
market matrix file for the merged peaks: `(example) pbmc1.merged_peaks.long_fmt.mtx.gz`   
barcodes file:	` (example) wtvr.barcodes`  
regions file:	` (example) wtvr.peaks`    
mtx file: ` (example) wtvr.mtx  `  
csr.npz file: `(example) wtvr.csr.npz`  
#### Directions:
WIP
## 5. Run Cicero
### Required software
`R`
#### R tools
`cicero`  
`stringr` 
`parallel`  
`data.table`
### Required Files 
Either :   
    Paula's cicero execution script: `cicero_islet1_remote.R`  
Or :   
    Paula's cicero execution jupyter notebook:  `/notebooks/pbmcNotebooksOrigCopy/Cicero_pbmc1_test.ipynb`
    
The annotated merged bed file: `(example) /home/jnewsome/pipeline_islet/peakCalls/islet.sorted.merged.bed`  
The cluster Labels file: `/home/jnewsome/pipeline_islet/Islet.cluster_labels.txt` 
The peaks long fmt market matrix file: `/home/jnewsome/pipeline_islet/peakCalls/islet.merged_peaks.long_fmt.mtx.gz`   
### Output:
rds file for each cell type: `(example) pbmc1.CELL_TYPE.1MB_cicero_conns.rds`  
txt file for each cell type:  `(example) pbmc1.CELL_TYPE.cicero_conns.txt`  
deduplicated txt file for each cell type: `(example) pbmc1.CELL_TYPE.cicero_conns_dedup.txt`  
### Directions:
 replace relevant information in the notebook or script   
## 6. Annotate Data, sort co-accessible pairs
### A. Intersect Reference promoter file with All of the peaks present in the data set
#### Required files:
Referance 1KB promoter reference file: `gencode.v19.1kb_all_possible_transcripts.bed`  
the merged peaks cell-type-annotated sorted bed file: `(example) /mnt/4d60fd49-d4ad-42d2-ac64-5b3f0265b9c1/snATAC_lung/lung.merged_peaks.anno.sorted.bed`  
#### Required Software:
`bedtools`  
#### Output Files:  
A bedfile with the intersected data set peaks, reference promoter names, and reference promoter peaks: `(example) /mnt/4d60fd49-d4ad-42d2-ac64-5b3f0265b9c1/r4_snATAC_islet/peakCalls/islet.sorted.merged.promoterAnnotated.refLocIncluded.bed`  
A bedfile with the reference promoter names and data set peaks: `(example) /mnt/4d60fd49-d4ad-42d2-ac64-5b3f0265b9c1/r4_snATAC_islet/peakCalls/islet.sorted.merged.promoterAnnotated.bed`  
#### Directions:  
cmd: `bedtools intersect -a [MERGED BED FILE] -b [REFERENCE FILE] -wa -wb > [REFERENCE / DATASET INTERSECTED BED FILE]`  
cmd: `awk '{print $1,$2,$3,$4,$8}' [REFERENCE / DATASET INTERSECTED BED FILE]  > [ANNOTATED DATASET PEAK BED FILE]`  
### B. Annotated Cicero Output File
#### Required files:  
Cicero output annotation script (JN): `/home/jacklyn/PycharmProjects/r4/annotateAndSortCiceroOutput.py`  
Cicero output deduplicated txt file for each cell type: `(example) pbmc1.CELL_TYPE.cicero_conns_dedup.txt`  
#### Output:
Annotated cicero output files:
```(example):
/mnt/4d60fd49-d4ad-42d2-ac64-5b3f0265b9c1/snATAC_islet/ciceroAnnotated/islet.Alpha_cell.cicero_conns_dedup.promAnno.bed  
/mnt/4d60fd49-d4ad-42d2-ac64-5b3f0265b9c1/snATAC_islet/ciceroAnnotated/islet.Alpha_cell.cicero_conns_dedup.promAnno.bed.AllPairAnnotations.bed  
/mnt/4d60fd49-d4ad-42d2-ac64-5b3f0265b9c1/snATAC_islet/ciceroAnnotated/islet.Alpha_cell.cicero_conns_dedup.promAnno.bed.CC_CoA.bed  
/mnt/4d60fd49-d4ad-42d2-ac64-5b3f0265b9c1/snATAC_islet/ciceroAnnotated/islet.Alpha_cell.cicero_conns_dedup.promAnno.bed.CC_NA.bed  
/mnt/4d60fd49-d4ad-42d2-ac64-5b3f0265b9c1/snATAC_islet/ciceroAnnotated/islet.Alpha_cell.cicero_conns_dedup.promAnno.bed.CC_Neg.bed  
/mnt/4d60fd49-d4ad-42d2-ac64-5b3f0265b9c1/snATAC_islet/ciceroAnnotated/islet.Alpha_cell.cicero_conns_dedup.promAnno.bed.CC_Zero.bed  
/mnt/4d60fd49-d4ad-42d2-ac64-5b3f0265b9c1/snATAC_islet/ciceroAnnotated/islet.Alpha_cell.cicero_conns_dedup.promAnno.bed.CP_CoA.bed  
/mnt/4d60fd49-d4ad-42d2-ac64-5b3f0265b9c1/snATAC_islet/ciceroAnnotated/islet.Alpha_cell.cicero_conns_dedup.promAnno.bed.CP_NA.bed  
/mnt/4d60fd49-d4ad-42d2-ac64-5b3f0265b9c1/snATAC_islet/ciceroAnnotated/islet.Alpha_cell.cicero_conns_dedup.promAnno.bed.CP_Neg.bed  
/mnt/4d60fd49-d4ad-42d2-ac64-5b3f0265b9c1/snATAC_islet/ciceroAnnotated/islet.Alpha_cell.cicero_conns_dedup.promAnno.bed.CP_Zero.bed  
/mnt/4d60fd49-d4ad-42d2-ac64-5b3f0265b9c1/snATAC_islet/ciceroAnnotated/islet.Alpha_cell.cicero_conns_dedup.promAnno.bed.idententicalPromoter_PP.bed  
/mnt/4d60fd49-d4ad-42d2-ac64-5b3f0265b9c1/snATAC_islet/ciceroAnnotated/islet.Alpha_cell.cicero_conns_dedup.promAnno.bed.PP_CoA.bed  
/mnt/4d60fd49-d4ad-42d2-ac64-5b3f0265b9c1/snATAC_islet/ciceroAnnotated/islet.Alpha_cell.cicero_conns_dedup.promAnno.bed.PP_NA.bed  
/mnt/4d60fd49-d4ad-42d2-ac64-5b3f0265b9c1/snATAC_islet/ciceroAnnotated/islet.Alpha_cell.cicero_conns_dedup.promAnno.bed.PP_Neg.bed  
/mnt/4d60fd49-d4ad-42d2-ac64-5b3f0265b9c1/snATAC_islet/ciceroAnnotated/islet.Alpha_cell.cicero_conns_dedup.promAnno.bed.PP_Zero.bed  
/mnt/4d60fd49-d4ad-42d2-ac64-5b3f0265b9c1/snATAC_islet/ciceroAnnotated/islet.Alpha_cell.cicero_conns_dedup.promAnno.bed.STATS.txt  
```

