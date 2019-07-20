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
### B. Merge the bed files
`/home/jacklyn/PycharmProjects/r4/peakCallScripts_orig/2_make_merged_bed_ORIG.sh`
### C. Create the Market matrix file
`/home/jacklyn/PycharmProjects/r4/peakCallScripts_orig/3_create_matrix_ORIG.sh`
### D. build csr int
#### Required files: 
Python script: `/home/jacklyn/PycharmProjects/r4/peakCallScripts_orig/4_build_csr_int_ORIG.py`
Merged Bed file: `Islet_123.combined.merged.bed`


## 5. Run Cicero

## 6. Annotate Data, sort co-accessible pairs




