# Gaulton Lab snATAC Pipeline Directions
### Required software for entire pipeline:
	Python 3
	Anaconda, probably
## 1. Run Cell Ranger on the raw FASTQ snATAC files
(https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger)  
### Required files: 
Directory with Input FASTQ Files
Cell Ranger Application: `/home/ysun/cellranger-atac-1.1.0/cellranger-atac count`
Reference Genome file: `/home/ysun/refdata-cellranger-atac-hg19-1.1.0`
### Other requirements: 
sample name (must match beginning of all FASTQ file names before the S1_ thing), example: `JYH_PBMC1`
id name to use, example: `PBMC1`
number of cores to use
### Commands: 
`[cell ranger application] --id=[id name] --fastq=[FASTQ directory] --sample=[sample name] --reference=[reference genome file] --localcores=[number of cores to use]`
example: `/home/ysun/cellranger-atac-1.1.0/cellranger-atac count --id=PBMC1 --fastq=/nfs/lab/projects/pbmc_snATAC/pipeline/snATAC/completeSamples/dataBackup/PBMC1 --sample=JYH_PBMC1 --reference=/home/ysun/refdata-cellranger-atac-hg19-1.1.0 --localcores=12`

### Outputs:
```
- Per-barcode fragment counts & metrics:        /nfs/lab/projects/pbmc_snATAC/pipeline/snATAC/completeSamples/cellRanger/PBMC3/outs/singlecell.csv
- Position sorted BAM file:                     /nfs/lab/projects/pbmc_snATAC/pipeline/snATAC/completeSamples/cellRanger/PBMC3/outs/possorted_bam.bam
- Position sorted BAM index:                    /nfs/lab/projects/pbmc_snATAC/pipeline/snATAC/completeSamples/cellRanger/PBMC3/outs/possorted_bam.bam.bai
- Summary of all data metrics:                  /nfs/lab/projects/pbmc_snATAC/pipeline/snATAC/completeSamples/cellRanger/PBMC3/outs/summary.json
- HTML file summarizing data & analysis:        /nfs/lab/projects/pbmc_snATAC/pipeline/snATAC/completeSamples/cellRanger/PBMC3/outs/web_summary.html
- Bed file of all called peak locations:        /nfs/lab/projects/pbmc_snATAC/pipeline/snATAC/completeSamples/cellRanger/PBMC3/outs/peaks.bed
- Raw peak barcode matrix in hdf5 format:       /nfs/lab/projects/pbmc_snATAC/pipeline/snATAC/completeSamples/cellRanger/PBMC3/outs/raw_peak_bc_matrix.h5
- Raw peak barcode matrix in mex format:        /nfs/lab/projects/pbmc_snATAC/pipeline/snATAC/completeSamples/cellRanger/PBMC3/outs/raw_peak_bc_matrix
- Directory of analysis files:                  /nfs/lab/projects/pbmc_snATAC/pipeline/snATAC/completeSamples/cellRanger/PBMC3/outs/analysis
- Filtered peak barcode matrix in hdf5 format:  /nfs/lab/projects/pbmc_snATAC/pipeline/snATAC/completeSamples/cellRanger/PBMC3/outs/filtered_peak_bc_matrix.h5
- Filtered peak barcode matrix in mex format:   /nfs/lab/projects/pbmc_snATAC/pipeline/snATAC/completeSamples/cellRanger/PBMC3/outs/filtered_peak_bc_matrix
- Barcoded and aligned fragment file:           /nfs/lab/projects/pbmc_snATAC/pipeline/snATAC/completeSamples/cellRanger/PBMC3/outs/fragments.tsv.gz
- Fragment file index:                          /nfs/lab/projects/pbmc_snATAC/pipeline/snATAC/completeSamples/cellRanger/PBMC3/outs/fragments.tsv.gz.tbi
- Filtered tf barcode matrix in hdf5 format:    /nfs/lab/projects/pbmc_snATAC/pipeline/snATAC/completeSamples/cellRanger/PBMC3/outs/filtered_tf_bc_matrix.h5
- Filtered tf barcode matrix in mex format:     /nfs/lab/projects/pbmc_snATAC/pipeline/snATAC/completeSamples/cellRanger/PBMC3/outs/filtered_tf_bc_matrix
- Loupe Cell Browser input file:                /nfs/lab/projects/pbmc_snATAC/pipeline/snATAC/completeSamples/cellRanger/PBMC3/outs/cloupe.cloupe
- csv summarizing important metrics and values: /nfs/lab/projects/pbmc_snATAC/pipeline/snATAC/completeSamples/cellRanger/PBMC3/outs/summary.csv
- Annotation of peaks with genes:               /nfs/lab/projects/pbmc_snATAC/pipeline/snATAC/completeSamples/cellRanger/PBMC3/outs/peak_annotation.tsv

```


## 2. Run the 10xGenomics processing script that Josh wrote
### Required files: 

Josh's 10x pipeline script: `/nfs/lab/projects/pbmc_snATAC/scripts/snATAC_pipeline_10X.py`
Cell ranger output Position sorted BAM file: `(example) /nfs/lab/projects/pbmc_snATAC/pipeline/snATAC/completeSamples/cellRanger/PBMC3/outs/possorted_bam.bam`
BAI file for the Cell ranger output Position sorted BAM file: `(example) /nfs/lab/projects/pbmc_snATAC/pipeline/snATAC/completeSamples/cellRanger/PBMC3/outs/possorted_bam.bam.bai`
### Required software:   
##### Python Packages:
	pysam
	numpy
	pandas
	scipy
	multiprocessing
### Other requirements: 
sample name to output with
output folder name
### Commands: 
`python [10x pipeline script] -b [bam file] -o [output folder name] -n [sample namer] -t 24 -m 2`
example: `/home/data/anaconda3/bin/python /nfs/lab/projects/pbmc_snATAC/scripts/snATAC_pipeline_10X.py -b /nfs/lab/projects/pbmc_snATAC/pipeline/snATAC/completeSamples/cellRanger/PBMC6/outs/possorted_bam.bam -o lab_pipeline -n pbmc6_completeSample -t 24 -m 2`
### Outputs:
```
pbmc11_completeSample.barcodes
pbmc11_completeSample.filt.md.bam
pbmc11_completeSample.filt.rmdup.bam
pbmc11_completeSample.long_fmt_mtx.txt.gz
pbmc11_completeSample_peaks.narrowPeak
pbmc11_completeSample.regions
pbmc11_completeSample.compiled.filt.bam
pbmc11_completeSample.filt.md.bam.bai
pbmc11_completeSample.filt.rmdup.tagAlign.gz
pbmc11_completeSample.mtx.gz
pbmc11_completeSample.qc_metrics.txt
```


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
    Paola's cicero execution script: `cicero_islet1_remote.R`  
Or :   
    Paola's cicero execution jupyter notebook:  `/notebooks/pbmcNotebooksOrigCopy/Cicero_pbmc1_test.ipynb`
    
The annotated merged bed file: `(example) /home/jnewsome/pipeline_islet/peakCalls/islet.sorted.merged.bed`  
The cluster Labels file: `/home/jnewsome/pipeline_islet/Islet.cluster_labels.txt` 
The peaks long fmt market matrix file: `/home/jnewsome/pipeline_islet/peakCalls/islet.merged_peaks.long_fmt.mtx.gz`   
### Output:
rds file for each cell type: `(example) pbmc1.CELL_TYPE.1MB_cicero_conns.rds`  
txt file for each cell type:  `(example) pbmc1.CELL_TYPE.cicero_conns.txt`  
deduplicated txt file for each cell type: `(example) pbmc1.CELL_TYPE.cicero_conns_dedup.txt`  
### Directions:
 replace relevant information in the notebook or script   
 Alternative command, implement later: `Rscript cicero_islet1_remote.R  /home/jnewsome/pipeline_islet/ciceroOutput islet.sorted.merged.fixed.bed Islet.cluster_labels.txt islet.merged_peaks.long_fmt.mtx islet`   
## 6. Annotate Data, sort co-accessible pairs
### A.v2. Make reference file for all peaks in dataset
#### 1. concatenate all of the by-chromosome cicero files together for each cell type
#### 2. use python script to make a file to intersect with the 1KB promoter reference file and a reference file with all of the pair ids for that file
##### Inputs:
merged cicero file, example: `/nfs/lab/projects/pbmc_snATAC/pipeline/snATAC/cicero/macro/merged/pbmc1to12.Macro.ALL_CELL_TYPES_MERGED.cicero_conns_dedup.txt`
##### Outputs:
ID reference file: `[specified output prefix]_IDReference.bed`
ID individual peaks with id: `[specified output prefix]_IndividualPeaks_withID.bed`

##### Script:
peakRefMakerScript = `/nfs/lab/projects/pbmc_snATAC/pipeline/snATAC/cicero/macro/createCiceroPairReferenceFileForBedtoolsIntersect.py`
##### Command:
`${pyBin} ${peakRefMakerScript} -i ${catAllCiceroFile} -o ${peakRefFileAlmost} -p ${refPrefix}`
### A. Intersect Reference promoter file with All of the peaks present in the data set
#### Required files:
Reference 1KB promoter reference file: `gencode.v19.1kb_all_possible_transcripts.bed`  
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
bed file from step 6.A: `(example) /mnt/4d60fd49-d4ad-42d2-ac64-5b3f0265b9c1/r4_snATAC_islet/peakCalls/islet.sorted.merged.promoterAnnotated.bed`  
#### Other: 
a threshold for Cicero coaccessibility score, as a decimal. eg, 5% = 0.05
#### Output:
General output: `/mnt/4d60fd49-d4ad-42d2-ac64-5b3f0265b9c1/snATAC_islet/ciceroAnnotated/islet.Alpha_cell.cicero_conns_dedup.promAnno.bed`  
Sorted annotated cicero output files:
```  
/mnt/4d60fd49-d4ad-42d2-ac64-5b3f0265b9c1/snATAC_islet/ciceroAnnotated/islet.Alpha_cell.cicero_conns_dedup.promAnno.bed.AllPairAnnotations.bed  
/mnt/4d60fd49-d4ad-42d2-ac64-5b3f0265b9c1/snATAC_islet/ciceroAnnotated/islet.Alpha_cell.cicero_conns_dedup.promAnno.bed.PP_NA.bed  
/mnt/4d60fd49-d4ad-42d2-ac64-5b3f0265b9c1/snATAC_islet/ciceroAnnotated/islet.Alpha_cell.cicero_conns_dedup.promAnno.bed.PP_Neg.bed  
/mnt/4d60fd49-d4ad-42d2-ac64-5b3f0265b9c1/snATAC_islet/ciceroAnnotated/islet.Alpha_cell.cicero_conns_dedup.promAnno.bed.PP_Zero.bed   
```   
Statistics file : 
`/mnt/4d60fd49-d4ad-42d2-ac64-5b3f0265b9c1/snATAC_islet/ciceroAnnotated/islet.Alpha_cell.cicero_conns_dedup.promAnno.bed.STATS.txt `  
#### Directions:  
cmd: `python annotateAndSortCiceroOutput.py -i [INPUT ANNOTATED CICERO FILE] -p [ANNOTATED PEAK REFERENCE BED FILE] -o [NAME OF GENERAL OUTPUT BED FILE] -t [SCORE THRESHOLD, eg '0.05']`  
Run this command for each cell type Cicero output deduplicated txt file  
#### Annotation Definitions:    
CC = CRE/CRE relationship (neither peak in cicero pair falls within a promoter region, as defined by the 1 kb reference promoter file)   
CP = CRE/Promoter relationship (1 peak in cicero pair falls within a promoter region)  
PP = Promoter/Promoter relationship (both peaks in cicero pair fall within promoter region(s))  
Scores:  
Zero = Cicero Coacessibility score for this pair of peaks is less than the selected threshold, but more than the negative threshold  
Neg = Cicero Coacessibility score for this pair of peaks is less or equal to the negative threshold  
CoA = Cicero Coacessibility score for this pair of peaks is greater than or equal to the selected threshold  
NA = The Cicero Coaccessibility score for this pair was marked as 'NA'  
## 7. Merge Bed files, a step that should be done before this but it's not, because that would mean I'd have to change the annotation script
#### Required files:
All of the "AllPairAnnotations.bed" files from the step before this one   
annotated merging cicero script : `/home/jacklyn/PycharmProjects/r4/mergeCiceroBedFilesAfterAnnotation.py` 
#### Other:
file prefix for the input files you need : `(example) islet.`   
file suffix for the input files you need : `(example) .cicero_conns_dedup.promAnno.bed.AllPairAnnotations.bed`  
the path to the directory with these files
#### Output Files:  
a merged "AllPairAnnotations.bed" file for the entire data set `(example) /mnt/4d60fd49-d4ad-42d2-ac64-5b3f0265b9c1/snATAC_islet/ciceroAnnotated/islet.ALL_CELL_TYPES.cicero_conns_dedup.promPairCell_annotated.bed`

#### Directions:  
cmd: python mergeCiceroBedFilesAfterAnnotation.py -i 

example: `python [MERGE SCRIPT] -i [INPUT DIRECTORY with the last forward slash] -o /mnt/4d60fd49-d4ad-42d2-ac64-5b3f0265b9c1/snATAC_islet/ciceroAnnotated/islet.ALL_CELL_TYPES.cicero_conns_dedup.promPairCell_annotated.bed -p islet. -s .cicero_conns_dedup.promAnno.bed.AllPairAnnotations.bed`

## 8. Get Unique Gene List
#### Required Files:
The merged promoter annotated bed file for the data set: `(example) /mnt/4d60fd49-d4ad-42d2-ac64-5b3f0265b9c1/snATAC_lung/snATAC_lung.cicero_conns_dedup.promPairCellAnno_ALL.bed`  
#### Output Files:
A list of genes from the 1st peak in the Cicero pairs: `(example) /mnt/4d60fd49-d4ad-42d2-ac64-5b3f0265b9c1/snATAC_lung/snATAC_lung.genesColA.txt`  
A list of genes from the 2nd peak in the Cicero pairs: `(example) /mnt/4d60fd49-d4ad-42d2-ac64-5b3f0265b9c1/snATAC_lung/snATAC_lung.genesColB.txt`  
A list of all genes in data set, with duplicates: `(example) /mnt/4d60fd49-d4ad-42d2-ac64-5b3f0265b9c1/snATAC_lung/snATAC_lung.genesColBoth.txt`  
A list of all genes in data set, without duplicates: `(example) /mnt/4d60fd49-d4ad-42d2-ac64-5b3f0265b9c1/snATAC_lung/snATAC_lung.genesColBoth_unique.txt`  
#### Directions:
cmd: `awk '{a[$4]++;} END{for(i in a) print a[i]"  "i}' [ANNOTATED_MERGED_BED_FILE] > [GENES_A_TXT]`  
`awk '{a[$5]++;} END{for(i in a) print a[i]"  "i}' [ANNOTATED_MERGED_BED_FILE] > [GENES_B_TXT]`
`cat [GENES_A_TXT] [GENES_B_TXT] > [GENES_BOTH_TXT]`
`awk '{a[$2]++;} END{for(i in a) print a[i]"  "i}' [GENES_BOTH_TXT] > [GENES_BOTH_UNIQUE_TXT]`
`awk '{print $2}' [GENES_BOTH_UNIQUE_TXT]  > [GENES_BOTH_UNIQUE_TXT]`  

## 9. (Optional) Trim Gene List to Relevant Entries

## 9. Annotate with pathways
#### Required Files:
A list of all promoter genes in data set: `(example) /mnt/4d60fd49-d4ad-42d2-ac64-5b3f0265b9c1/snATAC_lung/snATAC_lung.genesColBoth_unique.txt`   
Script for generating gene-pathway file from REACTOME analysis:
`/home/jacklyn/PycharmProjects/r4/createReferenceFromReactomeAnalysis.py`  
Script for annotating annotated Cicero bed file with pathways: `/home/jacklyn/PycharmProjects/r4/annotateWithPathways.py`  OR `/home/jacklyn/PycharmProjects/r4/annotateCiceroPromotersWithInformation.py`  
#### Other: 
access to the REACTOME web tool: `https://reactome.org/PathwayBrowser/#TOOL=AT`  
#### Output: 
REACTOME analysis output file: `(example) snATAC_lung.reactome_analysis_result.csv`
Cicero bed file, annotated with the pathways, and whatever else you've already done. it's like 9 PM. I'm tired : `(example) snATAC_lung.cicero_conns_dedup.promPairCellAnno_ALL.bed`
#### Directions:
Copy and paste the gene list into the REACTOME analysis tool. Click continue, Select "Project to human", and click "Analyze!" In the table at the bottom of the screen, click on the tab that says "Downloads" in the table (below "Results" and "Not Found"). Download the "Pathway analysis results" and whatever else ( I haven't used anything else yet).
cmd: `python [REACTOME TO GENE PATHWAYS SCRIPT] -i [REACTOME ANALYSIS FILE] -o [GENE-PATHWAY LIST]`
`python [PATHWAY ANNOTATION SCRIPT] -i [CICERO BED FILE] -o [CICERO BED FILE, NOW WITH PATHWAYS]`

## 10. Annotated with noted sites of interest (Islet Only?)
### A. Intersect Differential Site Bed file with Data set Peak Bed file:
#### Required files:
Differentially expressed sites of interest : https://github.com/anthony-aylward/islet-cytokines-outline/blob/master/differential-sites.bed
the merged peaks cell-type-annotated sorted bed file: `(example) /mnt/4d60fd49-d4ad-42d2-ac64-5b3f0265b9c1/snATAC_islet/peakCalls/islet.sorted.merged.bed`  
#### Required Software:
`bedtools`  
#### Output Files:  
A bedfile with the intersected data set peaks, reference promoter names, and reference promoter peaks: `(example) /mnt/4d60fd49-d4ad-42d2-ac64-5b3f0265b9c1/r4_snATAC_islet/peakCalls/islet.sorted.merged.promoterAnnotated.refLocIncluded.bed`  
A bedfile with the reference promoter names and data set peaks: `(example) /mnt/4d60fd49-d4ad-42d2-ac64-5b3f0265b9c1/r4_snATAC_islet/peakCalls/islet.sorted.merged.promoterAnnotated.bed`  
#### Directions:  
cmd: `bedtools intersect -a [MERGED BED FILE] -b [REFERENCE FILE] -wa -wb > [REFERENCE / DATASET INTERSECTED BED FILE]`  
cmd: `awk '{print $1,$2,$3,$4,$8}' [REFERENCE / DATASET INTERSECTED BED FILE]  > [ANNOTATED DATASET PEAK BED FILE]`    

## 12. Annotate with gene families?
#### Required Files:
Merged, annotated cicero bed file
Unique gene list
annotation script
#### output: 
a gene-family list text file
a gene- definition list text file
yet another version of the annotated cicero file
#### Other: 
Access to Panther Db: `http://pantherdb.org/`
#### Directions:
upload the gene list to Panther Db. download the table. 
cmd: `awk '{print $2,$5}' [PANTHER DB LIST]  > [GENE-FAMILY LIST]  awk '{print $2,$3}' [PANTHER DB LIST]  > [GENE-DEFINITION LIST]`
`python [ANNOTATION SCRIPT] -i [ANNOTATED CICERO BED FILE] -p [GENE-FAMILY LIST] -c GeneFamily -o [OUTPUT ANNOTATED FILE]`  
`python [ANNOTATION SCRIPT] -i [ANNOTATED CICERO BED FILE] -p [GENE-DEFINITION LIST] -c GeneDefinition -o [OUTPUT ANNOTATED FILE]`  
## 13. Get disease relevant co-accessible sites
#### Required files:
annotated Cicero File
annotation script
Gene-Disease association table from `http://www.disgenet.org/static/disgenet_ap1/files/downloads/curated_gene_disease_associations.tsv.gz`
OR 
`http://www.disgenet.org/static/disgenet_ap1/files/downloads/befree_gene_disease_associations.tsv.gz`
OR 
`http://www.disgenet.org/static/disgenet_ap1/files/downloads/all_gene_disease_associations.tsv.gz`
OR 
`http://www.disgenet.org/static/disgenet_ap1/files/downloads/all_gene_disease_pmid_associations.tsv.gz`
#### output: 
gene-disease association list
annotated Cicero file
#### Directions: 
pare down the gene-disease list to the diseases of interest
cmd: `awk '{print $2,$6}' [DISGENET TABLE]  > [GENE-DISEASE LIST]`
`python [ANNOTATION SCRIPT] -i [ANNOTATED CICERO BED FILE] -p [GENE-DISEASE] -c DiseaseAssociate -o [OUTPUT ANNOTATED FILE]`  


## . Make swoopy figures for sites of interest, for different cell populations







## Not organized yet 
lab general python binary: `/home/data/anaconda3/bin/python`
script to try and do the stupid peak calls and get the stupid, stupid bigwig files. i really hate that I've failed to generate these like 5 times. why does this take so long. why the heck can't the authors of macs2 update their gd software to make it work on python 3 so I can run it on my own computer. it would by done in 10 minutes:
`/home/data/anaconda3/bin/python /home/joshchiou/scripts/call_cluster_peaks.py -t /home/jnewsome/pipeline_islet/SAMN1120500120-hi-cyto.filt.rmdup.tagAlign.gz -c /home/jnewsome/pipeline_islet/peakCalls/Islet.cluster_labels.txt -o SAMN1120500120-hi-cyto`  

`/home/data/anaconda3/bin/jupyter notebook`

human ref genome fasta:
`/home/joshchiou/references/male.hg19.fa`
hiv ref genome fasta
`/nfs/lab/joshchiou/HIV_pbmc_snATAC/HIV1_refseq/HIV-1.fa`



`/home/joshchiou/references/gencode.v19.cicero_gene_annotation.txt`

current snATAC lung data:
`/home/joshchiou/joshchiou-data2/lung_snATAC/merged_samples/cicero/`



remote jupyter server startup command: `jupyter notebook --no-browser --port 1113  --notebook-dir /`

remote jupyter server connection command from local machine: `ssh -p 2022 -NL 1114:localhost:1114 jnewsome@gatsby.ucsd.edu`

local jupyter server setup command: `jupyter notebook --notebook-dir /mnt/4d60fd49-d4ad-42d2-ac64-5b3f0265b9c1
`

__lists for bash:__
__macro populations__
`macroList=(CD4_T-cell Monocyte CD8_T-cell NK_cell B-cell Megakaryocyte pDC)`
__micro populations__
`microList=(Adaptive_NK_cell Activated_CD4_T-cell Activated_CD8_T-cell Conventional_Dendritic_cell Memory_B-cell \
Non-Classical_Monocyte Regulatory_CD4_T-cell Classical_Monocyte Naive_B-cell pDC Cytotoxic_NK_cell Naive_CD4_T-cell \
Megakaryocyte Naive_CD8_T-cell)`
__subclusters__
`subclusterList=(5,6 8,8 2,6 1,11 2,7 11,3 8,5 5,4 0,3 6,2 1,10 0,2 5,9 4,11 4,6 16,5 10,4 0,6 6,6 1,8 2,4 8,0 6,5 3,1 \
6,0 10,6 6,8 6,10 9,3 1,6 3,2 6,11 1,7 4,3 2,0 10,7 0,1 0,7 1,1 4,10 10,1 1,4 6,12 2,5 11,6 2,1 2,2 5,3 4,9 13,0 9,5 5,7 \
14,21 1,13 1,0 0,0 8,1 3,8 2,8 3,9 5,2 2,16 9,15 12,0 0,4 2,13 14,1 1,9 0,11 8,6 5,8 6,7 2,14 13,1 13,10 0,9 9,0 6,3 7,3 \
11,9 13,3 4,13 1,12 4,8 5,12 8,4 11,2 1,2 0,14 6,4 7,0 0,8 6,1 14,16 5,1 16,0 12,1 3,7 11,4 3,0 4,5 3,5 7,7 9,13 3,6 \
14,15 2,3 5,10 0,12 6,14 14,0 16,2 7,2 8,7 1,5 3,13 12,4 16,24 3,4 0,10 1,3 16,13 0,5 5,5 8,9 7,11 12,2 15,8 2,11 9,1 \
0,13 4,2 4,14 3,3 4,0 6,16 2,9 8,2 7,8 12,9 14,5 2,10 4,7 16,6 8,3 11,8 3,10 3,12 6,9 10,14 13,2 8,10 8,11 11,0 9,12 \
5,11 14,9 9,7 6,13 9,2 9,8 9,11 7,6 1,14 14,17 2,15 10,5 10,9 6,18 11,7 2,12 14,11 12,10 14,6 9,6 3,11 0,15 1,16 4,4 \
12,5 5,0 7,21 12,13 7,1 14,3 7,15 6,15 4,12 4,15 16,3 14,2 11,15 16,8 16,19 6,19 13,5 10,10 15,12 10,16 15,11 10,11 \
9,10 12,8 11,1 15,4 16,7 11,5 12,3 7,4 11,12 14,19 9,4 9,9 12,12 14,13 6,17 15,2 5,14 16,4 4,1 12,11 13,7 15,9 14,7 \
13,13 16,9 7,14 16,17 10,12 7,5 7,10 4,18 7,20 15,1 15,0 14,8 16,1 11,13 5,13 15,7 11,10 7,12 14,10 8,12 14,18 11,11 \
2,17 14,12 14,20 9,16 3,15 4,16 14,4 10,3 15,5 15,3 13,9 10,0 16,10 5,15 7,22 12,6 14,14 7,9 14,22 16,22 12,7 7,24 13,6 \
4,17 7,17 13,12 13,8 13,4 3,14 10,8 9,14 8,14 15,14 15,13 7,16 16,14 15,15 7,13 1,15 10,13 16,11 15,6 16,15 16,12 7,18 \
13,11 7,23 16,16 15,10 8,13 16,20 14,24 16,27 16,21 16,23 16,18 10,2 4,19 11,14 7,19 4,20 10,15 5,16 16,25 14,23 16,28 \
16,26 )`




> Written with [StackEdit](https://stackedit.io/).


