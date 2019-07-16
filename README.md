###### Gaulton Lab snATAC Pipeline Directions

# 1. Run Cell Ranger on the raw FASTQ snATAC files
(https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger)

Commands: 
	``` I don't actually know what this pipeline entails. I have never run this ```

# 2. Run the 10xGenomics processing script that Josh wrote
Commands:
```python snATAC_pipeline_10X.py -b [cell ranger output bam file] -o [output directory ? ] -n [project name prefix] -t 24 -m 2 ```
Example:
``` /home/data/anaconda3/bin/python /nfs/lab/projects/pbmc_snATAC/scripts/snATAC_pipeline_10X.py -b cellranger_v1.1/outs/possorted_bam.bam -o lab_pipeline -n pbmc1 -t 24 -m 2 ```

Required files: 
Cell ranger output bam file, example: `cellranger_v1.1/outs/possorted_bam.bam `
Josh's 10x pipeline script ('snATAC_pipeline_10X.py'), example: `/nfs/lab/projects/pbmc_snATAC/scripts/snATAC_pipeline_10X.py`

Required software: 
Python Packages:
	pysam
	numpy
	pandas
	scipy
	multiprocessing
