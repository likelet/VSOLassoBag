## Parameters 
> Those parameters would cover the setting from `nextflow.config` file .  

- [Parameters](#parameters)
    + [Mandatory](#mandatory)
    + [Configuration](#configuration)
    + [Optional](#optional)
    + [Detailed instruction of parameters](#detailed-instruction-of-parameters)
- [Config](#config)

#### Mandatory   

> (plz configure those options in *nextflow.config* or *docker.config* file) .  


| Name | Example/Default value | Description |
|-----------|--------------:|-------------|
|--reads | `./*{1,2}.fq.gz` | input raw paired reads |
|--designfile     | `FALSE` | a txt file that stored experimental design information, plz see details from `--designfile` section above |
|--comparefile     | `FALSE` | a txt file that stored experimental compare information, plz see details from `--comparefile` section above | . 

#### Configuration   
> (paths to references, softwares and special environments. Only need to be set at the first time you run the pipeline) . 

| Name | Example/Default value | Description |
|-----------|--------------|-------------|
|--redir | `/path/to/refdir` | the folder containing all reference files and index |
|--genomefile  | `path/to/refdir/chr2.fa` | Path to Fasta reference(required if not set in config file)|
|--gtffile/--bedfile/--annotationfile  | `path/to/refdir/hg38_gencode.txt` |Different annotation files from GENCODE database for annotating circRNAs. e.g. [gencode.v25.annotation.gtf]/[gencode.v25.annotation.bed]/[hg38_gencode.txt]|
|--ciridir/--find_circdir/--mapsdir/--knifedir  | `/path/to/tools/directory` | Home folder of ciri/find_circ/mapsplice/knife installed location | . 

#### Optional  

 
| Name | Default value | Description |
|-----------|--------------|-------------|
|--starindex/--bowtie2index/--bwaindex/--segindex/--bowtieindex/--refmapsplice/--knifeindex  | `-`| Path to STAR/bowtie2/segemehl/bowtie/bwa/mapsplice/knife index. If not set, the pipeline will create the index itself.  |
|--singleEnd  | `false` | specify that the reads are single ended  |
|--merge | `true` | merge the different matrixes produce by different tools and draw the venn graph|
|--separate | `false` | annotate the results separately|
|--selectTools     | `1` | specify which tools should be use. `1` for circexplorer2, `2` for ciri, `3` for find_circ, `4` for mapsplice, `5` for segemehl, `6` for knife. For example, you can set `1,2,3,4,5` for running five tools in the same time. |
|--outdir |  `./Result` | the output directory of the results|
|--mRNA |  `path/to/gencode.rsem.fpkm_m6Astatus_11_29.mat` | Path to the mRNA expression matrix. Only need to be set when you want to do the correlation.|

#### Detailed instruction of parameters    


* `--reads`  
    
    suffix of your raw reads file. For example, `*_{1,2}.fq.gz` for paired end reads file `sampleA_1.fq.gz` and `sampleA_2.fq.gz `  
    
* `--designfile`  
    
    design file  
    
* `--comparefile`  
    
    compare file 


## Config
As a nextflow-based analysis pipeline, CircPipe allow users edit configure file `nextflow.config` to set the index files and default file path parameters instead of typing them into the command line.

To configure, please go to `params` line, and set the following information of various file locations and system environment settings

[configFile](https://github.com/likelet/RNAseqPipe/blob/master/nextflow.config)