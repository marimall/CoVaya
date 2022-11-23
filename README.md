# CoVaya
Do you feel that a specific short read company's SARS Cov 2 aligner sucks?Well, welcome to the club! This python code alignes your fastqs to the SARS Cov 2 genome and returns several statistics with a number of very well knonw freely available tools. 

In order to use it please download and transfer on the same directory freebayes: https://github.com/freebayes/freebayes

## Dependencies

This pipeline depends on the hard work of several teams that created truly amazing and helpful tools free of depending environments. This tools are 

  fastp
  bwa-mem
  samtools
  freebayes
  
  
 For freebayes please download and transfer on the same directory: https://github.com/freebayes/freebayes
 For the rest a simple sudo apt install {tool} command is more than enough

## Usage

CoVaya can process multiple fastq files contained in a directory. The simple command is

    python CoVaya --input_path {DIR_NAME} --primers V4_primer_set.bed --reference reference/GCA_009858895.3.fasta
    
    
You can always specify youw own reference and primer set as long as they are in the same format of the default files

The code has two more options
1. For variant calling add --vcf True to the command:
  
    python CoVaya --input_path {DIR_NAME} --primers V4_primer_set.bed --reference reference/GCA_009858895.3.fasta --vcf True
  
2. To calculate which site are heterozygous add both --vcf True and --hetero_report True:
  
    python CoVaya --input_path {DIR_NAME} --primers V4_primer_set.bed --reference reference/GCA_009858895.3.fasta --vcf True --hetero_report True

Bonus:
  Depending on your machine you can speed up the process by selecting --pools which defines the number of subprocesses the program performs. 

## Results

This pipeline creates a result directory in your given directory. Inside you will find a seperate directory for each sample in which you will find the cleaned fastqs, the aligned and sorted bam, statistics containing read lenghth, distributions, number of aligned reads etc. More over you will find a pdf file where the depth per site per site is illustrated. Selecting the --vcf command will result in producing a .vcf for each sample, while selecting the --hetero_report will result in the creation of a heterozygous_stats directory where you can find text files per sample in which the heterozygous site are reported. Result interpatation is up to you!  

