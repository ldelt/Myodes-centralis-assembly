# Hybrid assembly of the Myodes centralis genome

## Introduction

The lack of a suitable reference is one of the key problems faced by phylogenetic studies. Existing reference assemblies are often evolutionarily too far from the studied organism, which complicates the study of the phylogenetics of the latter. To study the hybridization of the *Arvicolinae* subfamily, the Laboratory of Evolutionary Genomics and Paleogenomics of ZIN RAN decided to assemble its own reference of *Myodes centralis*, based on Oxford Nanopore and Ilumina RNASeq sequencing data.

## Goals and objectives
**Goal: Perform hybrid de novo assembly of the Myodes centralis genome using data from Oxford Nanopore Illumina RNAseq**

Objectives:
- Carry out basecalling of raw Nanopore signal
- Assemble Nanopore reads
- Polish assembly using Nanopore reads
- Improve assembly quality using Illumina RNAseq reads
- Scaffold contigs using a closely related species as a reference
- Finalize assembly and evaluate quality

## Materials and methods

#### 1. Available data
- Raw Myodes centralis genome sequencing data in fast5 format obtained using the ONT MinION sequencer. Flow cell is FLO-MIN106D. Library kit is SQK-PSK004. 
- Illumina RNASeq sequencing results for Myodes centralis
- [*Myodes glareolus*](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_902806735.1/) genome assembly from NCBI RefSeq

#### 2. Methods
- The methods are described in detail below in the Workflow chapter. 
- A brief outline of the study is shown in the diagram

![Alt text](Experiment_scheme.png?raw=true)

## Results ans conclusion
- A hybrid assembly of Myodes centralis was carried out. 
- The total build length was 214Mb. Average coverage - 7. N50 for continues - 2382. N50 for scaffolds - 1.6Mb. 
- The polishing in the Racon didn't seem to improve the build quality at first glance. 
- The results of Nanopolish's work and their further processing will be included in the final report.
- Comparison of the assembly at various stages is presented in the table below

![Alt text](Quast_Summary_Table.png?raw=true)

## Workflow

### 1. Basecalling with Guppy
##### 1.1 Data preprocessing

- Due to the fact that basecoling is performed on Google Colab's servers, which have limitations both in terms of the amount of ROM and GPU usage time, we divided the raw fast5 sequencing results into ~7GB packages (100 fast5 files per package, taking into account their size). In total, we got 21 packages
```bash
cp *_{start..end}.fast5 ../../guppy_results/skip_XX/raw_data
```
	- XX - package number;
	- start..end - numbers of files, starting from which and ending with which we copy files to the current package.

- Next, we download the packages from the server. We upload packages to the GoogleDrive drive, since it will be mounted to the Google Colab server. Packages are loaded one at a time due to the limited amount of disk space in GoogleDrive.
```bash
scp -i BIkeys -P 8122 -r bgerda@84.204.46.23:/home/bgerda/guppy_results/skip_XX  GoogleDrive/Guppy/
```

##### 1.2 Basecalling in Google Colab

- The basecalling process is based on the [Nanopore basecalling on Google Colab](https://gist.github.com/sirselim/13f70ae69f2a512e7d9e1f00f9704f53) protocol
- Create a new notebook and change the runtime to GPU: Runtime > Change runtime type > Hardware accelerator > GPU
- Mount Google Drive with the current package. It will also store the results.
```python
from google.colab import drive
drive.mount('/content/gdrive', force_remount=True)
```

- Download the latest version of Guppy for GPU. We unpack it and check the performance. You can get a link to download Guppy on the ONT website from under a verified business account: [Guppy basecaller - Nanopore Community](https://id.customers.nanoporetech.com/app/nanoporetech-customers_myaccount_1/exk2kkmfwpBAaT3WI697/sso/saml?RelayState=https://community.nanoporetech.com/downloads)
```bash
%%shell
wget url_to_download_guppy_GPU_6.4.2
tar -xf ont-guppy_6.4.2_linux64.tar.gz
./ont-guppy/bin/guppy_basecaller --version
```

- Define the Guppy learning model used to predict outcomes. The model depends on two key parameters: the flow cell and the library preparation kit. In our case, the FLO-MIN106D cell and the SQK-PSK004 kit were used. You can find a match between the models and the cell/whale by entering the command below
```bash
%shell
guppy_basecaller --print_workflows > models_list.txt
```

- It is important to note that models_list.txt file will only contain high accuracy models (HAC). To use fast or super high accuracy (SUP) models, replace the _HAC suffix in the model name with _fast or _sup. In this work, we used SUP models.

- Start the basecalling on the current package. Depending on the provided GPU, the command execution time can vary from ~40 minutes to ~6 hours. The parameters were selected based on the  [Jetson Xavier basecalling notes](https://gist.github.com/sirselim/2ebe2807112fae93809aa18f096dbb94)
```bash
%%shell
./ont-guppy/bin/guppy_basecaller -i /content/gdrive/MyDrive/Guppy/skip_XX/raw_data -s /content/gdrive/MyDrive/Guppy/skip_XX/SUP -c dna_r9.4.1_450bps_sup.cfg --device auto --recursive --gpu_runners_per_device 16 --cpu_threads_per_caller 2 --chunks_per_runner 512 --min_qscore 10
```
	-i - path to the folder with source data in fast5 format;
	-s - folder for saving results;
	-c - the learning model defined in the previous paragraph;
	--device - used GPU-accelerator;
	--recursive - recursively search for child directories in the source data folder;
	--gpu_runnders_per_device - number of parallel neural network launches per GPU. Can significantly increase the speed of computing on GPUs with a large number of cores;
	--cpu_threads_per_caller - number of processor threads per running basecalling process. Google Colab provides only 2 threads. In addition, the main calculations take place on the GPU, which is why this parameter has little effect;
	--chunks_per_runner - the maximum number of chunks transmitted to one neural network before the start of calculations. This parameter significantly increases the speed of execution, however, it reaches a plateau after the value of 1024 and is limited by the size of the GPU memory. In our case, the optimal value was 512
	 --min_qscore - minimum qscore for filtering reads to the pass directory. Reads that do not pass the threshold will be written to the fail directory. The value 10 is chosen based on the requirements of the flye assembler.

- Download the SUP directory as a zip archive from GoogleDrive (we did this manually after each basecolling iteration).
- Send the archive with the results for the current package to the server
```bash
scp -i BIkeys -P 8122 ARCHIVE login@server:/home/bgerda/guppy_results/skip_XX
```

- Remove the current package and results from GoogleDrive. 
- Repeat the whole process with the remaining packages. If Google Colab refuses to provide a GPU - register a new account 

##### 1.3 Combining basecalling results
- Copy all the archives with the results to a separate folder and unpack
```bash
cp guppy_results/skip_*/SUP-*.zip guppy_SUP
find . -name "*.zip" -exec unzip -B -d . {} \;
```
	-d - adds a ~N suffix to filenames, where N is the number of the file with the duplicate name. Required, since Guppy always gives the same names to the results of its work.

- At the moment, we have received two folders - pass and fail, which store all the results of basecalling. However, after the previous step, filenames with initially duplicate names end in .fastq~N format. To avoid future formatting errors, we wrote a short **rename.py** script that renames files to the \*-N.fastq format. The script is available in the scripts folder of the repository
- Combine the results into the files resulted_pass.fastq and resulting_fail.fastq
```bash
cat pass/*.fastq > resulted_pass.fastq
cat fail/*.fastq > resulted_fail.fastq
```

- Checking the correctness of the combined results in [biopet-validatefastq](https://github.com/biopet/validatefastq)
```bash
biopet-validatefastq -i resulted_pass.fastq
biopet-validatefastq -i resulted_fail.fastq
```

- The validatefastq results told us that there were duplicated ids. Duplicates arose due to incorrectly specified ranges when creating individual packages. As a result, the reads from the border files ended up in two packages at the same time.
- To remove duplicates, we used the remove_dup.py python script, the code for which is taken from the corresponding topic on [StackOverflow](https://stackoverflow.com/questions/69425555/how-to-remove-duplicate-sequences-in-a-fasta-file). To run the script, you need to install the [Biopython](https://biopython.org/) library. The script is available in the scripts folder of the repository.
- Count the number of reads in the final files:
```bash
echo $(cat resulted_pass_wd.fastq|wc -l)/4|bc 
echo $(cat resulted_fail_wd.fastq|wc -l)/4|bc
```
Number of reads which passed quality filter: 6550650
Number of unpassed reads: 1486168
**Percentage of reads that passed filtering: 0.8151**

### 2. Assembling using Flye

- We used the [Flye](https://github.com/fenderglass/Flye) to assemble the genome. First, Flye is recommended to use by ONT to assemble long genomes. Secondly, Medaka is trained to work with assemblies obtained with Flye. Using Medaka with results from Canu and other assemblers may result in less accurate results.
- Only reads that passed the filtering were used for assembly, since the Flye parameters, calculated for the assembly of reads obtained using SUP models, imply a relatively high quality of reads.
- Any additional quality control was not carried out, as it is not recommended and is performed automatically by Flye.
```bash
flye --nano-hq resulted_pass_wd.fastq --out-dir ./flye_results --threads 40
```
	--nano-hq - a parameter indicating that the input data is obtained using the SUP model and have high quality
	--out-dir - folder to save results
	--threads - number of CPU threads

- Basic information about assembly results
| Total length | Contigs | N50  | Largest con. | Scaffolds | Mean coverage |
|--------------|---------|------|--------------|-----------|---------------|
| 205909836    | 83399   | 2382 | 13717        | 0         | 7             |

- The analysis of the assembly results was carried out in [quast](https://github.com/ablab/quast)
```bash
quast -t 40 -o flye_results/quast_flye flye_results/assembly.fasta
```
	-t - number of CPU threads
	-o - folder for saving results

- Results of quast are available at the quast directory at the repository 
- A comparison of the results of the analysis of the raw assembly and the assembly after all further processing steps is given in the summary table at the end of this report.

### 3. Polishing results without using Illumina reads or other sequencers

- There are several tools for polishing assembly results of ONT long reads. The main ones are: [Nanopolish](https://github.com/jts/nanopolish), [Racon](https://github.com/isovic/racon) and [Medaka](https://github.com/nanoporetech/medaka). Medaka is the most up-to-date and on average gives the best results, but it is based on deep learning algorithms and efficiently performs calculations only on the GPU. In addition, it does not cope well with a large number of short contigs and cannot be effectively  parallelized between multiple processor threads.

##### 3.1 Polishing using Racon
- For older versions of Medaka, one or 4 rounds of polishing in Racon is a mandatory step, since the models were trained to work with such data. Since we originally planned to use Medaka, Racon was the first tool we used. It is important to note that the new versions of Medaka do not require polishing in Racon and assume the use of a raw assembly from Flye.
- We used [minimap2](https://github.com/lh3/minimap2) to align the nanopore reads to the assembly.
```bash
minimap2 -a -t 40 flye_results/assembly.fasta resulted_pass_md.fastq > racon_results/minimap/algnment.sam
```

- Polighing in Racon
```
racon -m 8 -x -6 -g -8 -w 500 -t 14 resulted_pass_md.fastq racon_results/minimap/algnment.sam flye_results/assembly.fasta > racon_results/assembly_racon.fasta
```
	-m -x -g -w - parameters adjusting scores for match, missmatch, gap penalty and window size. The values of the parameters are set in accordance with the manual for working with Medaka and are determined by the fact that Medaka was trained on them.

- The polishing results are presented in a summary table. In general, the results do not seem to have improved. The number of long contigs has decreased, N50 has decreased from 2382 to 2299.

##### 3.2 Polishing using Nanopolish

- While polishing with Nanopolish gives better results than polishing with Racon alone (no further processing with Medaka), it uses raw fast5 sequencing results and is therefore significantly more time consuming. Especially for eukaryotic genomes. It also doesn't handle multiple CPU threads well.
- We used the same alignment as for Racon. However, nanopolish requires a sorted bam file with indexes built so we used [samtools](https://github.com/samtools/samtools) to build indexes. In addition, the initial alignment indexes must also be built.
```bash
bwa index flye_results/assembly.fasta
samtools view racon_results/minimap/algnment.sam -Sb | samtools sort - -@40 -o nanopolish_results/mapping.sorted.bam
samtools index pilon_results/BWA_1/mapping.sorted.bam -@ 40
```

- Nanopolish requires for its work indexes obtained on the basis of both raw signal in fast5 format and fastq basecalling results. Also, to speed up the work, it needs sequencing_summary.txt, which stores information related to the basecalling process. Since we basecalled in a non-standard way, we also need to merge all sequencing_summaries into one file before indexing.
```bash
cat SUP/sequencing_summary* > sequencing_summary.txt
nanopolish index -d ../NanoporeData/fast5_skip -s sequencing_summary.txt resulted_pass_md.fastq
```

- Start the polishing process. Command is taken directly from the manual on the nanopolish github page.
```bash
python3 nanopolish/scripts/nanopolish_makerange.py flye_results/assembly.fa | parallel --results nanopolish_results/nanopolish.results -P 5 nanopolish variants --consensus -o nanopolish_results/polished.{1}.vcf -w {1} -r resulted_pass_wd.fastq -b nanopolish_results/mapping.sorted.bam -g flye_results/assembly.fasta -t 8 --min-candidate-frequency 0.1
```

- Сombine the received .vcf files with individual corrected contigs into the final .fasta assembly file
```bash
nanopolish vcf2fasta -g flye_results/assembly.fasta nanopolish_results/polished.*.vcf > nanopolish_results/nanopolish_assembly.fasta
```

- The polishing results are presented in a summary table. The quality of the assembly has not decreased, but the quality of the consensus has probably increased.

##### 3.3 Polishing using Medaka

- Judging by the available tutorials, Medaka is currently the best tool for polishing results, excluding polishing methods with higher quality reads than ONT. In addition, its creators claim that it is about 50 times faster than Nanopolish. However, the latter applies only to work on the GPU. In our case, Medaka is extremely slow due to the following reasons: (1) the problem of the lack of a GPU, which, due to the limitations of GoogleDrive disk space, cannot be bypassed using GoogleColab; (2) a large number of short contigs with which Medaka does not work well; (3) the ability to perform calculations on several processor threads only by selecting a separate thread for a separate contig, which is unrealizable in our case.
- Before running Medaka, you need to decide on the model to use, which, like in the case of Guppy, depends on the methods for obtaining input data. You can get a list of models by typing the command below. The model name has the following format: {pore}\_{device}\_{caller variant}\_{caller version}. The closest version of Guppy available in Medaka to the one we used is 5.0.7, which is why we chose it.
```bash
medaka tools list\_models > models_medaka.txt
```

- Run Medaka	
```bash
medaka_consensus -i resulted_pass_wd.fastq -d flye_results/assembly.fasta -o guppy_SUP/results_medaka_flye -m r941_min_sup_g507 -t XXX -b XXX
```
	-i - path to file with reads;
	-d - path to raw assembly from Flye;
	-o - folder for saving results;
	-m - model.

- Unfortunately, Medaka did not manage to complete its work by the deadline for this project and therefore its results will not be included in this report.

### 3. Polishing results with Illumina RNAseq reads

- Since we had the RNAseq results, we decided to use them to polish the resulting assembly. Polishing was performed using a [Pilon](https://github.com/broadinstitute/pilon). Running Pilon on RNAseq with standard parameters could lead to disruption of the assembly structure, so we limited the correctable errors to only SNPs and short indels. 
- We polished both the raw Flye assembly and the finished Racon assembly to compare the results. It is recommended to use Pilon after Medaka, but we did not have such an opportunity
- Building the assembly index
```bash
bwa index ../results_racon_v2/assembly_racon.fasta
```

- Aligning RNAseq reads to assemblies, converting results to .bam, sorting and building indexes
```bash
bwa mem -t 40 results_racon_v2/assembly_racon.fasta ../RNAseq/RNA_merged_R1.fastq ../RNAseq/RNA_merged_R2.fastq | samtools view - -Sb | samtools sort - -@40 -o pilon_results_racon/mapping.sorted.bam
samtools index pilon_results_racon/mapping.sorted.bam -@ 40 
```

- Run Pilon
```bash
pilon --genome results_racon_v2/assembly_racon.fasta --fix bases --changes --frags pilon_results_racon/mapping.sorted.bam --output pilon_results_racon/pilon_1 -Xmx200G
```
	--genome - path to assembly; 
	--frags - path to alignment;
	--output - folder to save results;
	-Xmx200G - the parameter responsible for the RAM allocated by java;
	--fix - type of allowed corrections. The bases value limits pilon to only SNPs and short indels.

- The results were evaluated in quast and presented in a summary table 

### 4. Scaffolding using RagTag

- Scaffoldind was performed with [RagTag](https://github.com/malonge/RagTag) using the [*Myodes glareolus*](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_902806735.1/) genome assembly as the closest relative. All parameters were left as default, except for the number of processor threads.
```bash
ragtag.py scaffold GCF_902806735.1_Bank_vole1_10x_genomic.fna pilon_results/pilon_1.fasta -o ragtag_results -t 40
```

- Assembly quality assessment was carried out in Quast and presented in a summary table
