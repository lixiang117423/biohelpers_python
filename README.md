Some useful script tools written during data processing.

## Installation

```bash
git clone git@github.com:lixiang117423/biohelpers_python.git

cd biohelpers_python

pip install -e .
```

or:

```bash
pip install biohelpers
```

## Dependencies

### Install miniforge

Install [miniforge](https://github.com/conda-forge/miniforge) according to the instructions on the website.

### Install dependencies

```bash
mamba install conda-forge::biopython=1.85
mamba install bioconda::gffread=0.12.7
mamba install bioconda::seqkit=2.10.0
mamba install bioconda::bcftools=1.21

```

## Usage

### example data

The demo data were downloaded from [RiceSuperPIRdb](http://www.ricesuperpir.com/web/download) from the paper, [A super pan-genomic landscape of rice](https://www.nature.com/articles/s41422-022-00685-z).

```bash
wget http://www.ricesuperpir.com/uploads/common/gene_annotation/NIP-T2T.gff3.gz
wget http://www.ricesuperpir.com/uploads/common/genome_sequence/NIP-T2T.fa.gz

gunzip NIP-T2T.gff3.gz
gunzip NIP-T2T.fa.gz

mv NIP-T2T.gff3 Nipponbare.gff3
mv NIP-T2T.fa Nipponbare.fa
```

### get the meta information of a project or run from ENA

```bash
get_fq_meta -h

usage: get_fq_meta [-h] [-id ACCESSION] [-o OUTPUT] [-s SAVE [SAVE ...]]

Download sequencing metadata TSV from ENA API

options:
  -h, --help            show this help message and exit
  -id ACCESSION, --accession ACCESSION
                        ENA accession number (required) (e.g. PRJNA123456)
  -o OUTPUT, --output OUTPUT
                        Output path (supports .tsv/.csv/.txt/.xlsx extensions, default: ./tmp/[accession].meta.tsv)
  -s SAVE [SAVE ...], --save SAVE [SAVE ...]
                        Fields to save (all|field1 field2), available fields: secondary_study_accession,sample_accession,secondary_sample_acces  
                        sion,experiment_accession,study_accession,submission_accession,tax_id,scientific_name,instrument_model,nominal_length,l  
                        ibrary_layout,library_source,library_selection,base_count,first_public,last_updated,study_title,experiment_alias,run_al  
                        ias,fastq_bytes,fastq_md5,fastq_ftp,fastq_aspera,fastq_galaxy,submitted_bytes,submitted_md5,submitted_ftp,submitted_gal  
                        axy,submitted_format,sra_bytes,sra_md5,sra_ftp,sample_alias,broker_name,sample_title,nominal_sdev,bam_ftp,bam_bytes 
```

```bash
get_fq_meta -id PRJNA510920 -o PRJNA510920.meta.txt
# Meta information file saved to: PRJNA510920.meta.txt
```

### download FASTQ format data from ENA

```bash
get_fq_file -h

usage: get_fq_file [-h] --accession ACCESSION --type {ftp,aspera} [--key KEY] [--method {run,save}] [--output OUTPUT]

Download FASTQ files from ENA

options:
  -h, --help            show this help message and exit
  --accession ACCESSION, -id ACCESSION
                        Accession number (required) Format example: PRJNA661210/SRP000123 Supports ENA/NCBI standard accession formats
                        (default: None)
  --type {ftp,aspera}, -t {ftp,aspera}
                        Download protocol type ftp: Standard FTP download aspera: High-speed transfer protocol (requires private key) (default:  
                        None)
  --key KEY, -k KEY     Path to aspera private key Required when using aspera protocol Default location:
                        ~/.aspera/connect/etc/asperaweb_id_dsa.openssh (default: None)
  --method {run,save}, -m {run,save}
                        Execution mode run: Execute download commands directly save: Generate download script (default) (default: save)
  --output OUTPUT, -o OUTPUT
                        Output directory Default format: [accession].fastq.download Auto-create missing directories (default: None)
```

`-m run` will directly download the FASTQ files.  **But we strongly recommend using ` -m save` to save download script then get the FASTQ data.**

#### wget

```bash
get_fq_file -id PRJNA510920 -m save -t ftp -o ./fastq

# Please run the next command to download the FASTQ data:
# bash download_PRJNA510920_fastq_by_wget.sh
```

#### aspera

```bash
get_fq_file -id PRJNA510920 -m save -t aspera -k ./asperaweb_id_dsa.openssh  -o ./fastq

# Please run the next command to download the FASTQ data:
# bash download_PRJNA510920_fastq_by_aspera.sh
```

### download the HMM file

```bash
download_hmm -h
usage: download_hmm [-h] -id HMM_ID -o OUTPUT

Download HMM profile from InterPro

options:
  -h, --help            show this help message and exit
  -id, --hmm_id HMM_ID  Pfam HMM ID (e.g. PF00010)
  -o, --output OUTPUT   Output directory path

```

```bash
download_hmm -id PF00010 -o example/
```

### get the longest transcript for each gene

```bash
parse_longest_mrna -h
usage: parse_longest_mrna [-h] -g GENOME -f GFF3 -o OUTPUT

Extract longest mRNA transcripts

options:
  -h, --help           show this help message and exit
  -g, --genome GENOME  Input genome FASTA file
  -f, --gff3 GFF3      Input GFF3 annotation file
  -o, --output OUTPUT  Output FASTA file
```

```bash
parse_longest_mrna -g example/Nipponbare.fa -f example/Nipponbare.gff3 -o test/longest.pep.fa
```

```bash
################################################################
Total genes: 57359
Total transcripts: 67818
Genes with multiple transcripts: 6510
################################################################
Successfully extracted 57359 longest transcripts
Longest transcript protein sequences saved to: test/longest.pep.fa
Gene and transcript information saved to: example/Nipponbare.gene.info.txt
```

### process blast results

```bash
process_blast -h
usage: process_blast [-h] -i INPUT [-e EVALUE] -o OUTPUT [-n NUMBER]

Process BLAST results and filter by E-value

options:
  -h, --help           show this help message and exit
  -i, --input INPUT    Input BLAST result file
  -e, --evalue EVALUE  E-value threshold (default: 1e-5)
  -o, --output OUTPUT  Output file path
  -n, --number NUMBER  Number of top hits to retain per query (default: 1)
```

```bash
process_blast -i example/diamond.blast.txt -e 1e-6 -o test/filtered.blast.txt
```

```bash
Successfully processed 85500 query sequences, retained 85500 records.
Results saved to: test/filtered.blast.txt
```

### run [hisat2](https://github.com/DaehwanKimLab/hisat2)

```bash
usage: run_hisat2.py [-h] -x INDEX [-t THREADS] -f FOLDER [-m {run,save}] -o OUTPUT

HISAT2 RNA-seq alignment pipeline

options:
  -h, --help            show this help message and exit
  -x, --index INDEX     Reference genome index path
  -t, --threads THREADS
                        Number of threads (default: all cores)
  -f, --folder FOLDER   Input directory containing FASTQ files
  -m, --method {run,save}
                        Execution method: run immediately or save to script
  -o, --output OUTPUT   Output directory for BAM files
```

```bash
run_hisat2 -x 03.genome/acuce.genome.hisat2.index -t 60 -f 01.data -m save -o 04.mapping 
```

```bash
Please run the following command to execute the alignment:
bash run_hisat2.sh
```

### get haplotype information

```bash
tabix example/chr1.36545388.snp.vcf
```

```bash 
get_hap -h
usage: get_hap.py [-h] -v VCF -c CHR -p POSITION [-s START] [-e END] -o OUTPUT

Extract haplotype information from VCF files

options:
  -h, --help            show this help message and exit
  -v VCF, --vcf VCF     Input VCF file path
  -c CHR, --chr CHR     Chromosome identifier
  -p POSITION, --position POSITION
                        Target SNP position
  -s START, --start START
                        Upstream window size
  -e END, --end END     Downstream window size
  -o OUTPUT, --output OUTPUT
                        Output file path
```

```bash
get_hap -v example/chr1.36545388.snp.vcf -c Chr1 -p 36545388 -o test.vcf.txt
```

```bash
Chr	Position	REF	ALT	Sample	GT	Alleles	Frequency	Biological_Meaning
Chr1	36545388	C	T	100	./.	./.	85.86%	Missing
Chr1	36545388	C	T	101	0/1	C/T	8.08%	Heterozygous
Chr1	36545388	C	T	10	./.	./.	85.86%	Missing
Chr1	36545388	C	T	11	./.	./.	85.86%	Missing
Chr1	36545388	C	T	12	./.	./.	85.86%	Missing
Chr1	36545388	C	T	13	0/1	C/T	8.08%	Heterozygous
Chr1	36545388	C	T	14	./.	./.	85.86%	Missing
Chr1	36545388	C	T	15	0/1	C/T	8.08%	Heterozygous
Chr1	36545388	C	T	16	./.	./.	85.86%	Missing
```

### get_gene_pairs

Parse gene pairs like NLR-pairs from gff file and gene id file.

```bash
get_gene_pairs -h
usage: get_gene_pairs [-h] --gff GFF --id ID [--type {gene,mrna}] [--distance DISTANCE] --output OUTPUT

Find gene pairs near target genes

options:
  -h, --help            show this help message and exit
  --gff GFF, -g GFF     Path to GFF file
  --id ID, -i ID        File containing target gene IDs
  --type {gene,mrna}, -t {gene,mrna}
  --distance DISTANCE, -d DISTANCE
                        The number of other genes between pairs of genes. The default value is 3.
  --output OUTPUT, -o OUTPUT
                        Output filename
```

```bash
get_gene_pairs -g data/gff3/534M.gff3 -i result/03.nlr-pairs/534M.nlr.id.txt -t mrna -o result/03.nlr-pairs/534M.NLR-pairs.txt
```

```bash
chr	gene_id_1	gene_start_1	gene_end_1	strand_1	gene_id_2	gene_start_2	gene_end_2	strand_2
chr1	Gla4_010.100	708894	713156	+	Gla4_010.97	683583	686891	-
chr1	Gla4_010.3267	31640732	31642210	+	Gla4_010.3270	31655373	31658766	+
chr1	Gla4_010.3270	31655373	31658766	+	Gla4_010.3268	31648144	31649616	+
chr1	Gla4_010.3270	31655373	31658766	+	Gla4_010.3269	31651858	31653312	+
chr1	Gla4_010.3270	31655373	31658766	+	Gla4_010.3273	31672690	31674213	+
chr1	Gla4_010.3271	31661494	31663020	+	Gla4_010.3268	31648144	31649616	+
chr1	Gla4_010.3271	31661494	31663020	+	Gla4_010.3269	31651858	31653312	+
chr1	Gla4_010.3271	31661494	31663020	+	Gla4_010.3270	31655373	31658766	+
chr1	Gla4_010.3271	31661494	31663020	+	Gla4_010.3272	31668519	31670025	+
chr1	Gla4_010.3271	31661494	31663020	+	Gla4_010.3273	31672690	31674213	+
chr1	Gla4_010.3271	31661494	31663020	+	Gla4_010.3274	31675773	31676637	+
```
### new gff file from BRAKER result

```bash
usage: new_gff_braker.py [-h] -i INOUT -s SPECIES [-d DISTANCE] -o OUTPUT

Process Braker GTF to GFF3 with customized gene/mRNA/feature annotation.

options:
  -h, --help            show this help message and exit
  -i, --inout INOUT     Input GTF file (default: None)
  -s, --species SPECIES
                        Species name string (for prefix, e.g. Os) (default: Os)
  -d, --distance DISTANCE
                        Gene id multiplier/distance, like Os01g000010 and Os01g000020. (default: 10)
  -o, --output OUTPUT   Output GFF3 file (default: None)
```

```bash
new_gff_braker -i example/braker.gtf -s Os -d 10 -o example/braker.gff3
```

### VCF file information

```bash
stat_vcf -h
```

```bash
usage: stat_vcf [-h] [-o OUTPUT] [-v] [--sample-details] [--total-only] vcf_file

Count SNPs and INDELs per chromosome in VCF files (total and per-sample)

positional arguments:
  vcf_file              Input VCF file path

options:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        Output file path (optional, prints to stdout if not specified)
  -v, --verbose         Enable verbose output
  --sample-details      Show detailed per-sample statistics for each chromosome
  --total-only          Show only total statistics, skip per-sample analysis

Examples:
  test.py input.vcf
  test.py sample.vcf -o output.txt
  test.py variants.vcf --output results.txt --sample-details

Notes:
  - SNP: Single nucleotide polymorphism (REF and ALT both length 1)
  - INDEL: Insertion/deletion variant (REF and ALT different lengths)
  - MNV: Multi-nucleotide variant (REF and ALT same length but >1)
  - Per-sample analysis requires sample columns in VCF (columns 10+)
  - Only variants present in each sample (non-reference genotypes) are counted
  - Compressed VCF files (.vcf.gz) need to be decompressed first
```

```bash
stat_vcf input.vcf --sample-details -o clean_data.tsv
```

### run fastp

```bash
# Basic usage
run_fastp -i 01.data -o 02.clean_data

# View all parameters
run_fastp --help

# Custom parameters
run_fastp \
    -i 01.data \
    -o 02.clean_data \
    -t 8 \
    -q 25 \
    -l 100 \
    --verbose \
    --log

# Different file patterns
run_fastp \
    -i raw_data \
    -o clean_data \
    -p "*.R1.fastq.gz" \
    -t 12
```

### run RNA-Seq using HISAT2

```bash
run_rnaseq -h   

usage: run_rnaseq [-h] -g GENOME -f GTF -i INPUT -o OUTPUT [-p PATTERN] [-r {yes,y,no,n}] [-t THREADS]

RNA-seq analysis pipeline: HISAT2 + StringTie

options:
  -h, --help            show this help message and exit
  -g, --genome GENOME   Genome fasta file path
  -f, --gtf GTF         Gene annotation GTF file path
  -i, --input INPUT     Input fastq file directory or sample information file
  -o, --output OUTPUT   Output directory
  -p, --pattern PATTERN
                        Fastq file naming pattern, e.g., "*.R1.fastq.gz" or "*_1.fq.gz", * represents sample name
  -r, --remove {yes,y,no,n}
                        Remove BAM files after processing (default: no)
  -t, --threads THREADS
                        Number of threads (default: 8)
```

```bash
run_rnaseq -g 03.genome/genome.fa -f 03.genome/genome.gtf -i 01.data/raw -o output -t 60 -p "*_1.fq.gz"
```

## Requirements

- Python 3.7+
- requests>=2.31.0

