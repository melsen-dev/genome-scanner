# genome-scanner

`genome-scanner` is an open source application to scan for Single Nucleotide Polymorphisms (SNPs) from Variant Call File (VCF) files that result from Whole Genome Sequencing (WGS) data analysis (e.g. sequencing.com, dantelabs.com).
It scans the user-provided VCF file (can be downloaded from sequencing provider) for annotated genes from the SNP database and outputs a report table with associated risk or each found SNP.
This application is an implementation of the findings in referenced papers for research purposes only. It is not recommended as guidance in medical diagnosis and treatment.

SNP Database\
[<img src="./doc/snpdb.png" width="800"/>](./doc/snpdb.png)

Report based on VCF file\
[<img src="./doc/results.png" width="1000"/>](./doc/results.png)

# VCF File Compatibility
This application has been tested with VCF v4.2 files, see https://samtools.github.io/hts-specs/VCFv4.2.pdf for details on file format. You can find example VCF files for individual genomes from the `International Genome Sample Resource (IGSR)` and the `1000 Genomes Project` at https://www.internationalgenome.org/faq/can-i-get-phased-genotypes-and-haplotypes-for-the-individual-genomes/.

# Python Code
You can clone this repository to execute the application on your local computer with the following
[Requirements](./code/requirements.txt).
You can edit the
[SNP database](./db/snpdb.csv)
to add more annotated genes.
To execute this code, run [genome-scanner](./code/genome-scanner.py).

# Web Service
Add hosted web service address
