# Neoantigen-Discovery-in-the-HCC1395-Breast-Cancer-Cell-Line
Pipeline for integrated whole-genome and transcriptome analysis to identify neoantigens in the HCC1395 breast cancer cell line. Combines somatic variant calling, RNA expression quantification, HLA typing, and neoantigen prediction for immunotherapy research.
Comprehensive Neo-antigen Discovery Pipeline
Project: HCC1395 / HCC1395BL multi-omics proteogenomics
Author: Ahmed Hassan · Last updated: 2025-08-02

1. Project Rationale
Why HCC1395 / HCC1395BL?	Details
Clinically relevant model	Triple-negative breast-cancer (TNBC) cell line with a high somatic SNV/indel burden (~8 mut/Mb) and a clonal TP53 driver – ideal for neo-antigen discovery.
Ground-truth benchmarking	Matched normal cell line (HCC1395BL) plus gold-standard variant calls released by SEQC2 / GIAB allow objective evaluation of pipeline accuracy.

2. Objectives
Call and benchmark somatic variants from matched tumour/normal WGS data.

Translate variants to candidate peptides and prioritise those that are expressed, processed, and presented on MHC-I.

Validate neo-antigens with mass-spectrometry evidence and 3-D structural modelling.

Produce publication-ready figures and a reproducible, containerised workflow suitable for reuse or extension.

3. Data Overview
Modality	Source	File examples
WGS (2 × 101 bp) | Normal – HCC1395BL | ERR194147 | ENA | ERR194147_1.fastq.gz, ERR194147_2.fastq.gz 
EMBL-EBI
WGS (2 × 101 bp) | Tumour – HCC1395 | ERR194146 | ENA | ERR194146_1.fastq.gz, ERR194146_2.fastq.gz
RNA-seq  | Tumour  |SRX8401273 - SRX8401274 – SRX8401275 |  SRA (BioProject PRJNA635123) 
RNA-seq  | Normal  |SRX5908263 - SRX5908260 – SRX5908218 |  SRA (BioProject PRJNA504037)
somatic hg38 af only gnomad genomics public data
Mills and 1000g gold standard incdels hg38 genomics public data
PON hg38 gs://gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz



4. Workflow Summary
Biological question	Workflow step (tool → container)	Key output
Which bases differ?	nf-core/sarek (BWA-MEM → Mutect2, Strelka2, …)	somatic.filtered.vcf.gz
Does DNA → protein change?	pVACseq generate_protein_fasta	neo.fasta (ref + alt peptides)
Is the mutation used?	RNA filter (--exp-rna), MS PSM filter (MaxQuant / Comet)	filtered_neo.tsv
Will HLA present it?	NetMHCpan 4.2 wrapped by pVACtools	binding_predictions.tsv


A schematic of the full pipeline is provided in docs/workflow.png.

5. Software Stack
Component	Version	Container	Citation
nf-core/sarek	-r 3.14.0	ghcr.io/nfcore/sarek:3.3	10.1038/sdata.2020.28
pVACtools	4.0.15	griffithlab/pvactools:4.0.15	10.1186/s12859-018-2693-2
NetMHCpan	4.2	local install (academic)	10.1007/978-1-4939-9173-0_4


All containers are pulled automatically by Nextflow; no manual installation required.

6. Quick Start
# Clone repository
git clone [https://github.com/ahmedhassan/proteogenomic-neoantigen.git](https://github.com/AhmedHassan-bioinfo/Neoantigen-Discovery-in-the-HCC1395-Breast-Cancer-Cell-Line/edit/main/README.md)

# Launch the workflow (Nextflow ≥23.04)
nextflow run main.nf \
  --tumour_fq data/WGS/HCC1395_*.fastq.gz \
  --normal_fq data/WGS/HCC1395BL_*.fastq.gz \
  --rna_fq    data/RNAseq/*.fastq.gz \
  --ms_raw    data/Proteomics/*.raw \
  --outdir    results \
  -profile docker
Hardware note: end-to-end run (WGS ≈ 60×, RNA ≈ 50 M reads, 4 RAW files) finishes in ~ 18 h on a 32-core machine with 32 GB RAM and < 1 TB scratch space.

7. Output Structure
Copy code
results/
├── 01_variants/
│   └── somatic.filtered.vcf.gz
├── 02_pvacseq/
│   ├── neo.fasta
│   └── binding_predictions.tsv
├── 03_proteomics/
│   └── validated_neo.tsv
├── 04_structures/
│   └── *.pdb
└── 05_figures/
    ├── volcano_affinity_vs_RNA.png
    ├── mutation_class_bar.png
    ├── peptide_length_hist.png
    └── HLA_heatmap.png
8. Interpreting the Plots
Plot	Biological insight	Actionable take-away
Volcano (affinity vs RNA)	Peptides that are both expressed and strong binders	Candidates in top-right = vaccine leads
Mutation-class bar	Frameshift vs missense contribution	Frameshifts often yield more unique epitopes
Length histogram	9-mer optimum distribution	12–13 mers → consider class II or artefacts
HLA heat-map	Allele coverage	Focus validation on alleles with sparse hits
AlphaFold cartoon	3-D plausibility of pMHC	Disordered peptide tails suggest trimming issues

9. Results in Brief
Total somatic variants: nnn

Neo-antigen candidates (NetMHCpan ≤ 500 nM & RNA TPM > 1): nn


Top vaccine lead: TP53 p.R248Q → LLGRNSFEVR (IC50 = 42 nM, TPM = 38, PSM = 3)

Full result tables are in results/03_proteomics/validated_neo.tsv.

10. How to Cite
Ahmed Hassan. Comprehensive Neo-antigen Discovery Pipeline for HCC1395/HCC1395BL. GitHub repository, 2025. https://github.com/AhmedHassan/proteogenomic-neoantigen

11. Acknowledgements
SEQC2 / GIAB consortium for truth-set variant files.

CPTAC for public proteomics data.

nf-core community for the Sarek pipeline.

12. License
This project is released under the MIT License. See LICENSE for details.

Contact
For questions, open an issue or email ahmednasser1378@gmail.com
