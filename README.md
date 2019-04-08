# Notes on SNP-related tools and genome variation analysis

These notes are not intended to be comprehensive. They include notes about methods, packages and tools I would like to explore. For a comprehensive overview of the subject, consider [other bioinformatics resources](https://github.com/mdozmorov/blogs/tree/master/Bioinformatics). Issues with suggestions and pull requests are welcome!


# Table of content

* [SNP tools](#snp-tools)
* [SNP callers](#snp-callers)
  * [Deep learning SNP callers](#deep-learning-snp-callers)
* [SNP annotations](#snp-annotations)
* [SNP signatures](#snp-signatures)
* [SNP pathogenicity scores](#snp-pathogenicity-scores)
* [SNP databases](#snp-databases)

## SNP tools

- `ABRA` - Assembly Based ReAligner, https://github.com/mozack/abra

- `BAMsurgeon` tools for adding mutations to existing .bam files, used for testing mutation callers, https://github.com/adamewing/bamsurgeon

- `gwasTools` - A collection of R scripts that might be useful for exploring and plotting GWAS results. https://github.com/bnwolford/gwasTools

- `maftools` - Summarize, Analyze and Visualize MAF files from TCGA or in house studies. Bioconductor, https://bioconductor.org/packages/release/bioc/vignettes/maftools/inst/doc/maftools.html, and GitHub, https://github.com/PoisonAlien/maftools

- `manhattanly` - Interactive Manhattan plots, https://cran.r-project.org/web/packages/manhattanly/

- `mutcraft` - R tools to mine & craft somatic mutations from cancer genomes, https://github.com/EmilieT/mutcraft

- `MutScan` - Detect and visualize target mutations by scanning FastQ files directly. https://github.com/OpenGene/MutScan

- `PyVCF` - A Variant Call Format Parser for Python. https://pyvcf.readthedocs.io/en/latest/. Has `vcf_melt` tool to reformat a VCF into long format.

- `samplot` - Plot structural variant signals from many BAMs and CRAMs. https://github.com/ryanlayer/samplot

- `ttplot` - Tao Yan's Plot Toolkit, plots LD Heatmap, Manhattan plot. https://github.com/YTLogos/ttplot

- `vcfR` - Manipulate and Visualize VCF Data. https://cran.r-project.org/web/packages/vcfR/index.html


## SNP callers

- Xu, Chang. “A Review of Somatic Single Nucleotide Variant Calling Algorithms for Next-Generation Sequencing Data.” Computational and Structural Biotechnology Journal 16 (2018): 15–24. https://doi.org/10.1016/j.csbj.2018.01.003. - Overview of 46 somatic Single Nucleotide Variant (SNV) caller tools. Pre-processing, variant evaluation, and post-filtering steps. Four categories of algorithms, description of each, and the corresponding tools: matched tumor-normal (position-, haplotype-, mathine learning-based methods, [Table 1](https://www.sciencedirect.com/science/article/pii/S2001037017300946?via%3Dihub#t0005)), single-sample ([Table 2](https://www.sciencedirect.com/science/article/pii/S2001037017300946?via%3Dihub#t0010), some offer somatic-germline classification), UMI-based (UMI technology, [Figure 1](https://www.sciencedirect.com/science/article/pii/S2001037017300946?via%3Dihub#f0005), [Table 3](https://www.sciencedirect.com/science/article/pii/S2001037017300946?via%3Dihub#t0015)), and RNA-seq (Technology, issues, [Table 4](https://www.sciencedirect.com/science/article/pii/S2001037017300946?via%3Dihub#t0020)) variant calling. Benchmarking using tools for generating synthetic reads, spike-ins, GiAB, melanoma-normal samples, performance evaluation metrics. Issues in representing complex variants and tools for variant normalization. Deep neural network-based algorithms perform best.

- `MutationSeq` - somatic SNV detection from tumor-normal pairs. http://compbio.bccrc.ca/software/mutationSeq/




### Deep learning SNP callers

- `NeuSomatic` - Convolutional neural network (9 layers) for somatic mutations calling. Reads within 7bp window around candidate SNP are extracted, realigned, summarized into matrices for tumor-normal samples,  used for classifying mutation type, length, position. Tested on GiB samples and on DREAM datasets. Comparison with other SNP callers. https://github.com/bioinform/neusomatic
    - Sahraeian, Sayed Mohammad Ebrahim, Ruolin Liu, Bayo Lau, Marghoob Mohiyuddin, and Hugo Y. K. Lam. “Deep Convolutional Neural Networks for Accurate Somatic Mutation Detection,” September 4, 2018. https://doi.org/10.1101/393801.


## SNP annotations

- `slivar` - variant expressions, annotation, and filtering, by Brent Pedersen, https://github.com/brentp/slivar


## SNP signatures

- https://github.com/ictic-bioinformatics/CANCERSIGN - CANCERSIGN: a user-friendly and robust tool for identification and classification of mutational signatures and patterns in cancer genomes. Masroor Bayati, Hamid Reza Rabiee, Mehrdad Mehrbod, Fatemeh Vafaee, Diako Ebrahimi, Alistair Forrest, Hamid Alinejad-Rokny. bioRxiv 424960; doi: https://doi.org/10.1101/424960

- https://github.com/danro9685/SparseSignatures, https://bioconductor.org/packages/release/bioc/html/SparseSignatures.html - De Novo Mutational Signature Discovery in Tumor Genomes using SparseSignatures. Daniele Ramazzotti, Avantika Lal, Keli Liu, Robert Tibshirani, Arend Sidow. bioRxiv 384834; doi: https://doi.org/10.1101/384834

- https://bioconductor.org/packages/release/bioc/html/YAPSA.html - Yet Another Package for Signature Analysis, functionality used in L. Alexandrov et al., Nature 2013


## SNP pathogenicity scores

- ClinPred - pathogenicity prediction for all nonsynonymous SNPs. Trained on ClinVar, validated on nine other databases. Random forest and gradient boosted decision tree, comparison with other machine learning algorithms. Downloadable scores for all nonsynonymous SNPs, https://sites.google.com/site/clinpred/home
    - Alirezaie, Najmeh, Kristin D. Kernohan, Taila Hartley, Jacek Majewski, and Toby Dylan Hocking. “ClinPred: Prediction Tool to Identify Disease-Relevant Nonsynonymous Single-Nucleotide Variants.” The American Journal of Human Genetics 103, no. 4 (October 2018): 474–83. https://doi.org/10.1016/j.ajhg.2018.08.005.

## SNP databases

- `clinvar` - This repo provides tools to convert ClinVar data into a tab-delimited flat file, and also provides that resulting tab-delimited flat file. https://github.com/macarthur-lab/clinvar

- Clinical Interpretation of Variants in Cancer database, http://www.civicdb.org/. CIViC interface public API, http://griffithlab.org/civic-api-docs/

