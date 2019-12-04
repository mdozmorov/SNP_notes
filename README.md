# Notes on SNP-related tools and genome variation analysis

These notes are not intended to be comprehensive. They include notes about methods, packages and tools I would like to explore. For a comprehensive overview of the subject, consider [other bioinformatics resources](https://github.com/mdozmorov/Bioinformatics_notes) and [collections of links to various resources](https://github.com/mdozmorov/MDmisc_notes). Issues with suggestions and pull requests are welcome!


# Table of content

* [Variant calling pipelines](#variant-calling-pipelines)
  * [Preprocessing tools](#preprocessing-tools)
  * [Depth](#depth)
* [SNP callers](#snp-callers)
  * [Deep learning SNP callers](#deep-learning-snp-callers)
* [SNP annotations](#snp-annotations)
* [SNP signatures](#snp-signatures)
* [SNP pathogenicity scores](#snp-pathogenicity-scores)
* [SNP visualization](#snp-visualization)
* [SNP databases](#snp-databases)
* [InDels](#indels)
* [CNV, SV](#cnv--sv)
* [Miscellaneous](#miscellaneous)

## Variant calling pipelines

- `MapCaller` - map reads and identify genomic variants, indels, inversions, translocations. Based on KART method. Outperforms GATK, Freebayes, Mpileup. https://github.com/hsinnan75/MapCaller
    - Lin, Hsin-Nan, and Wen-Lian Hsu. “MapCaller - An Integrated and Efficient Tool for Short-Read Mapping and Variant Calling Using High-Throughput Sequenced Data.” Preprint. Bioinformatics, September 26, 2019. https://doi.org/10.1101/783605.

- TOPMed Variant Calling Pipeline, https://github.com/statgen/topmed_variant_calling

- `DNAscan` - a pipeline for DNA-seq analysis to call SNPs, indels, SVs, repeat expansions and viral genetic material. Four parts - alignment (HISAT2, BWA mem), analysis (Freebayes, GATC HC, Manta, Expansion Hunter), annotation (Annovar), report generation. https://github.com/KHP-Informatics/DNAscan
    - Iacoangeli, A., A. Al Khleifat, W. Sproviero, A. Shatunov, A. R. Jones, S. L. Morgan, A. Pittman, R. J. Dobson, S. J. Newhouse, and A. Al-Chalabi. “DNAscan: Personal Computer Compatible NGS Analysis, Annotation and Visualisation.” BMC Bioinformatics 20, no. 1 (December 2019): 213. https://doi.org/10.1186/s12859-019-2791-8.

### Preprocessing tools

- `ABRA` - Assembly Based ReAligner, https://github.com/mozack/abra

- `alleleCount` - Takes a file of locations and a [cr|b]am file and generates a count of coverage of each allele [ACGT] at that location. https://github.com/cancerit/alleleCount

- `BAMsurgeon` tools for adding mutations to existing .bam files, used for testing mutation callers, https://github.com/adamewing/bamsurgeon

- `PyVCF` - A Variant Call Format Parser for Python. https://pyvcf.readthedocs.io/en/latest/. Has `vcf_melt` tool to reformat a VCF into long format.

- `SURVIVOR` - Toolset for SV simulation, comparison and filtering. https://github.com/fritzsedlazeck/SURVIVOR/tree/1.0.7

- `VariantQC` - VCF quality control tool, part of DISCVRseq toolkit. Uses GATK4 engine. Java wrapper of GATK$'s VariantEval  tool. Input - VCF file and an indexed genome FASTA file. Output - MultiQC-templated report. https://github.com/BimberLab/DISCVRSeq/
    - Yan, Melissa Y, Betsy Ferguson, and Benjamin N Bimber. “VariantQC: A Visual Quality Control Report for Variant Evaluation.” Edited by Jonathan Wren. Bioinformatics, July 16, 2019, btz560. https://doi.org/10.1093/bioinformatics/btz560.

- `vcfR` - Manipulate and Visualize VCF Data. https://cran.r-project.org/web/packages/vcfR/index.html

### Depth

- `BAMscale` - BAMscale is a one-step tool for either 1) quantifying and normalizing the coverage of peaks or 2) generated scaled BigWig files for easy visualization of commonly used DNA-seq capture based methods.
- `mosdepth` - fast BAM/CRAM depth calculation for WGS, exome, or targetted sequencing., https://github.com/brentp/mosdepth
- `indexcov` - fast genome coverage, aberrant coverage detection, infer sex. Visualization. https://github.com/brentp/goleft
- `bamCoverage` - BAM to bigWig conversion, https://deeptools.readthedocs.io/en/latest/content/tools/bamCoverage.html
- `histoneSig` - R package for working with genome files as continuous representations or "signals". https://github.com/semibah/histonesig



## SNP callers

- Liu, Fenglin, Yuanyuan Zhang, Lei Zhang, Ziyi Li, Qiao Fang, Ranran Gao, and Zemin Zhang. "Systematic comparative analysis of single-nucleotide variant detection methods from single-cell RNA sequencing data." Genome Biology 20, no. 1 (2019): 1-15. - https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1863-4 - comparison of seven tools for SNP detection in scRNA-seq data. SAMtools, Strelka2, FreeBayes, and CTAT are best.

- Xu, Chang. “A Review of Somatic Single Nucleotide Variant Calling Algorithms for Next-Generation Sequencing Data.” Computational and Structural Biotechnology Journal 16 (2018): 15–24. https://doi.org/10.1016/j.csbj.2018.01.003. - Overview of 46 somatic Single Nucleotide Variant (SNV) caller tools. Pre-processing, variant evaluation, and post-filtering steps. Four categories of algorithms, description of each, and the corresponding tools: matched tumor-normal (position-, haplotype-, mathine learning-based methods, [Table 1](https://www.sciencedirect.com/science/article/pii/S2001037017300946?via%3Dihub#t0005)), single-sample ([Table 2](https://www.sciencedirect.com/science/article/pii/S2001037017300946?via%3Dihub#t0010), some offer somatic-germline classification), UMI-based (UMI technology, [Figure 1](https://www.sciencedirect.com/science/article/pii/S2001037017300946?via%3Dihub#f0005), [Table 3](https://www.sciencedirect.com/science/article/pii/S2001037017300946?via%3Dihub#t0015)), and RNA-seq (Technology, issues, [Table 4](https://www.sciencedirect.com/science/article/pii/S2001037017300946?via%3Dihub#t0020)) variant calling. Benchmarking using tools for generating synthetic reads, spike-ins, GiAB, melanoma-normal samples, performance evaluation metrics. Issues in representing complex variants and tools for variant normalization. Deep neural network-based algorithms perform best.

- `MutationSeq` - somatic SNV detection from tumor-normal pairs. http://compbio.bccrc.ca/software/mutationSeq/

- `RNA-MuTect` method to detect somatic mutations from tumor DNA - matched normal RNA samples. Applied to TCGA, GTeX data. Most mutated are sun-exposed skin, esophagus mucosa, lung. Number of mutations in RNA is 5-fold larger than in DNA, mutations depend on coverage. Only half of DNA mutations have sufficient coverage in RNA. Filtering using 1) two aligners, 2) removal of errors based on a site-specific error model, 3) removal of RNA editing sites. https://zenodo.org/record/2620062#.XPz9e29KhQI
    - Yizhak, Keren, François Aguet, Jaegil Kim, Julian M Hess, Kirsten Kübler, Jonna Grimsby, Ruslana Frazer, et al. “RNA Sequence Analysis Reveals Macroscopic Somatic Clonal Expansion across Normal Tissues.” HUMAN GENETICS, 2019, 11.



### Deep learning SNP callers

- `DeepVariant` is an analysis pipeline that uses a deep neural network to call genetic variants from next-generation DNA sequencing data. https://github.com/google/deepvariant [Tweet 1](https://twitter.com/acarroll_ATG/status/1194786628759216128?s=20), [Tweet 2](https://twitter.com/brent_p/status/1194792729907093505?s=20)

- `NeuSomatic` - Convolutional neural network (9 layers) for somatic mutations calling. Reads within 7bp window around candidate SNP are extracted, realigned, summarized into matrices for tumor-normal samples,  used for classifying mutation type, length, position. Tested on GiB samples and on DREAM datasets. Comparison with other SNP callers. https://github.com/bioinform/neusomatic
    - Sahraeian, Sayed Mohammad Ebrahim, Ruolin Liu, Bayo Lau, Marghoob Mohiyuddin, and Hugo Y. K. Lam. “Deep Convolutional Neural Networks for Accurate Somatic Mutation Detection,” September 4, 2018. https://doi.org/10.1101/393801.


## SNP annotations

Various genome annotations, [Source: ConsHMM Data availability section](https://www.nature.com/articles/s42003-019-0488-1#data-availability): **25-state chromatin state annotations**: http://compbio.mit.edu/roadmap; **CADD score v1.0**: http://krishna.gs.washington.edu/download/CADD/v1.0/whole_genome_SNVs.tsv.gz; **CADD score v1.4**: http://krishna.gs.washington.edu/download/CADD/v1.4/GRCh37/whole_genome_SNVs.tsv.gz; **CDTS score**: http://www.hli-opendata.com/noncoding/coord_CDTS_percentile_N7794unrelated.txt.gz, http://www.hli-opendata.com/noncoding/SNVusedForCDTScomputation_N7794unrelated_allelicFrequency0.001truncated.txt.gz; **CNEEs** from ref. 22: http://www.stanford.edu/~lowec/data/threePeriods/hg19cnee.bed.gz; **DANN score**: https://cbcl.ics.uci.edu/public_data/DANN/data/; **EIGEN and Eigen-PC score**: https://xioniti01.u.hpc.mssm.edu/v1.1/; **ENCODE DHS**: http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwDnase/; **FATHMM-XF score**: http://fathmm.biocompute.org.uk/fathmm-xf/; **FIRE score**: https://sites.google.com/site/fireregulatoryvariation/; **fitCons score**: http://compgen.cshl.edu/fitCons/0downloads/tracks/i6/scores/; **FunSeq2 score**: http://org.gersteinlab.funseq.s3-website-us-east-1.amazonaws.com/funseq2.1.2/hg19_NCscore_funseq216.tsv.bgz; **GENCODE v19**: https://www.gencodegenes.org/releases/19.html; **GERP++ scores and constrained element calls**: http://mendel.stanford.edu/SidowLab/downloads/gerp/; **GWAS catalog variants**: https://www.ebi.ac.uk/gwas/; **LINSIGHT score**: http://compgen.cshl.edu/~yihuang/tracks/LINSIGHT.bw; **Motif instances and background**: http://compbio.mit.edu/encode-motifs/; REMM score: https://zenodo.org/record/1197579/files/ReMM.v0.3.1.tsv.gz; **Roadmap Epigenomics DHS**: http://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/; **SiPhy-omega and SiPhy-pi constrained element calls (hg19 liftOver)**: https://www.broadinstitute.org/mammals-models/29-mammals-project-supplementary-info

- XGBoost to predict the impact of non-coding variants. Uses functional essentiality features, 3D genome organization, enhancer reporter data, existing deleteriousness metrics (CADD, ncEigen, FATHMM, FunSeq2, LINSIGHT, ORION, ReMM, ncRVIS). Methods how to train/optimize XGBoost. Provides a score, ncER (non-coding essential regulation) nucleotide-resolution and average over 10bp bins.  https://github.com/TelentiLab/ncER_datasets,  https://www.ai-omni.com/
    - Wells, Alex, David Heckerman, Ali Torkamani, Li Yin, Jonathan Sebat, Bing Ren, Amalio Telenti, and Julia di Iulio. “Ranking of Non-Coding Pathogenic Variants and Putative Essential Regions of the Human Genome.” Nature Communications 10, no. 1 (December 2019): 5241. https://doi.org/10.1038/s41467-019-13212-3.

- `ConsHMM` - genome segmentation into 100 conservation states based on a 100 species DNA sequence alignment. Hidden Markov Model, extension of ChromHMM. Compared with 12 other scores (CADD, CDTS, DANN, Eigen, Eigen-PC, FATHMM-XF, FIRE, fitCons, GERP++, PhastCons, PhyloP, REMM, also LINSIGHT, FunSeq2). Clustering, GO enrichment reveals distinct functionality. hg19 single nucleotide scores https://github.com/ernstlab/ConsHMM,  https://figshare.com/articles/ConsHMM_100-state_Segmentation_of_hg19_Human_Genome/8162036/1
    - Arneson, Adriana, and Jason Ernst. “Systematic Discovery of Conservation States for Single-Nucleotide Annotation of the Human Genome.” Communications Biology 2, no. 1 (December 2019): 248. https://doi.org/10.1038/s42003-019-0488-1.

- Predicting pathogenic vs. non-pathogenic SNPs and then regulatory status of each base using 38 functional and structural features. XGBoost model, parameter tuning. All predictions are at https://omni-variants.herokuapp.com/. Existing noncoding deleteriousness metrics: CADD, ncEigen, FATHMM, FunSeq2, LINSIGHT, ORION, ReMM, ncRVIS
    - telenti, amalio, Alexander C Wells, David Heckerman, Ali Torkamani, Bing Ren, and Julia di Iulio. “Identification of Essential Regulatory Elements in the Human Genome.” Preprint. Genomics, October 16, 2018. https://doi.org/10.1101/444562.


- `slivar` - variant expressions, annotation, and filtering, by Brent Pedersen, https://github.com/brentp/slivar

- The Variant Interpretation for Cancer Consortium Meta-Knowledgebase. Aggregate interpretations covering 3,437 unique variants in 415 genes, 357 diseases, and 791 drugs. Validation using GENIE database. https://search.cancervariants.org
    - Wagner, Alex Handler, Brian Walsh, Georgia Mayfield, David Tamborero, Dmitriy Sonkin, Kilannin Krysiak, Jordi Deu Pons, et al. “A Harmonized Meta-Knowledgebase of Clinical Interpretations of Cancer Genomic Variants,” July 11, 2018. https://doi.org/10.1101/366856.

- `vcfanno` - annotate a VCF with other VCFs/BEDs/tabixed files. https://github.com/brentp/vcfanno
    - Pedersen, Brent S., Ryan M. Layer, and Aaron R. Quinlan. “Vcfanno: Fast, Flexible Annotation of Genetic Variants.” Genome Biology 17, no. 1 (December 2016). https://doi.org/10.1186/s13059-016-0973-5.

- `atSNP` - search for effects of SNPs on transcription factor binding. DB of 37 billion variant-motif pairs. Search by SNP IDs, window of SNPs, genomic location, gene, transcription factor. http://atsnp.biostat.wisc.edu/
    - Shin, Sunyoung, Rebecca Hudson, Christopher Harrison, Mark Craven, and Sündüz Keleş. “AtSNP Search: A Web Resource for Statistically Evaluating Influence of Human Genetic Variation on Transcription Factor Binding.” Edited by John Hancock. Bioinformatics, December 8, 2018. https://doi.org/10.1093/bioinformatics/bty1010.


## SNP signatures

- https://github.com/ictic-bioinformatics/CANCERSIGN - CANCERSIGN: a user-friendly and robust tool for identification and classification of mutational signatures and patterns in cancer genomes. Masroor Bayati, Hamid Reza Rabiee, Mehrdad Mehrbod, Fatemeh Vafaee, Diako Ebrahimi, Alistair Forrest, Hamid Alinejad-Rokny. bioRxiv 424960; doi: https://doi.org/10.1101/424960

- https://github.com/danro9685/SparseSignatures, https://bioconductor.org/packages/release/bioc/html/SparseSignatures.html - De Novo Mutational Signature Discovery in Tumor Genomes using SparseSignatures. Daniele Ramazzotti, Avantika Lal, Keli Liu, Robert Tibshirani, Arend Sidow. bioRxiv 384834; doi: https://doi.org/10.1101/384834

- `SignatureAnalyzer` - finding mutation patterns in multiple samples, NMF. https://software.broadinstitute.org/cancer/cga/msp

- https://bioconductor.org/packages/release/bioc/html/YAPSA.html - Yet Another Package for Signature Analysis, functionality used in L. Alexandrov et al., Nature 2013


## SNP pathogenicity scores

- `regBase` - Prediction of regulatory impact of variants outside of protein-coding regions, human. Trained on prediction scores from 23 tools, Gradient Tree Boosting, thorough training and evaluation. hg19 predictions are available for download. Python implementation https://github.com/mulinlab/regBase
    - Zhang, Shijie, Yukun He, Huanhuan Liu, Haoyu Zhai, Dandan Huang, Xianfu Yi, Xiaobao Dong, et al. “RegBase: Whole Genome Base-Wise Aggregation and Functional Prediction for Human Non-Coding Regulatory Variants.” Nucleic Acids Research, September 12, 2019, gkz774. https://doi.org/10.1093/nar/gkz774.

- `ClinPred` - pathogenicity prediction for all nonsynonymous SNPs. Trained on ClinVar, validated on nine other databases. Random forest and gradient boosted decision tree, comparison with other machine learning algorithms. Downloadable scores for all nonsynonymous SNPs, https://sites.google.com/site/clinpred/home
    - Alirezaie, Najmeh, Kristin D. Kernohan, Taila Hartley, Jacek Majewski, and Toby Dylan Hocking. “ClinPred: Prediction Tool to Identify Disease-Relevant Nonsynonymous Single-Nucleotide Variants.” The American Journal of Human Genetics 103, no. 4 (October 2018): 474–83. https://doi.org/10.1016/j.ajhg.2018.08.005.

## SNP visualization

- `CoMutPlotter` - plotting cancer mutational profiles. Supports VCF, MAF, TSV. http://tardis.cgu.edu.tw/comutplotter/
    - Huang, Po-Jung, Hou-Hsien Lin, Chi-Ching Lee, Ling-Ya Chiu, Shao-Min Wu, Yuan-Ming Yeh, Petrus Tang, Cheng-Hsun Chiu, Ping-Chiang Lyu, and Pei-Chien Tsai. “CoMutPlotter: A Web Tool for Visual Summary of Mutations in Cancer Cohorts.” BMC Medical Genomics 12, no. S5 (July 2019): 99. https://doi.org/10.1186/s12920-019-0510-y.

- `gwasTools` - A collection of R scripts that might be useful for exploring and plotting GWAS results. https://github.com/bnwolford/gwasTools

- `gpart` - R package for defining LD blocks (Big-LD algorithm), and visualizing them. https://bioconductor.org/packages/release/bioc/html/gpart.html
    - Ah Kim, Sun, Myriam Brossard, Delnaz Roshandel, Andrew D Paterson, Shelley B Bull, and Yun Joo Yoo. “Gpart: Human Genome Partitioning and Visualization of High-Density SNP Data by Identifying Haplotype Blocks.” Edited by Alfonso Valencia. Bioinformatics, May 9, 2019, btz308. https://doi.org/10.1093/bioinformatics/btz308.

- `maftools` - Summarize, Analyze and Visualize MAF files from TCGA or in house studies. Bioconductor, https://bioconductor.org/packages/release/bioc/vignettes/maftools/inst/doc/maftools.html, and GitHub, https://github.com/PoisonAlien/maftools

- `manhattanly` - Interactive Manhattan plots, https://cran.r-project.org/web/packages/manhattanly/

- `mutcraft` - R tools to mine & craft somatic mutations from cancer genomes, https://github.com/EmilieT/mutcraft

- `MutScan` - Detect and visualize target mutations by scanning FastQ files directly. https://github.com/OpenGene/MutScan

- `samplot` - Plot structural variant signals from many BAMs and CRAMs. https://github.com/ryanlayer/samplot

- `ttplot` - Tao Yan's Plot Toolkit, plots LD Heatmap, Manhattan plot. https://github.com/YTLogos/ttplot

- `VIVA` - VCF visualization tool, written in Julia. Competing tools - vcfR, IGVZ, Genome Browser, Genome Savant, svviz, jvarkit - JfxNgs. Input - VCF file and, optionally, variant list, sample list, sample metadata. Filtering. Heatmap visualization. https://github.com/compbiocore/VariantVisualization.jl
    - Tollefson, George A, Jessica Schuster, Fernando Gelin, Ashok Ragavendran, Isabel Restrepo, Paul Stey, James Padbury, and Alper Uzun. “VIVA (VIsualization of VAriants): A VCF File Visualization Tool.” BioRxiv, March 28, 2019. https://doi.org/10.1101/589879.


## SNP databases

- `GWASatlas` resource, analysis of pleiotropy, genetic architecture of complex traits.  https://atlas.ctglab.nl/
    - Watanabe, Kyoko, Sven Stringer, Oleksandr Frei, Masa Umićević Mirkov, Tinca J.C. Polderman, Sophie van der Sluis, Ole A. Andreassen, Benjamin M. Neale, and Danielle Posthuma. “A Global View of Pleiotropy and Genetic Architecture in Complex Traits.” BioRxiv, January 1, 2018, 500090. https://doi.org/10.1101/500090.

- `GWAScentral` - central GWAS repository. Browser, download. https://www.gwascentral.org
    - Beck, Tim, Tom Shorter, and Anthony J Brookes. “GWAS Central: A Comprehensive Resource for the Discovery and Comparison of Genotype and Phenotype Data from Genome-Wide Association Studies.” Nucleic Acids Research, October 15, 2019, gkz895. https://doi.org/10.1093/nar/gkz895.

- `clinvar` - This repo provides tools to convert ClinVar data into a tab-delimited flat file, and also provides that resulting tab-delimited flat file. https://github.com/macarthur-lab/clinvar

- Clinical Interpretation of Variants in Cancer database, http://www.civicdb.org/. CIViC interface public API, http://griffithlab.org/civic-api-docs/

## InDels

- Kosugi, Shunichi, Yukihide Momozawa, Xiaoxi Liu, Chikashi Terao, Michiaki Kubo, and Yoichiro Kamatani. “Comprehensive Evaluation of Structural Variation Detection Algorithms for Whole Genome Sequencing.” Genome Biology 20, no. 1 (December 2019): 117. https://doi.org/10.1186/s13059-019-1720-5. - Benchmarking of structural variant detection tools. Introduction to types of structural variants. No tool detects all. Table 1 prioritizes best tools for deletion, duplication, insertion, invertion detection.

- `Pindel` - breakpoints of large deletions, medium sized insertions, inversions, tandem duplications and other structural variants. http://gmt.genome.wustl.edu/packages/pindel/

- `Dindel` - Accurate indel calls from short-read data. http://www.sanger.ac.uk/science/tools/dindel

- `Destruct` - joint prediction of rearrangement breakpoints from single or multiple tumour samples. https://bitbucket.org/dranew/destruct.git

- `MindTheGap` - detection and assembly of DNA insertion variants, https://gatb.inria.fr/software/mind-the-gap/, https://github.com/GATB/MindTheGap

- `Breakfast` - a software for detecting genomic structural variants from DNA sequencing data, https://github.com/annalam/breakfast

- `SVAFotate` - Annotate a (lumpy) structual variant (SV) VCF with allele frequencies (AFs) from large population SV cohorts. https://github.com/fakedrtom/SVAFotate


## CNV, SV

- `CNVnator` - a tool for CNV discovery and genotyping from depth-of-coverage by mapped reads. https://github.com/abyzovlab/CNVnator

- `TITAN` - a tool for predicting subclonal copy number alterations (CNA) and loss of heterozygosity (LOH) from tumour whole genome sequencing data. http://compbio.bccrc.ca/software/titan/

- `QDNASeq` - Quantitative DNA sequencing for chromosomal aberrations. https://bioconductor.org/packages/release/bioc/html/QDNAseq.html

- `samplot` - Plot structural variant signals from many BAMs and CRAMs. https://github.com/ryanlayer/samplot

- `smoove` - structural variant calling and genotyping with existing tools, but, smoothly. https://github.com/brentp/smoove

- `SynthEx` - CNV detection from exome and whole genome sequencing.

- `HATCHet` (Holistic Allele-specific Tumor Copy-number Heterogeneity) is an algorithm that infers allele and clone-specific CNAs and WGDs jointly across multiple tumor samples from the same patient, and that leverages the relationships between clones in these samples. https://github.com/raphael-group/hatchet
    - Zaccaria, Simone, and Benjamin J. Raphael. “Accurate Quantification of Copy-Number Aberrations and Whole-Genome Duplications in Multi-Sample Tumor Sequencing Data.” BioRxiv, January 1, 2018, 496174. https://doi.org/10.1101/496174.

- `CNVkit` - capturing CNVs in on-target and off-target genomic regions. Existing tools (CNVer, ExomeCNV, exomeCopy, CONTRA, CoNIFER, ExomeDepth, VarScan2, XHMM, ngCGH, EXCAVATOR, CANOES, PatternCNV, CODEX, Control-FREEC, cn.MOPS, cnvOffSeq, CopyWriteR). Account for GC content, mappability. Python 2.7 implementation. https://github.com/etal/cnvkit, https://github.com/etal/cnvkit-examples
    - Talevich, Eric, A. Hunter Shain, Thomas Botton, and Boris C. Bastian. “CNVkit: Genome-Wide Copy Number Detection and Visualization from Targeted DNA Sequencing.” PLOS Computational Biology 12, no. 4 (April 21, 2016): e1004873. https://doi.org/10.1371/journal.pcbi.1004873.

- `Manta` - SV detection in single- and tumor-normal samples. parallelized for within-sample performance. Fast, detects more variants of different types. https://github.com/Illumina/manta
    - Chen, Xiaoyu, Ole Schulz-Trieglaff, Richard Shaw, Bret Barnes, Felix Schlesinger, Morten Källberg, Anthony J. Cox, Semyon Kruglyak, and Christopher T. Saunders. “Manta: Rapid Detection of Structural Variants and Indels for Germline and Cancer Sequencing Applications.” Bioinformatics 32, no. 8 (April 15, 2016): 1220–22. https://doi.org/10.1093/bioinformatics/btv710.

- `Control-FREEC` - assess copy number and genotype information in whole genome and exome sequencing data. Corrects for contamination by normal cells and variable sample ploidy. With a matched normal sample, distinguishes somatic from germline events. http://boevalab.com/tools.html
    - Boeva, Valentina, Tatiana Popova, Kevin Bleakley, Pierre Chiche, Julie Cappo, Gudrun Schleiermacher, Isabelle Janoueix-Lerosey, Olivier Delattre, and Emmanuel Barillot. “Control-FREEC: A Tool for Assessing Copy Number and Allelic Content Using next-Generation Sequencing Data.” Bioinformatics (Oxford, England) 28, no. 3 (February 1, 2012): 423–25. https://doi.org/10.1093/bioinformatics/btr670.



## Miscellaneous

 Awesome papers and projects about CNV and SV using NGS data. https://github.com/geocarvalho/sv-cnv-studies

- DNA sequencing analysis notes from Ming Tang. https://github.com/crazyhottommy/DNA-seq-analysis

- `SNPhylo` - A pipeline to generate a phylogenetic tree from huge SNP data, https://github.com/thlee/SNPhylo, http://chibba.pgml.uga.edu/snphylo/

- Application for making ENCODE Blacklists, and links to canonical blacklists, https://github.com/Boyle-Lab/Blacklist
    - Blacklist citation: Amemiya, Haley M., Anshul Kundaje, and Alan P. Boyle. “The ENCODE Blacklist: Identification of Problematic Regions of the Genome.” Scientific Reports 9, no. 1 (December 2019): 9354. https://doi.org/10.1038/s41598-019-45839-z.


- HOT/XOT regions. The high occupancy target (HOT) and extreme occupancy target (XOT) regions in all contexts were downloaded through the ENCODE data portal at http://encode-ftp.s3.amazonaws.com/modENCODE_VS_ENCODE/Regulation/Human/hotRegions/maphot_hs_selection_reg_cx_simP05_all.bed and http://encode-ftp.s3.amazonaws.com/modENCODE_VS_ENCODE/Regulation/Human/hotRegions/maphot_hs_selection_reg_cx_simP01_all.bed (hg38 ?). Potential source

- `GEM` - mappability calculations for each genomic region, accounting for mismatches. Pre-calculated UCSC genome browser tracks for human and mouse. Mappability of genes, both protein-coding and non-protein coding. RPKUM - unique exons for quantifying gene expression. https://sourceforge.net/projects/gemlibrary/files/gem-library/
    - Derrien, Thomas, Jordi Estellé, Santiago Marco Sola, David G. Knowles, Emanuele Raineri, Roderic Guigó, and Paolo Ribeca. “Fast Computation and Applications of Genome Mappability.” PloS One 7, no. 1 (2012): e30377. https://doi.org/10.1371/journal.pone.0030377.

- `refGenie` - reference genome manager. http://refgenie.databio.org/en/latest/

- `genomepy` - Download genomes the easy way. https://github.com/simonvh/genomepy

- `GeneticsDesign` - GWAS power analysis, functions for designing genetics studies, https://www.bioconductor.org/packages/release/bioc/html/GeneticsDesign.html

- Sample swap check. https://github.com/parklab/NGSCheckMate, https://github.com/brentp/somalier
