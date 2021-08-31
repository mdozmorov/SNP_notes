# SNP-related notes

[![MIT License](https://img.shields.io/apm/l/atomic-design-ui.svg?)](https://github.com/tterb/atomic-design-ui/blob/master/LICENSEs) [![PR's Welcome](https://img.shields.io/badge/PRs-welcome-brightgreen.svg?style=flat)](http://makeapullrequest.com) 

SNP- and genome variation-related tools and genomics data analysis resources. Please, [contribute and get in touch](CONTRIBUTING.md)! See [MDmisc notes](https://github.com/mdozmorov/MDmisc_notes) for other programming and genomics-related notes.

# Table of content

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->


- [Variant calling pipelines](#variant-calling-pipelines)
  - [Preprocessing tools](#preprocessing-tools)
  - [Depth](#depth)
- [SNP callers](#snp-callers)
  - [Deep learning SNP callers](#deep-learning-snp-callers)
- [SNP annotations](#snp-annotations)
- [SNP signatures](#snp-signatures)
- [SNP pathogenicity scores](#snp-pathogenicity-scores)
- [SNP visualization, clustering](#snp-visualization-clustering)
- [SNP, GWAS databases](#snp--gwas-databases)
- [InDels](#indels)
- [CNV, SV](#cnv-sv)
- [Power](#power)
- [Miscellaneous](#miscellaneous)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## Variant calling pipelines

- Koboldt, Daniel C. “[Best Practices for Variant Calling in Clinical Sequencing](https://doi.org/10.1186/s13073-020-00791-w).” Genome Medicine, (December 2020) - Introduction in genomic variant calling, panel/exome/whole genome sequencing technologies (Table 1), preprocessing, analysis (SNVs/indels, mutations, CNVs, SVs, gene fusions, Table 2), gold standard datasets (GIAB), best practices, filtering for each type of genomic variant.

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

- [Rust-bio-tools](https://github.com/rust-bio/rust-bio-tools) - VCF matching, conversion to text, report, FASTQ split/filter, BAM depth, merging.

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

- [Varlociraptor](https://varlociraptor.github.io) - a unifying statistical model allowing for the detection of SNVs, MNVs, InDels, Inversions, Duplications, Breakends. Explicitly controls FDR. Defines variant calling grammar, allowing for defining all types of variants. Tumor-normal comparison allows for classifying variants into germline, somatic variants. Evaluated on simulated (artificial clones) and real (Venter's genome) data. Outperforms six other variant callers. Can call variants in single samples, in RNA-seq data, FFPE attifact detection. [Tweet](https://twitter.com/johanneskoester/status/1331218560367128577?s=20)
    - Köster, Johannes, Louis J. Dijkstra, Tobias Marschall, and Alexander Schönhuth. “[Varlociraptor: Enhancing Sensitivity and Controlling False Discovery Rate in Somatic Indel Discovery](https://doi.org/10.1186/s13059-020-01993-6).” Genome Biology 21, no. 1 (December 2020)

- Liu, Fenglin, Yuanyuan Zhang, Lei Zhang, Ziyi Li, Qiao Fang, Ranran Gao, and Zemin Zhang. "Systematic comparative analysis of single-nucleotide variant detection methods from single-cell RNA sequencing data." Genome Biology 20, no. 1 (2019): 1-15. - https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1863-4 - comparison of seven tools for SNP detection in scRNA-seq data. SAMtools, Strelka2, FreeBayes, and CTAT are best.

- Xu, Chang. “A Review of Somatic Single Nucleotide Variant Calling Algorithms for Next-Generation Sequencing Data.” Computational and Structural Biotechnology Journal 16 (2018): 15–24. https://doi.org/10.1016/j.csbj.2018.01.003. - Overview of 46 somatic Single Nucleotide Variant (SNV) caller tools. Pre-processing, variant evaluation, and post-filtering steps. Four categories of algorithms, description of each, and the corresponding tools: matched tumor-normal (position-, haplotype-, mathine learning-based methods, [Table 1](https://www.sciencedirect.com/science/article/pii/S2001037017300946?via%3Dihub#t0005)), single-sample ([Table 2](https://www.sciencedirect.com/science/article/pii/S2001037017300946?via%3Dihub#t0010), some offer somatic-germline classification), UMI-based (UMI technology, [Figure 1](https://www.sciencedirect.com/science/article/pii/S2001037017300946?via%3Dihub#f0005), [Table 3](https://www.sciencedirect.com/science/article/pii/S2001037017300946?via%3Dihub#t0015)), and RNA-seq (Technology, issues, [Table 4](https://www.sciencedirect.com/science/article/pii/S2001037017300946?via%3Dihub#t0020)) variant calling. Benchmarking using tools for generating synthetic reads, spike-ins, GiAB, melanoma-normal samples, performance evaluation metrics. Issues in representing complex variants and tools for variant normalization. Deep neural network-based algorithms perform best.

- `MutationSeq` - somatic SNV detection from tumor-normal pairs. http://compbio.bccrc.ca/software/mutationSeq/

- `Samovar` - mosaic SNV calling from single samples sequenced using WGS 10X linked reads technology. Intro into mosaic SNVs, linked reads. Outperforms MuTect2, MosaicHunter. https://github.com/cdarby/samovar
    - Darby, Charlotte A, James R Fitch, Patrick J Brennan, Benjamin J Kelly, Natalie Bir, Vincent Magrini, Jeffrey Leonard, et al. “Samovar: Single-Sample Mosaic SNV Calling with Linked Reads.” Preprint. Genomics, February 25, 2019. https://doi.org/10.1101/560532.

- `RNA-MuTect` method to detect somatic mutations from tumor DNA - matched normal RNA samples. Applied to TCGA, GTeX data. Most mutated are sun-exposed skin, esophagus mucosa, lung. Number of mutations in RNA is 5-fold larger than in DNA, mutations depend on coverage. Only half of DNA mutations have sufficient coverage in RNA. Filtering using 1) two aligners, 2) removal of errors based on a site-specific error model, 3) removal of RNA editing sites. https://zenodo.org/record/2620062#.XPz9e29KhQI
    - Yizhak, Keren, François Aguet, Jaegil Kim, Julian M Hess, Kirsten Kübler, Jonna Grimsby, Ruslana Frazer, et al. “RNA Sequence Analysis Reveals Macroscopic Somatic Clonal Expansion across Normal Tissues.” HUMAN GENETICS, 2019, 11.



### Deep learning SNP callers

- [DeepVariant](https://github.com/google/deepvariant) is an analysis pipeline that uses a deep neural network to call genetic variants from next-generation DNA sequencing data. CNN trained to detect SNPs. Inception v2 architecture. Outperforms GATK and others in FDA challenge, GIB data. Docker/Singularity CPU and GPU versions, local and cloud compatible. [Tweet 1](https://twitter.com/acarroll_ATG/status/1194786628759216128?s=20), [Tweet 2](https://twitter.com/brent_p/status/1194792729907093505?s=20), [GitHub repo](https://github.com/google/deepvariant/)
    - Poplin, Ryan, Pi-Chuan Chang, David Alexander, Scott Schwartz, Thomas Colthurst, Alexander Ku, Dan Newburger, et al. “[Creating a Universal SNP and Small Indel Variant Caller with Deep Neural Networks,](https://doi.org/10.1101/092890)” March 20, 2018

- `NeuSomatic` - Convolutional neural network (9 layers) for somatic mutations calling. Reads within 7bp window around candidate SNP are extracted, realigned, summarized into matrices for tumor-normal samples,  used for classifying mutation type, length, position. Tested on GiB samples and on DREAM datasets. Comparison with other SNP callers. https://github.com/bioinform/neusomatic
    - Sahraeian, Sayed Mohammad Ebrahim, Ruolin Liu, Bayo Lau, Marghoob Mohiyuddin, and Hugo Y. K. Lam. “Deep Convolutional Neural Networks for Accurate Somatic Mutation Detection,” September 4, 2018. https://doi.org/10.1101/393801.

- [DeepSNV](https://bioconductor.org/packages/deepSNV/) - an R package for detecting clonal and subclonal variants from deeply sequenced clonally heterogeneous cancer data. The deepSNV algorithm is used for a comparative setup with a control experiment of the same loci and uses a beta-binomial model and a likelihood ratio test to discriminate sequencing errors and subclonal SNVs. The shearwater algorithm computes a Bayes classifier based on a beta-binomial model for variant calling with multiple samples for precisely estimating model parameters - such as local error rates and dispersion - and prior knowledge, e.g. from variation data bases such as COSMIC. 
    - Gerstung, M., E. Papaemmanuil, and P. J. Campbell. “[Subclonal Variant Calling with Multiple Samples and Prior Knowledge](https://doi.org/10.1093/bioinformatics/btt750).” Bioinformatics, (May 1, 2014)

## SNP annotations

Various genome annotations, [Source: ConsHMM Data availability section](https://www.nature.com/articles/s42003-019-0488-1#data-availability): [25-state chromatin state annotations](http://compbio.mit.edu/roadmap); [CADD score v1.0](http://krishna.gs.washington.edu/download/CADD/v1.0/whole_genome_SNVs.tsv.gz); [CADD score v1.4](http://krishna.gs.washington.edu/download/CADD/v1.4/GRCh37/whole_genome_SNVs.tsv.gz); [CDTS score](http://www.hli-opendata.com/noncoding/coord_CDTS_percentile_N7794unrelated.txt.gz), [another link](http://www.hli-opendata.com/noncoding/SNVusedForCDTScomputation_N7794unrelated_allelicFrequency0.001truncated.txt.gz); [CNEEs](http://www.stanford.edu/~lowec/data/threePeriods/hg19cnee.bed.gz); [DANN score](https://cbcl.ics.uci.edu/public_data/DANN/data/); [EIGEN and Eigen-PC score](https://xioniti01.u.hpc.mssm.edu/v1.1/); [ENCODE DHS](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwDnase/); [FATHMM-XF score](http://fathmm.biocompute.org.uk/fathmm-xf/); [FIRE score](https://sites.google.com/site/fireregulatoryvariation/); [fitCons score](http://compgen.cshl.edu/fitCons/0downloads/tracks/i6/scores/); [FunSeq2 score](http://org.gersteinlab.funseq.s3-website-us-east-1.amazonaws.com/funseq2.1.2/hg19_NCscore_funseq216.tsv.bgz); [GENCODE v19](https://www.gencodegenes.org/releases/19.html); [GERP++ scores and constrained element calls](http://mendel.stanford.edu/SidowLab/downloads/gerp/); [GWAS catalog variants](https://www.ebi.ac.uk/gwas/); [LINSIGHT score](http://compgen.cshl.edu/~yihuang/tracks/LINSIGHT.bw); [Motif instances and background](http://compbio.mit.edu/encode-motifs/); [REMM score](https://zenodo.org/record/1197579/files/ReMM.v0.3.1.tsv.gz); [Roadmap Epigenomics DHS](http://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/); [SiPhy-omega and SiPhy-pi constrained element calls (hg19 liftOver)](https://www.broadinstitute.org/mammals-models/29-mammals-project-supplementary-info)

- [VCF-plotein](https://vcfplotein.liigh.unam.mx/#/) - graphical, interactive interptetation of exome sequencing data in VCF format. Includes major databases, from gnomAD, COSMIC to Human Phenotype Ontology and GO terms. [Shiny app](https://vcfplotein.liigh.unam.mx/#/), [GitHub](https://github.com/raulossio/VCF-plotein)
    - Ossio, Raul, O Isaac Garcia-Salinas, Diego Said Anaya-Mancilla, Jair S Garcia-Sotelo, Luis A Aguilar, David J Adams, and Carla Daniela Robles-Espinoza. “[VCF/Plotein: Visualisation and Prioritisation of Genomic Variants from Human Exome Sequencing Projects](https://doi.org/10.1093/bioinformatics/btz458).” Bioinformatics, June 4, 2019

- [The Personal Cancer Genome Reporter (PCGR)](https://github.com/sigven/pcgr) is a stand-alone software package for functional annotation and translation of individual cancer genomes for precision cancer medicine. Currently, it interprets both somatic SNVs/InDels and copy number aberrations. The software extends basic gene and variant annotations from the [Ensembl’s Variant Effect Predictor (VEP)](http://www.ensembl.org/info/docs/tools/vep/index.html) with oncology-relevant, up-to-date annotations retrieved flexibly through vcfanno, and produces interactive HTML reports intended for clinical interpretation. NOTE: If you also want to interrogate the clinical impact of germline variants in the same individual, try the accompanying tool [Cancer Predisposition Sequencing Reporter (CPSR)](https://github.com/sigven/cpsr).
    - Sigve Nakken, Ghislain Fournous, Daniel Vodák, Lars Birger Aaasheim, Ola Myklebost, and Eivind Hovig. [Personal Cancer Genome Reporter: variant interpretation report for precision oncology](https://doi.org/10.1093/bioinformatics/btx817) (2017). Bioinformatics

- `FUMA` - Functional Mapping and Annotation of GWAS using 43 scRNA-seq datasets (human and mouse). MAGMA Cell type-specific enrichment. Applied to 26 GWAS disorders,  https://fuma.ctglab.nl/. Processed data and instructions for self-download: https://github.com/Kyoko-wtnb/FUMA_scRNA_data
    - Watanabe, Kyoko, Maša Umićević Mirkov, Christiaan A. de Leeuw, Martijn P. van den Heuvel, and Danielle Posthuma. “Genetic Mapping of Cell Type Specificity for Complex Traits.” Nature Communications 10, no. 1 (December 2019): 3222. https://doi.org/10.1038/s41467-019-11181-1.

- XGBoost to predict the impact of non-coding variants. Uses functional essentiality features, 3D genome organization, enhancer reporter data, existing deleteriousness metrics (CADD, ncEigen, FATHMM, FunSeq2, LINSIGHT, ORION, ReMM, ncRVIS). Methods how to train/optimize XGBoost. Provides a score, ncER (non-coding essential regulation) nucleotide-resolution and average over 10bp bins.  https://github.com/TelentiLab/ncER_datasets,  https://www.ai-omni.com/
    - Wells, Alex, David Heckerman, Ali Torkamani, Li Yin, Jonathan Sebat, Bing Ren, Amalio Telenti, and Julia di Iulio. “Ranking of Non-Coding Pathogenic Variants and Putative Essential Regions of the Human Genome.” Nature Communications 10, no. 1 (December 2019): 5241. https://doi.org/10.1038/s41467-019-13212-3.

- `ConsHMM` - genome segmentation into 100 conservation states based on a 100 species DNA sequence alignment. Hidden Markov Model, extension of ChromHMM. Compared with 12 other scores (CADD, CDTS, DANN, Eigen, Eigen-PC, FATHMM-XF, FIRE, fitCons, GERP++, PhastCons, PhyloP, REMM, also LINSIGHT, FunSeq2). Clustering, GO enrichment reveals distinct functionality. hg19 single nucleotide scores https://github.com/ernstlab/ConsHMM,  https://figshare.com/articles/ConsHMM_100-state_Segmentation_of_hg19_Human_Genome/8162036/1
    - Arneson, Adriana, and Jason Ernst. “Systematic Discovery of Conservation States for Single-Nucleotide Annotation of the Human Genome.” Communications Biology 2, no. 1 (December 2019): 248. https://doi.org/10.1038/s42003-019-0488-1.

- Predicting pathogenic vs. non-pathogenic SNPs and then regulatory status of each base using 38 functional and structural features. XGBoost model, parameter tuning. All predictions are at https://omni-variants.herokuapp.com/. Existing noncoding deleteriousness metrics: CADD, ncEigen, FATHMM, FunSeq2, LINSIGHT, ORION, ReMM, ncRVIS
    - telenti, amalio, Alexander C Wells, David Heckerman, Ali Torkamani, Bing Ren, and Julia di Iulio. “Identification of Essential Regulatory Elements in the Human Genome.” Preprint. Genomics, October 16, 2018. https://doi.org/10.1101/444562.

- [Seattle](https://snp.gs.washington.edu/SeattleSeqAnnotation154/) - The SeattleSeq Annotation server provides annotation of SNVs (single-nucleotide variations) and small indels, both known and novel. This annotation includes dbSNP rs IDs, gene names and accession numbers, variation functions (e.g. missense), protein positions and amino-acid changes, conservation scores, HapMap frequencies, PolyPhen predictions, and clinical association

- [snpEff](https://pcingola.github.io/SnpEff/se_introduction/) - a variant annotation and effect prediction tool

- `slivar` - variant expressions, annotation, and filtering, by Brent Pedersen, https://github.com/brentp/slivar

- The Variant Interpretation for Cancer Consortium Meta-Knowledgebase. Aggregate interpretations covering 3,437 unique variants in 415 genes, 357 diseases, and 791 drugs. Validation using GENIE database. https://search.cancervariants.org
    - Wagner, Alex Handler, Brian Walsh, Georgia Mayfield, David Tamborero, Dmitriy Sonkin, Kilannin Krysiak, Jordi Deu Pons, et al. “A Harmonized Meta-Knowledgebase of Clinical Interpretations of Cancer Genomic Variants,” July 11, 2018. https://doi.org/10.1101/366856.

- `vcfanno` - annotate a VCF with other VCFs/BEDs/tabixed files. https://github.com/brentp/vcfanno
    - Pedersen, Brent S., Ryan M. Layer, and Aaron R. Quinlan. “Vcfanno: Fast, Flexible Annotation of Genetic Variants.” Genome Biology 17, no. 1 (December 2016). https://doi.org/10.1186/s13059-016-0973-5.

- `atSNP` - search for effects of SNPs on transcription factor binding. DB of 37 billion variant-motif pairs. Search by SNP IDs, window of SNPs, genomic location, gene, transcription factor. http://atsnp.biostat.wisc.edu/
    - Shin, Sunyoung, Rebecca Hudson, Christopher Harrison, Mark Craven, and Sündüz Keleş. “AtSNP Search: A Web Resource for Statistically Evaluating Influence of Human Genetic Variation on Transcription Factor Binding.” Edited by John Hancock. Bioinformatics, December 8, 2018. https://doi.org/10.1093/bioinformatics/bty1010.


## SNP signatures

- `sigminer` - Genomic Alteration Signature Analysis in R, developed by Shixiang Wang, novel and known signature extraction, visualization. Code: https://github.com/ShixiangWang/sigminer, Documentation: https://shixiangwang.github.io/sigminer-doc/
    - `sigminer.prediction` - Train and Predict Cancer Subtype with Keras Model based on Mutational Signatures, https://github.com/ShixiangWang/sigminer.prediction

- https://github.com/ictic-bioinformatics/CANCERSIGN - CANCERSIGN: a user-friendly and robust tool for identification and classification of mutational signatures and patterns in cancer genomes. Masroor Bayati, Hamid Reza Rabiee, Mehrdad Mehrbod, Fatemeh Vafaee, Diako Ebrahimi, Alistair Forrest, Hamid Alinejad-Rokny. bioRxiv 424960; doi: https://doi.org/10.1101/424960

- https://github.com/danro9685/SparseSignatures, https://bioconductor.org/packages/release/bioc/html/SparseSignatures.html - De Novo Mutational Signature Discovery in Tumor Genomes using SparseSignatures. Daniele Ramazzotti, Avantika Lal, Keli Liu, Robert Tibshirani, Arend Sidow. bioRxiv 384834; doi: https://doi.org/10.1101/384834

- `SignatureAnalyzer` - finding mutation patterns in multiple samples, NMF. https://software.broadinstitute.org/cancer/cga/msp

- https://bioconductor.org/packages/release/bioc/html/YAPSA.html - Yet Another Package for Signature Analysis, functionality used in L. Alexandrov et al., Nature 2013


## SNP pathogenicity scores

- `SEMpl` - predict the impact of SNPs on TF binding. Uses ChIP-seq data, DNAse-seq, and PWMs. Simulate all possible SNPs in PWMs, estimate effect on ChIP signal. https://github.com/Boyle-Lab/SEM_CPP
    - Nishizaki, Sierra S, Natalie Ng, Shengcheng Dong, Robert S Porter, Cody Morterud, Colten Williams, Courtney Asman, Jessica A Switzenberg, and Alan P Boyle. “Predicting the Effects of SNPs on Transcription Factor Binding Affinity.” Edited by John Hancock. Bioinformatics, August 2, 2019, btz612. https://doi.org/10.1093/bioinformatics/btz612.

- `regBase` - Prediction of regulatory impact of variants outside of protein-coding regions, human. Trained on prediction scores from 23 tools, Gradient Tree Boosting, thorough training and evaluation. hg19 predictions are available for download. Python implementation https://github.com/mulinlab/regBase
    - Zhang, Shijie, Yukun He, Huanhuan Liu, Haoyu Zhai, Dandan Huang, Xianfu Yi, Xiaobao Dong, et al. “RegBase: Whole Genome Base-Wise Aggregation and Functional Prediction for Human Non-Coding Regulatory Variants.” Nucleic Acids Research, September 12, 2019, gkz774. https://doi.org/10.1093/nar/gkz774.

- `ClinPred` - pathogenicity prediction for all nonsynonymous SNPs. Trained on ClinVar, validated on nine other databases. Random forest and gradient boosted decision tree, comparison with other machine learning algorithms. Downloadable scores for all nonsynonymous SNPs, https://sites.google.com/site/clinpred/home
    - Alirezaie, Najmeh, Kristin D. Kernohan, Taila Hartley, Jacek Majewski, and Toby Dylan Hocking. “ClinPred: Prediction Tool to Identify Disease-Relevant Nonsynonymous Single-Nucleotide Variants.” The American Journal of Human Genetics 103, no. 4 (October 2018): 474–83. https://doi.org/10.1016/j.ajhg.2018.08.005.

## SNP visualization, clustering

- [CONQUER](https://github.com/roderickslieker/CONQUER) - an R package for visualizing individual and multiple SNPs in epigenomic context, omics expression values, analysis of QTL genes for pathway enrichment. Interactive D3 graphics, circos plots. Similar functionality - [DEPICT](https://data.broadinstitute.org/mpg/depict/index.html) Input - rsIDs. 
    - Bouland, Gerard A, Joline W J Beulens, Joey Nap, and Arnaud Zaldumbide. “[CONQUER: An Interactive Toolbox to Understand Functional Consequences of GWAS Hits](https://doi.org/10.1093/nargab/lqaa085),” NAR Genomics and Bioinformatics, October 27, 2021, 7.

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

- [VarClust](https://github.com/fasterius/VarClust) - A Python package for clustering of single nucleotide variants from high-through seqencing data. Works on single-sample VCF files.

- `VIVA` - VCF visualization tool, written in Julia. Competing tools - vcfR, IGVZ, Genome Browser, Genome Savant, svviz, jvarkit - JfxNgs. Input - VCF file and, optionally, variant list, sample list, sample metadata. Filtering. Heatmap visualization. https://github.com/compbiocore/VariantVisualization.jl
    - Tollefson, George A, Jessica Schuster, Fernando Gelin, Ashok Ragavendran, Isabel Restrepo, Paul Stey, James Padbury, and Alper Uzun. “VIVA (VIsualization of VAriants): A VCF File Visualization Tool.” BioRxiv, March 28, 2019. https://doi.org/10.1101/589879.


## SNP, GWAS databases

- Havrilla, James M., Brent S. Pedersen, Ryan M. Layer, and Aaron R. Quinlan. “[A Map of Constrained Coding Regions in the Human Genome](https://doi.org/10.1038/s41588-018-0294-6).” Nature Genetics, (2019) - Constrained coding regions (CCRs), analysis of gnomAD. [GitHub](https://github.com/quinlan-lab/ccr). [Data availability](https://www.nature.com/articles/s41588-018-0294-6#data-availability) section contain links to many genomics datasets. [Supplementary data](https://www.nature.com/articles/s41588-018-0294-6#Sec28):
    - [Supplementary Table 1](https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-018-0294-6/MediaObjects/41588_2018_294_MOESM3_ESM.txt) - Genes with CCRs in the 99th percentile or higher
    - [Supplementary Table 2](https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-018-0294-6/MediaObjects/41588_2018_294_MOESM4_ESM.txt) - CCRs under purifying selection specifically in humans
    - [Supplementary Table 3](https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-018-0294-6/MediaObjects/41588_2018_294_MOESM5_ESM.txt) - CCR enrichment in Pfam domains
    - [Supplementary Table 4](https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-018-0294-6/MediaObjects/41588_2018_294_MOESM6_ESM.txt) - Highly constrained CCRs not covered by missense depletion
    - [CCR Browser](https://s3.us-east-2.amazonaws.com/ccrs/ccr.html), [CCR BED files, autosomes](https://s3.us-east-2.amazonaws.com/ccrs/ccrs/ccrs.autosomes.v2.20180420.bed.gz), [CCR BED file, X chromosome](https://s3.us-east-2.amazonaws.com/ccrs/ccrs/ccrs.xchrom.v2.20180420.bed.gz)


- [Publicly available cancer GWAS](https://docs.google.com/document/d/13jmOtCNKiB8GPN0oTMJJOSVCQ_MUuO5bIG3zDKFHU0c/edit) by Peter Kraft. [Tweet](https://twitter.com/GENES_PK/status/1418543289968611330?s=20)

- [Open access to gnomAD data on multiple cloud providers](https://gnomad.broadinstitute.org/blog/2020-10-open-access-to-gnomad-data-on-multiple-cloud-providers/)

- `CAUSALdb` - curated summary statistics for GWASs, mapped to MeSH terms, Manhattan plot visualization. Download available. http://mulinlab.tmu.edu.cn/causaldb/index.html
    - Wang, Jianhua, Dandan Huang, Yao Zhou, Hongcheng Yao, Huanhuan Liu, Sinan Zhai, Chengwei Wu, et al. “CAUSALdb: A Database for Disease/Trait Causal Variants Identified Using Summary Statistics of Genome-Wide Association Studies.” Nucleic Acids Research, November 6, 2019, gkz1026. https://doi.org/10.1093/nar/gkz1026.

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

- [Awesome papers and projects about CNV and SV using NGS data](https://github.com/geocarvalho/sv-cnv-studies) - Relevant studies with Structual Variants and Copy Number Variants in NGS (Genome, Exome and Amplicon Sequencing) pipelines

- Benchmarking of 15 WGS-based structural variant callers (focusing on deletions), selected out of 55 (Supplementary Table 1). Gold standard - homozygous deletions in inbred mouse chromosomes. Lumpy, Manta give most optimal performance under various test strategies. [GitHub](https://github.com/Mangul-Lab-USC/benchmarking_SV), [Tweet](https://twitter.com/brent_p/status/1265803854785855490?s=20)
    - Sarwal, Varuni, Sebastian Niehus, Ram Ayyala, Sei Chang, Angela Lu, Nicholas Darci-Maher, Russell Littman, et al. “[A Comprehensive Benchmarking of WGS-Based Structural Variant Callers](https://doi.org/10.1101/2020.04.16.045120).” Preprint. April 18, 2020. 

- Benchmark of 10 CNV callers. LUMPY performs best overall, Canvas is good for high specificity, CNVnator and RDXplorer are good for high sensitivity. Table 1 summarizes functionality of each caller. Used the Database of Genomic Variants as a gold standard, call CNVs from NA12878 genome
    - Zhang, Le, Wanyu Bai, Na Yuan, and Zhenglin Du. “[**Comprehensively Benchmarking Applications for Detecting Copy Number Variation.**](https://doi.org/10.1371/journal.pcbi.1007069)” PLOS Computational Biology 15, no. 5 (May 28, 2019)

- Review of structural variant callers. De novo-based approaches (graph- or scaffold-based), short-read DNA-seq and RNA-seq (gene fusion) mapping, long-read (PacBio, Oxford Nanopore) mapping, multimethods approaches. SV calling from newer technologies, such as optical mapping, strand-seq, 10X Genomics linked reads, Hi-C. Brief description of tools, their performance, references to reviews. [Table 1](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1828-7/tables/1) - categorized list of tools, brief description and links.
    - Mahmoud, Medhat, Nastassia Gobet, Diana Ivette Cruz-Dávalos, Ninon Mounier, Christophe Dessimoz, and Fritz J. Sedlazeck. “[Structural Variant Calling: The Long and the Short of It](https://doi.org/10.1186/s13059-019-1828-7).” Genome Biology 20, no. 1 (December 2019): 246. 

- Benchmarking of seven SV callers (BreakDancer, Pindel, Delly, Lumpy, Manta, GRIDSS, SvABA) to detect different SV types, sizes, the effect of SV abundance and sequencing coverage, sequence similarity, biases (GC content and homopolymers), and mapping quality. Overview of read-pair, read-depth, split-read, and local assembly. SV types and definitions. Manta, Lumpy, of GRIDSS perform well. [Supplementary material](https://academic.oup.com/bib/article/22/3/bbaa056/5831479#supplementary-data) - code examples
    - Gong, Tingting, Vanessa M Hayes, and Eva K F Chan. “[Detection of Somatic Structural Variants from Short-Read next-Generation Sequencing Data](https://doi.org/10.1093/bib/bbaa056),” Briefings in Bioinformatics, May 2021

- [CNVpytor](https://github.com/abyzovlab/CNVpytor) - CNV/CNA detection from read depth and SNPs (allelic imbalance). Python implementation of CNVnator, faster parsing (pysam) smaller .pytor files (H5). Works on a cloud, in Jupyter notebooks. Visualization using JBrowse
    - Suvakov, Milovan, Arijit Panda, Colin Diesh, Ian Holmes, and Alexej Abyzov. “[CNVpytor: A Tool for CNV/CNA Detection and Analysis from Read Depth and Allele Imbalance in Whole Genome Sequencing](https://doi.org/10.1101/2021.01.27.428472
),” biorXiv, January 27, 2021

- [HATCHet](https://github.com/raphael-group/hatchet) (Holistic Allele-specific Tumor Copy-number Heterogeneity) is an algorithm that infers allele and clone-specific CNAs and WGDs jointly across multiple tumor samples from the same patient, and that leverages the relationships between clones in these samples.
    - Zaccaria, Simone, and Benjamin J. Raphael. “[Accurate Quantification of Copy-Number Aberrations and Whole-Genome Duplications in Multi-Sample Tumor Sequencing Data](https://doi.org/10.1101/496174).” BioRxiv, January 1, 2018

- [SV-Bay](https://github.com/BoevaLab/SV-Bay) - structural variant detection, Bayesian, corrects for GC-content and mappability. Uses both normal and abnormal reads, paired-end and depth information. Somatic variants if a normal sample is available. Detailed methods, stats. Compared with GASVPro, Lumpy, BreakDancer, DELLY on simulated (TGsim) and experimental neuroblastoma datasets. Improve sensitivity and specificity of SV detection, less false positives. Reasonably fast.
    - Iakovishina, Daria, Isabelle Janoueix-Lerosey, Emmanuel Barillot, Mireille Regnier, and Valentina Boeva. “[SV-Bay: Structural Variant Detection in Cancer Genomes Using a Bayesian Approach with Correction for GC-Content and Read Mappability](https://doi.org/10.1093/bioinformatics/btv751).” Bioinformatics, (April 1, 2016)

- [CNVkit](https://github.com/etal/cnvkit) - capturing CNVs in on-target and off-target genomic regions. Existing tools (CNVer, ExomeCNV, exomeCopy, CONTRA, CoNIFER, ExomeDepth, VarScan2, XHMM, ngCGH, EXCAVATOR, CANOES, PatternCNV, CODEX, Control-FREEC, cn.MOPS, cnvOffSeq, CopyWriteR). Account for GC content, mappability. Python 2.7 implementation. [Examples](https://github.com/etal/cnvkit-examples)
    - Talevich, Eric, A. Hunter Shain, Thomas Botton, and Boris C. Bastian. “[CNVkit: Genome-Wide Copy Number Detection and Visualization from Targeted DNA Sequencing](https://doi.org/10.1371/journal.pcbi.1004873).” PLOS Computational Biology, (April 21, 2016)

- [Manta](https://github.com/Illumina/manta) - SV detection in single- and tumor-normal samples. parallelized for within-sample performance. Fast, detects more variants of different types.
    - Chen, Xiaoyu, Ole Schulz-Trieglaff, Richard Shaw, Bret Barnes, Felix Schlesinger, Morten Källberg, Anthony J. Cox, Semyon Kruglyak, and Christopher T. Saunders. “[Manta: Rapid Detection of Structural Variants and Indels for Germline and Cancer Sequencing Applications](https://doi.org/10.1093/bioinformatics/btv710).” Bioinformatics, (April 15, 2016)

- [LUMPY-SV](https://github.com/arq5x/lumpy-sv/) - a general probabilistic framework for structural variant discovery. Integrates multiple signals - read-pair, split-read, read-depth and prior knowledge. Operates on multiple samples. 
    - Layer, Ryan M, Colby Chiang, Aaron R Quinlan, and Ira M Hall. “[**LUMPY: A Probabilistic Framework for Structural Variant Discovery.**](https://doi.org/10.1186/gb-2014-15-6-r84)” Genome Biology 15, no. 6 (2014): R84. 

- [DELLY](https://tobiasrausch.com/delly/) - detection of structural variants, such as CNVs, duplications, inversions, translocations. Paired-end and split-read analysis 
    - Rausch, Tobias, Thomas Zichner, Andreas Schlattl, Adrian M. Stütz, Vladimir Benes, and Jan O. Korbel. “[DELLY: Structural Variant Discovery by Integrated Paired-End and Split-Read Analysis](https://doi.org/10.1093/bioinformatics/bts378).” Bioinformatics, (September 15, 2012)

- [Control-FREEC](http://boevalab.com/FREEC/) - assess copy number and genotype information in whole genome and exome sequencing data. Corrects for contamination by normal cells and variable sample ploidy. With a matched normal sample, distinguishes somatic from germline events.
    - Boeva, Valentina, Tatiana Popova, Kevin Bleakley, Pierre Chiche, Julie Cappo, Gudrun Schleiermacher, Isabelle Janoueix-Lerosey, Olivier Delattre, and Emmanuel Barillot. “[Control-FREEC: A Tool for Assessing Copy Number and Allelic Content Using next-Generation Sequencing Data](https://doi.org/10.1093/bioinformatics/btr670).” Bioinformatics, (February 1, 2012)

- [CNVnator](https://github.com/abyzovlab/CNVnator) - CNV using read depth. Bin the genome, mean-shift technique to quantify CNVs. Poor agreement in methods detection. Misses retrotransposon-CNVs but detects >50% other regions missed by other methods. https://github.com/abyzovlab/CNVnator
    - Abyzov, Alexej, Alexander E. Urban, Michael Snyder, and Mark Gerstein. “[**CNVnator: An Approach to Discover, Genotype, and Characterize Typical and Atypical CNVs from Family and Population Genome Sequencing.**](https://doi.org/10.1101/gr.114876.110)” Genome Research 21, no. 6 (June 2011)

- [CREST](https://www.stjuderesearch.org/site/lab/zhang) - mapping somatic structural variations in cancer genomes ad base-pair resolution. realignment of soft-clipped subsequences.Compared with BreakDancer, Pindel. Tested on experimental and simulated data, hg18, some SVs experimentally validated. 
    - Wang, Jianmin, Charles G Mullighan, John Easton, Stefan Roberts, Sue L Heatley, Jing Ma, Michael C Rusch, et al. “[CREST Maps Somatic Structural Variation in Cancer Genomes with Base-Pair Resolution](https://doi.org/10.1038/nmeth.1628).” Nature Methods 8, no. 8 (August 2011)

- [DNAcopy](https://bioconductor.org/packages/DNAcopy/) - an R package for DNA copy number data analysis using Circular Binary Segmentation (CBS) algorithm. Input - CHG-like data, regions with comparative copy number values. Smoothing, outlier detection, segmentation using changepoint detection, visualization

- [copynumber](https://bioconductor.org/packages/copynumber/) - Segmentation of single- and multi-track copy number data by penalized least squares regression. Similar functionality to [DNAcopy](https://bioconductor.org/packages/DNAcopy/)

- [instasv](https://github.com/venyao/intansv) - Integrative analysis of structural variations called by different tools.

- [TITAN](http://compbio.bccrc.ca/software/titan/) - a tool for predicting subclonal copy number alterations (CNA) and loss of heterozygosity (LOH) from tumour whole genome sequencing data.

- [QDNASeq](https://bioconductor.org/packages/release/bioc/html/QDNAseq.html) - Quantitative DNA sequencing for chromosomal aberrations. 

- [samplot](https://github.com/ryanlayer/samplot) - Plot structural variant signals from many BAMs and CRAMs. 

- [smoove](https://github.com/brentp/smoove) - structural variant calling and genotyping with existing tools, but, smoothly. 

- [gnomAD-SV](https://gnomad.broadinstitute.org/) - structural variants from deep WGS, added to gnomAD browser. 14,891 genomes, average statistics of SVs in general population.
    - Collins, Ryan L., Harrison Brand, Konrad J. Karczewski, Xuefang Zhao, Jessica Alföldi, Amit V. Khera, Laurent C. Francioli, et al. “[An Open Resource of Structural Variation for Medical and Population Genetics](https://doi.org/10.1101/578674).” BioRxiv, March 14, 2019.

## Power

- [powerEQTL](https://bwhbioinfo.shinyapps.io/powerEQTL/) - power analysis for eQTL studies. Includes two models, one-way unbalanced ANOVA (categorical genotypes) and linear regression (additive counting genotypes). Applicable to bulk and scRNA-seq. [Shiny app](https://bwhbioinfo.shinyapps.io/powerEQTL/), [CRAN](https://CRAN.R-project.org/package=powerEQTL)
    - Dong, Xianjun, Xiaoqi Li, Tzuu-Wang Chang, Scott T Weiss, and Weiliang Qiu. “[PowerEQTL: An R Package and Shiny Application for Sample Size and Power Calculation of Bulk Tissue and Single-Cell EQTL Analysis](https://doi.org/10.1101/2020.12.15.422954),” bioRxiv, December 16, 2020.

- `GeneticsDesign` - GWAS power analysis, functions for designing genetics studies, https://www.bioconductor.org/packages/release/bioc/html/GeneticsDesign.html

- Genetic Association Study (GAS) Power Calculator - online tool to compute statistical power for large one-stage genetic association studies. The underlying method is derived from the CaTS power calculator for two-stage association studies (2006). http://csg.sph.umich.edu/abecasis/gas_power_calculator/index.html

- `SEQPower` - GWAS power analysis for case/control and quantitative studies, rare variants. Command-line, http://bioinformatics.org/spower/start
    - Wang, Gao T., Biao Li, Regie P. Lyn Santos-Cortez, Bo Peng, and Suzanne M. Leal. “Power Analysis and Sample Size Estimation for Sequence-Based Association Studies.” Bioinformatics 30, no. 16 (August 15, 2014): 2377–78. https://doi.org/10.1093/bioinformatics/btu296.


## Miscellaneous

 - [Awesome papers and projects about CNV and SV using NGS data](https://github.com/geocarvalho/sv-cnv-studies)

- [Atlas of Variant Age](https://human.genome.dating/) - age estimate of \~45M SNPs. Method - Genealogical Estimation of Variant Age (GEVA), performs similar or better to PSMC. https://human.genome.dating/
    - Albers, Patrick K., and Gil McVean. “Dating Genomic Variants and Shared Ancestry in Population-Scale Sequencing Data.” Edited by Nick H. Barton. PLOS Biology 18, no. 1 (January 17, 2020): e3000586. https://doi.org/10.1371/journal.pbio.3000586.

- Long-range sequencing and software for data processing. SMRT by Pacific Bioscience, nanopore-based by Oxford Nanopore, genome partitioning and barcoding by 10X Genomics, Hi-C based, BioNano optical mapping (Table 1). Table 2 - Bioinformatics tools for de novo genome assembly, SNP/CNV etc. variant detection, phasing, RNA-seq, methylation. Text describes each tool.
    - Sedlazeck, Fritz J., Hayan Lee, Charlotte A. Darby, and Michael C. Schatz. “[Piercing the Dark Matter: Bioinformatics of Long-Range Sequencing and Mapping](https://doi.org/10.1038/s41576-018-0003-4).” Nature Reviews Genetics 19, no. 6 (June 2018)

- GWAS tutorial. Quality control with PLINK, population stratification (MDS), association tests (binary and quantitative using PLINK), polygenic risk scores гыштп ЗКЫ. PLINK file formats. Box 1 - GWAS definitions.https://github.com/MareesAT/GWA_tutorial/
    - Marees, Andries T., Hilde de Kluiver, Sven Stringer, Florence Vorspan, Emmanuel Curis, Cynthia Marie-Claire, and Eske M. Derks. “A Tutorial on Conducting Genome-Wide Association Studies: Quality Control and Statistical Analysis.” International Journal of Methods in Psychiatric Research 27, no. 2 (2018): e1608. https://doi.org/10.1002/mpr.1608.

- DNA sequencing analysis notes from Ming Tang. https://github.com/crazyhottommy/DNA-seq-analysis

- `SNPhylo` - A pipeline to generate a phylogenetic tree from huge SNP data, https://github.com/thlee/SNPhylo, http://chibba.pgml.uga.edu/snphylo/

- Application for making ENCODE Blacklists, and links to canonical blacklists, https://github.com/Boyle-Lab/Blacklist
    - Blacklist citation: Amemiya, Haley M., Anshul Kundaje, and Alan P. Boyle. “The ENCODE Blacklist: Identification of Problematic Regions of the Genome.” Scientific Reports 9, no. 1 (December 2019): 9354. https://doi.org/10.1038/s41598-019-45839-z.


- HOT/XOT regions. The high occupancy target (HOT) and extreme occupancy target (XOT) regions in all contexts were downloaded through the ENCODE data portal at http://encode-ftp.s3.amazonaws.com/modENCODE_VS_ENCODE/Regulation/Human/hotRegions/maphot_hs_selection_reg_cx_simP05_all.bed and http://encode-ftp.s3.amazonaws.com/modENCODE_VS_ENCODE/Regulation/Human/hotRegions/maphot_hs_selection_reg_cx_simP01_all.bed (hg38 ?). Potential source

- `GEM` - mappability calculations for each genomic region, accounting for mismatches. Pre-calculated UCSC genome browser tracks for human and mouse. Mappability of genes, both protein-coding and non-protein coding. RPKUM - unique exons for quantifying gene expression. https://sourceforge.net/projects/gemlibrary/files/gem-library/
    - Derrien, Thomas, Jordi Estellé, Santiago Marco Sola, David G. Knowles, Emanuele Raineri, Roderic Guigó, and Paolo Ribeca. “Fast Computation and Applications of Genome Mappability.” PloS One 7, no. 1 (2012): e30377. https://doi.org/10.1371/journal.pone.0030377.

- `refGenie` - reference genome manager. http://refgenie.databio.org/en/latest/

- `genomepy` - Download genomes the easy way. https://github.com/simonvh/genomepy

- Sample swap check. https://github.com/parklab/NGSCheckMate, https://github.com/brentp/somalier
