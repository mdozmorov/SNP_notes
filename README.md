# SNP-related notes

[![MIT License](https://img.shields.io/apm/l/atomic-design-ui.svg?)](https://github.com/tterb/atomic-design-ui/blob/master/LICENSEs) [![PR's Welcome](https://img.shields.io/badge/PRs-welcome-brightgreen.svg?style=flat)](http://makeapullrequest.com) 

SNP- and genome variation-related tools and genomics data analysis resources. Please, [contribute and get in touch](CONTRIBUTING.md)! See [MDmisc notes](https://github.com/mdozmorov/MDmisc_notes) for other programming and genomics-related notes.

# Table of content

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->


- [Variant calling pipelines](#variant-calling-pipelines)
  - [Preprocessing tools](#preprocessing-tools)
  - [Depth](#depth)
  - [Genome comparison](#genome-comparison)
- [SNP callers](#snp-callers)
  - [Deep learning SNP callers](#deep-learning-snp-callers)
- [SNP annotations](#snp-annotations)
- [SNP signatures](#snp-signatures)
- [SNP pathogenicity scores](#snp-pathogenicity-scores)
- [SNP visualization, clustering](#snp-visualization-clustering)
- [SNP, GWAS databases](#snp--gwas-databases)
  - [GWAS pipelines](#gwas-pipelines)
  - [Ancestry](#ancestry)
  - [eQTLs](#eqtls)
  - [Polygenic risk score](#polygenic-risk-score)
  - [Regulatory](#regulatory)
- [InDels](#indels)
- [CNV, SV](#cnv-sv)
- [Power](#power)
- [Miscellaneous](#miscellaneous)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## Variant calling pipelines

- [SAMtools, BCFtools, and HTSlib](https://www.htslib.org/) development history, functionality. Work directly with SAM, BAM, or CRAM formats. Improved functionality, e.g., indexing files as they written, speed improvement. [HTSlib landing page](https://www.htslib.org/) with download, workflows, documentation, support links. [HTSlib manual pages](https://www.htslib.org/doc/#manual-pages), [SAMtools manual](https://www.htslib.org/doc/samtools.html), [BCFtools manual](https://www.htslib.org/doc/bcftools.html). [Supplementary data](http://gigadb.org/dataset/100866), Table S1: SAMtools commands; Table S2: BCFtools commands. <details>
    <summary>Paper</summary>
    Danecek, Petr, James K Bonfield, Jennifer Liddle, John Marshall, Valeriu Ohan, Martin O Pollard, Andrew Whitwham, et al. “Twelve Years of SAMtools and BCFtools.” GigaScience 10, no. 2 (January 29, 2021): giab008. https://doi.org/10.1093/gigascience/giab008.
</details>

- Koboldt, Daniel C. “[Best Practices for Variant Calling in Clinical Sequencing](https://doi.org/10.1186/s13073-020-00791-w).” Genome Medicine, (December 2020) - Introduction in genomic variant calling, panel/exome/whole genome sequencing technologies (Table 1), preprocessing, analysis (SNVs/indels, mutations, CNVs, SVs, gene fusions, Table 2), gold standard datasets (GIAB), best practices, filtering for each type of genomic variant.

- [crg2](https://github.com/ccmbioinfo/crg2) is a research pipeline aimed at discovering clinically relevant variants (SNVs, SVs) in whole genome sequencing data. Snakemake and Conda, runs on PBS. Python, shell. Conda, local installation.

- [DNAscan](https://github.com/KHP-Informatics/DNAscan) - a pipeline for DNA-seq analysis to call SNPs, indels, SVs, repeat expansions and viral genetic material. Four parts - alignment (HISAT2, BWA mem), analysis (Freebayes, GATC HC, Manta, Expansion Hunter), annotation (Annovar), report generation. Python, shell. Docker, local installation. <details>
    <summary>Paper</summary>
    Iacoangeli, A., A. Al Khleifat, W. Sproviero, A. Shatunov, A. R. Jones, S. L. Morgan, A. Pittman, R. J. Dobson, S. J. Newhouse, and A. Al-Chalabi. “DNAscan: Personal Computer Compatible NGS Analysis, Annotation and Visualisation.” BMC Bioinformatics 20, no. 1 (December 2019): 213. https://doi.org/10.1186/s12859-019-2791-8.
</details>

- [MapCaller](https://github.com/hsinnan75/MapCaller) - map reads and identify genomic variants, indels, inversions, translocations. Based on KART method. Outperforms GATK, Freebayes, Mpileup. C, C++. Conda, local installation. <details>
    <summary>Paper</summary>
    Lin, Hsin-Nan, and Wen-Lian Hsu. “MapCaller - An Integrated and Efficient Tool for Short-Read Mapping and Variant Calling Using High-Throughput Sequenced Data.” Preprint. Bioinformatics, September 26, 2019. https://doi.org/10.1101/783605.
</details>

- [TOPMed Variant Calling Pipeline](https://github.com/statgen/topmed_variant_calling). C++, Python. Runs on Google cloud.


### Preprocessing tools

- [genomics_general](https://github.com/simonhmartin/genomics_general) - This is a collection of scripts for a range of genomic data processing and analysis. Processing VCF files; Diversity and divergence analyses in sliding windows; Distance matrix; ABBA-BABA statistics in sliding windows; Trees for sliding windows. Python.

- [alleleCount](https://github.com/cancerit/alleleCount) - Takes a file of locations and a [cr|b]am file and generates a count of coverage of each allele [ACGT] at that location. C, Perl, shell. Docker.

- [BAMsurgeon](https://github.com/adamewing/bamsurgeon) tools for adding mutations to existing .bam files, used for testing mutation callers. Python, shell.

- [QTLtools](https://qtltools.github.io/qtltools/) - A complete tool set for molecular QTL discovery and analysis. Has `mbv` mode to match BAM to VCF, PCA, functional enrichment, colocalization analyses, more. C++.

- [PyVCF](https://pyvcf.readthedocs.io) - A Variant Call Format Parser for Python. Has `vcf_melt` tool to reformat a VCF into long format. Python.

- [Rust-bio-tools](https://github.com/rust-bio/rust-bio-tools) - VCF matching, conversion to text, report, FASTQ split/filter, BAM depth, merging. Rust.

- [SURVIVOR](https://github.com/fritzsedlazeck/SURVIVOR/) - Toolset for SV simulation, comparison and filtering. C++.

- [VariantQC](https://github.com/BimberLab/DISCVRSeq/) - VCF quality control tool, part of [DISCVRseq](https://github.com/BimberLab/DISCVRSeq) toolkit. Uses GATK4 engine. Java wrapper of GATK's VariantEval tool. Input - VCF file and an indexed genome FASTA file. Output - MultiQC-templated report. Java. <details>
    <summary>Paper</summary>
    Yan, Melissa Y, Betsy Ferguson, and Benjamin N Bimber. “VariantQC: A Visual Quality Control Report for Variant Evaluation.” Edited by Jonathan Wren. Bioinformatics, July 16, 2019, btz560. https://doi.org/10.1093/bioinformatics/btz560.
</details>

- [vcf-diff](https://github.com/dfornika/vcf-diff) - Check for differences between two vcf files

- [vcf-split](https://github.com/auerlab/vcf-split) - Split combined-sample VCF stream into single-sample VCF files. C, shell.

- [vcfR](https://cran.r-project.org/web/packages/vcfR/index.html) - Manipulate and Visualize VCF Data. R.

- [vcf2tsvpy](https://github.com/sigven/vcf2tsvpy) - Genomic VCF to tab-separated values, by [Sigve Nakken](https://github.com/sigven). Python.

### Depth

- [bamCoverage](https://deeptools.readthedocs.io/en/latest/content/tools/bamCoverage.html) - BAM to bigWig conversion, part of [Deeptools](https://deeptools.readthedocs.io/). Python, shell.

- [BAMscale](https://github.com/ncbi/BAMscale) - BAMscale is a one-step tool for either 1) quantifying and normalizing the coverage of peaks or 2) generated scaled BigWig files for easy visualization of commonly used DNA-seq capture based methods. C, R.

- [histoneSig](https://github.com/semibah/histonesig) - R package for working with genome files as continuous representations or "signals". R.

- [goleft](https://github.com/brentp/goleft) - a collection of bioinformatics tools distributed under MIT license in a single static binary. covstats - estimate coverage and insert-size statistics on bams by sampling; depth - parallelize calls to samtools in user-defined windows; depthwed - matricize output from depth to n-sites * n-samples; indexcov - quick coverage estimate using only the bam index; indexsplit - generate regions of even data across a cohort (for parallelization); samplename - report samplename(s) from a bam's SM tag. Go.

- [mosdepth](https://github.com/brentp/mosdepth) - fast BAM/CRAM depth calculation for WGS, exome, or targetted sequencing. Nim, shell.

- [perbase](https://github.com/sstadick/perbase) - Per-base per-nucleotide depth analysis. Written in Rust, command line.

### Genome comparison

- [Ragout2](https://github.com/fenderglass/Ragout) - a reference-assisted assembly tool that works for large and complex genomes. Combines Cactus, a multiple whole-genome aligner, with a new iterative graph simplification algorithm that produces hierarchical synteny blocks on multiple scales. Utilizes information from multple genomes. Infers the evolutionary relationships between the genomes and builds the final assemblies using a genome rearrangement approach. Detects synteny blocks, resolves repeats, chimeric sequences. Accuracy comparable with assemblies from BioNano maps, Hi-C, BAC clones, FISH. Benchmarked against RACA on simulated and experimental (16 mouse genomes) data. See also [maf2synteny](https://github.com/fenderglass/maf2synteny) <details>
    <summary>Paper</summary>
    Kolmogorov, Mikhail, Joel Armstrong, Brian J. Raney, Ian Streeter, Matthew Dunn, Fengtang Yang, Duncan Odom, et al. “Chromosome Assembly of Large and Complex Genomes Using Multiple References.” Genome Research 28, no. 11 (November 2018): 1720–32. https://doi.org/10.1101/gr.236273.118.
</details>

- [Sibelia](https://github.com/bioinf/Sibelia) - syntheny finder in multiple microbial genomes using iterative de Bruijn graphs. Represents synteny blocks in a hierarchy. Outperforms Mugsy, Multiz, Mauve on synthetic and experimental data. [Webserver](http://bioinf.spbau.ru/sibelia). <details>
    <summary>Paper</summary>
    Minkin, Ilya, Anand Patel, Mikhail Kolmogorov, Nikolay Vyahhi, and Son Pham. “Sibelia: A Scalable and Comprehensive Synteny Block Generation Tool for Closely Related Microbial Genomes.” arXiv, July 30, 2013. http://arxiv.org/abs/1307.7941.
</details>

- [syntenyPlotteR](https://github.com/marta-fb/syntenyPlotteR) - R package to draw syntenic plots in 3 different styles (highway, inferCARs, linear pairwise style)

- [LASTZ](https://github.com/lastz/lastz) - pairwise DNA sequence aligner and plotter. [Documentation](https://lastz.github.io/lastz/).



## SNP callers

- Liu, Fenglin, Yuanyuan Zhang, Lei Zhang, Ziyi Li, Qiao Fang, Ranran Gao, and Zemin Zhang. "[Systematic comparative analysis of single-nucleotide variant detection methods from single-cell RNA sequencing data](https://doi.org/10.1186/s13059-019-1863-4)." Genome Biology, (19 November 2019) - comparison of seven tools for SNP detection in scRNA-seq data. SAMtools, Strelka2, FreeBayes, and CTAT are best.

- Xu, Chang. “[A Review of Somatic Single Nucleotide Variant Calling Algorithms for Next-Generation Sequencing Data](https://doi.org/10.1016/j.csbj.2018.01.003).” Computational and Structural Biotechnology Journal, (6 February 2018) - Overview of 46 somatic Single Nucleotide Variant (SNV) caller tools. Pre-processing, variant evaluation, and post-filtering steps. Four categories of algorithms, description of each, and the corresponding tools: matched tumor-normal (position-, haplotype-, mathine learning-based methods, [Table 1](https://www.sciencedirect.com/science/article/pii/S2001037017300946?via%3Dihub#t0005)), single-sample ([Table 2](https://www.sciencedirect.com/science/article/pii/S2001037017300946?via%3Dihub#t0010), some offer somatic-germline classification), UMI-based (UMI technology, [Figure 1](https://www.sciencedirect.com/science/article/pii/S2001037017300946?via%3Dihub#f0005), [Table 3](https://www.sciencedirect.com/science/article/pii/S2001037017300946?via%3Dihub#t0015)), and RNA-seq (Technology, issues, [Table 4](https://www.sciencedirect.com/science/article/pii/S2001037017300946?via%3Dihub#t0020)) variant calling. Benchmarking using tools for generating synthetic reads, spike-ins, GiAB, melanoma-normal samples, performance evaluation metrics. Issues in representing complex variants and tools for variant normalization. Deep neural network-based algorithms perform best.

- [ABRA](https://github.com/mozack/abra) - Assembly Based ReAligner. It uses localized assembly (de Bruijn graph) and global realignment to align reads more accurately, thus improving downstream analysis (detection of indels and complex variants in particular). For exome sequencing, requires targets in BED file. Extended with Cadabra, somatic indel variant caller. Java, C++. <details>
    <summary>Paper</summary>
    Mose, Lisle E., Matthew D. Wilkerson, D. Neil Hayes, Charles M. Perou, and Joel S. Parker. “ABRA: Improved Coding Indel Detection via Assembly-Based Realignment.” Bioinformatics 30, no. 19 (October 2014): 2813–15. https://doi.org/10.1093/bioinformatics/btu376.
</details>

- [GATK on Biowulf](https://hpc.nih.gov/training/gatk_tutorial/) - A practical introduction to GATK 4 on Biowulf (NIH HPC), by Qi Yu and Wolfgang Resch. Optimized scripts for HPC environment, [GitHub](https://github.com/DylanLawless/gatk4_tutorial)

- [HipSTR](https://hipstr-tool.github.io/HipSTR/) - short tandem repeat caller, resolving haplotypes. Learns locus-specific PCR stutter models using an EM algorithm, other methods so obtain more robust STR genotypes. Works on Illumina data. Input - BAM files, BED coordinates of pre-defined STRs, FASTA genome. Output - VCF of the STR genotypes. Visualization of genotypes. Works for a variety of data amount/quality. Outperforms lobstr, RepeatSeq, GATK HaplotypeCaller, Platypus, freebayes, SAMtools. [Tutorial](https://hipstr-tool.github.io/HipSTR-tutorial/). C, C++. <details>
    <summary>Paper</summary>
    Willems, Thomas, Dina Zielinski, Jie Yuan, Assaf Gordon, Melissa Gymrek, and Yaniv Erlich. “Genome-Wide Profiling of Heritable and de Novo STR Variations.” Nature Methods 14, no. 6 (June 2017): 590–92. https://doi.org/10.1038/nmeth.4267.
</details>

- [RNA-MuTect](https://zenodo.org/record/2620062#.XPz9e29KhQI) method to detect somatic mutations from tumor DNA - matched normal RNA samples. Applied to TCGA, GTeX data. Most mutated are sun-exposed skin, esophagus mucosa, lung. Number of mutations in RNA is 5-fold larger than in DNA, mutations depend on coverage. Only half of DNA mutations have sufficient coverage in RNA. Filtering using 1) two aligners, 2) removal of errors based on a site-specific error model, 3) removal of RNA editing sites. Java, shell. <details>
    <summary>Paper</summary>
    Yizhak, Keren, François Aguet, Jaegil Kim, Julian M Hess, Kirsten Kübler, Jonna Grimsby, Ruslana Frazer, et al. “RNA Sequence Analysis Reveals Macroscopic Somatic Clonal Expansion across Normal Tissues.” HUMAN GENETICS, 2019, 11. https://doi.org/10.1126/science.aaw0726
</details>

- [Samovar](https://github.com/cdarby/samovar) - mosaic SNV calling from single samples sequenced using WGS 10X linked reads technology. Intro into mosaic SNVs, linked reads. Outperforms MuTect2, MosaicHunter. Python, shell. <details>
    <summary>Paper</summary>
    Darby, Charlotte A, James R Fitch, Patrick J Brennan, Benjamin J Kelly, Natalie Bir, Vincent Magrini, Jeffrey Leonard, et al. “Samovar: Single-Sample Mosaic SNV Calling with Linked Reads.” Preprint. Genomics, February 25, 2019. https://doi.org/10.1101/560532.
</details>

- [slivar](https://github.com/brentp/slivar) - variant expressions, annotation, and filtering, by Brent Pedersen. expr - filter and/or annotate with INFO, trio, sample, group expressions; make-gnotate - make a compressed zip file of annotations for use by slivar; compound-hets - true compound hets using phase-by-inheritance within gene annotations. Nim, binary.

- [Strelka2](https://github.com/Illumina/strelka) - a fast and accurate small variant caller optimized for analysis of germline variation in small cohorts and somatic variation in tumor/normal sample pairs. Tiered haplotype model, selected according to the locus properties (Methods). Multiple improvements, speed-ups. Outperforms GATK4, FreeBayes on SNPs and indel ([PrecisionFDA](https://precision.fda.gov/), [GIAB](https://github.com/genome-in-a-bottle) datasets, in silico mixtures). C++, Python. <details>
    <summary>Paper</summary>
    Kim, Sangtae, Konrad Scheffler, Aaron L. Halpern, Mitchell A. Bekritsky, Eunho Noh, Morten Källberg, Xiaoyu Chen, et al. “Strelka2: Fast and Accurate Calling of Germline and Somatic Variants.” Nature Methods 15, no. 8 (August 2018): 591–94. https://doi.org/10.1038/s41592-018-0051-x.
</details>

- [UNCeqR](https://github.com/mwilkers/unceqr) - Somation mutation prediction by RNA and DNA integration. C.

- [UNMASC](https://github.com/pllittle/UNMASC) - R package for tumor-only variant calling. Filter and annotate tumor-only variant calls through the integration of public database annotations (pools of unmatched normals), clustering, and segmentation (data-driven detection of hard to map regions) to provide the user with a clear characterization of each variant when called against a set of unmatched normal controls (10 sufficient). Data and database at [dbGAP](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs001713.v1.p1). R, C++. <details>
    <summary>Paper</summary>
    Little, Paul, Heejoon Jo, Alan Hoyle, Angela Mazul, Xiaobei Zhao, Ashley H Salazar, Douglas Farquhar, et al. “UNMASC: Tumor-Only Variant Calling with Unmatched Normal Controls.” NAR Cancer 3, no. 4 (October 4, 2021): zcab040. https://doi.org/10.1093/narcan/zcab040.
</details>

- [Varlociraptor](https://varlociraptor.github.io) - a unifying statistical model allowing for the detection of SNVs, MNVs, InDels, Inversions, Duplications, Breakends. Explicitly controls FDR. Defines variant calling grammar, allowing for defining all types of variants. Tumor-normal comparison allows for classifying variants into germline, somatic variants. Evaluated on simulated (artificial clones) and real (Venter's genome) data. Outperforms six other variant callers. Can call variants in single samples, in RNA-seq data, FFPE attifact detection. [Tweet](https://twitter.com/johanneskoester/status/1331218560367128577?s=20). [dna-seq-varlociraptor](dna-seq-varlociraptor) - A Snakemake workflow for calling small and structural variants under any kind of scenario (tumor/normal, tumor/normal/relapse, germline, pedigree, populations) via the unified statistical model of Varlociraptor. Rust, Conda installation. <details>
    <summary>Paper</summary>
    Köster, Johannes, Louis J. Dijkstra, Tobias Marschall, and Alexander Schönhuth. “[Varlociraptor: Enhancing Sensitivity and Controlling False Discovery Rate in Somatic Indel Discovery](https://doi.org/10.1186/s13059-020-01993-6).” Genome Biology 21, no. 1 (December 2020)
</details>

### Deep learning SNP callers

- [DeepSNV](https://bioconductor.org/packages/deepSNV/) - an R package for detecting clonal and subclonal variants from deeply sequenced clonally heterogeneous cancer data. The deepSNV algorithm is used for a comparative setup with a control experiment of the same loci and uses a beta-binomial model and a likelihood ratio test to discriminate sequencing errors and subclonal SNVs. The shearwater algorithm computes a Bayes classifier based on a beta-binomial model for variant calling with multiple samples for precisely estimating model parameters - such as local error rates and dispersion - and prior knowledge, e.g. from variation data bases such as COSMIC. R. <details>
    <summary>Paper</summary>
    Gerstung, M., E. Papaemmanuil, and P. J. Campbell. “[Subclonal Variant Calling with Multiple Samples and Prior Knowledge](https://doi.org/10.1093/bioinformatics/btt750).” Bioinformatics, (May 1, 2014)
</details>

- [DeepVariant](https://github.com/google/deepvariant) is an analysis pipeline that uses a deep neural network to call genetic variants from next-generation DNA sequencing data. CNN trained to detect SNPs. Inception v2 architecture. Outperforms GATK and others in FDA challenge, GIB data. Docker/Singularity CPU and GPU versions, local and cloud compatible. [Tweet 1](https://twitter.com/acarroll_ATG/status/1194786628759216128?s=20), [Tweet 2](https://twitter.com/brent_p/status/1194792729907093505?s=20), [GitHub repo](https://github.com/google/deepvariant/). Python, C++. Docker. <details>
    <summary>Paper</summary>
    Poplin, Ryan, Pi-Chuan Chang, David Alexander, Scott Schwartz, Thomas Colthurst, Alexander Ku, Dan Newburger, et al. “[Creating a Universal SNP and Small Indel Variant Caller with Deep Neural Networks,](https://doi.org/10.1101/092890)” March 20, 2018
</details>

- [NeuSomatic](https://github.com/bioinform/neusomatic) - Convolutional neural network (9 layers) for somatic mutations calling. Reads within 7bp window around candidate SNP are extracted, realigned, summarized into matrices for tumor-normal samples,  used for classifying mutation type, length, position. Tested on GiB samples and on DREAM datasets. Comparison with other SNP callers. Python, C++. Conda. <details>
    <summary>Paper</summary>   
    Sahraeian, Sayed Mohammad Ebrahim, Ruolin Liu, Bayo Lau, Marghoob Mohiyuddin, and Hugo Y. K. Lam. “Deep Convolutional Neural Networks for Accurate Somatic Mutation Detection,” September 4, 2018. https://doi.org/10.1101/393801.
</details>

## SNP annotations

- [awesome-cancer-variant-databases](https://github.com/seandavi/awesome-cancer-variant-databases) - A community-maintained repository of cancer clinical knowledge bases and databases focused on cancer and normal variants, by Sean Davis.

- Various genome annotations, [Source: ConsHMM Data availability section](https://www.nature.com/articles/s42003-019-0488-1#data-availability): [25-state chromatin state annotations](http://compbio.mit.edu/roadmap); [CADD score v1.0](http://krishna.gs.washington.edu/download/CADD/v1.0/whole_genome_SNVs.tsv.gz); [CADD score v1.4](http://krishna.gs.washington.edu/download/CADD/v1.4/GRCh37/whole_genome_SNVs.tsv.gz); [CDTS score](http://www.hli-opendata.com/noncoding/coord_CDTS_percentile_N7794unrelated.txt.gz), [another link](http://www.hli-opendata.com/noncoding/SNVusedForCDTScomputation_N7794unrelated_allelicFrequency0.001truncated.txt.gz); [CNEEs](http://www.stanford.edu/~lowec/data/threePeriods/hg19cnee.bed.gz); [DANN score](https://cbcl.ics.uci.edu/public_data/DANN/data/); [EIGEN and Eigen-PC score](https://xioniti01.u.hpc.mssm.edu/v1.1/); [ENCODE DHS](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwDnase/); [FATHMM-XF score](http://fathmm.biocompute.org.uk/fathmm-xf/); [FIRE score](https://sites.google.com/site/fireregulatoryvariation/); [fitCons score](http://compgen.cshl.edu/fitCons/0downloads/tracks/i6/scores/); [FunSeq2 score](http://org.gersteinlab.funseq.s3-website-us-east-1.amazonaws.com/funseq2.1.2/hg19_NCscore_funseq216.tsv.bgz); [GENCODE v19](https://www.gencodegenes.org/releases/19.html); [GERP++ scores and constrained element calls](http://mendel.stanford.edu/SidowLab/downloads/gerp/); [GWAS catalog variants](https://www.ebi.ac.uk/gwas/); [LINSIGHT score](http://compgen.cshl.edu/~yihuang/tracks/LINSIGHT.bw); [Motif instances and background](http://compbio.mit.edu/encode-motifs/); [REMM score](https://zenodo.org/record/1197579/files/ReMM.v0.3.1.tsv.gz); [Roadmap Epigenomics DHS](http://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/); [SiPhy-omega and SiPhy-pi constrained element calls (hg19 liftOver)](https://www.broadinstitute.org/mammals-models/29-mammals-project-supplementary-info)

- [ANNOVAR](https://annovar.openbioinformatics.org/en/latest/user-guide/download/) - functional annotation of genetic variants (hg18, hg19, hg38, mouse, worm, fly, yeast, etc.). Gene-based, region-based, filter-based annotations. [wANNOVAR](http://wannovar.wglab.org/) - web interface. <details>
    <summary>Paper</summary>
    Wang K, Li M, Hakonarson H. ANNOVAR: functional annotation of genetic variants from high-throughput sequencing data. Nucleic acids research. 2010 Sep 1;38(16):e164-. https://doi.org/10.1093/nar/gkq603
</details>

- [atSNP](http://atsnp.biostat.wisc.edu/) - search for effects of SNPs on transcription factor binding. DB of 37 billion variant-motif pairs. Search by SNP IDs, window of SNPs, genomic location, gene, transcription factor. <details>
    <summary>Paper</summary>
    Shin, Sunyoung, Rebecca Hudson, Christopher Harrison, Mark Craven, and Sündüz Keleş. “AtSNP Search: A Web Resource for Statistically Evaluating Influence of Human Genetic Variation on Transcription Factor Binding.” Edited by John Hancock. Bioinformatics, December 8, 2018. https://doi.org/10.1093/bioinformatics/bty1010.
</details>

- [Cancer Genome Interpreter](https://www.cancergenomeinterpreter.org/) - annotation and interpretation of cancer genome variants. A compiled [Catalog of Cancer Genes](https://www.cancergenomeinterpreter.org/genes) involved in the onset and progression of different types of cancer, obtained via different methods and from different sources. 837 genes with evidence of a tumorigenic role in 193 different cancer types. C[atalog of Validated Oncogenic Mutations](https://www.cancergenomeinterpreter.org/mutations) - protein-affecting mutations (PAMs) that occur in genes of the Catalog of Cancer Genes. The [Cancer Biomarkers database](https://www.cancergenomeinterpreter.org/biomarkers), an expanded collection of genomic biomarkers of anti-cancer drug response on 1624 genomic biomarkers of response (sensitivity, resistance, or toxicity) to 310 drugs across 130 types of cancer. The [Cancer Bioactivities database](https://www.cancergenomeinterpreter.org/bioactivities) of 20,243 chemical compound–protein product interactions that may support novel research applications. <details>
    <summary>Paper</summary>
    Tamborero, David, Carlota Rubio-Perez, Jordi Deu-Pons, Michael P. Schroeder, Ana Vivancos, Ana Rovira, Ignasi Tusquets et al. "Cancer Genome Interpreter annotates the biological and clinical relevance of tumor alterations." Genome medicine, (28 March 2018) https://doi.org/10.1186/s13073-018-0531-8
</details>

- [CIViC](https://civicdb.org/welcome) - Clinical Interpretation of Variants in Cancer database, community curated resource. [Data releases](https://civicdb.org/releases/main) Supports variant interpretation guidelines, interfaces with other databases. [CIViCmine](http://bionlp.bcgsc.ca/civicmine/) - Text mined biomarkers in cancer for curation into the CIViC database. [Documentation](https://civic.readthedocs.io/en/latest/model/variants/name.html), [API](https://griffithlab.github.io/civic-v2/), [Youtube](https://www.youtube.com/channel/UCtZ14M7ZNmMDivSJzaHO8BQ),  [GitHub](https://github.com/griffithlab/civic-v2) <details>
    <summary>Paper</summary>
    Krysiak, Kilannin, Arpad M Danos, Jason Saliba, Joshua F McMichael, Adam C Coffman, Susanna Kiwala, Erica K Barnell, et al. “CIViCdb 2022: Evolution of an Open-Access Cancer Variant Interpretation Knowledgebase.” Nucleic Acids Research, November 14, 2022, gkac979. https://doi.org/10.1093/nar/gkac979.
</details>

- [clinvar](https://github.com/macarthur-lab/clinvar) -  tools to convert ClinVar data into a tab-delimited flat file, and also provides that resulting tab-delimited flat file.

- [ConsHMM](https://github.com/ernstlab/ConsHMM) - genome segmentation into 100 conservation states based on a 100 species DNA sequence alignment. Hidden Markov Model, extension of ChromHMM. Compared with 12 other scores (CADD, CDTS, DANN, Eigen, Eigen-PC, FATHMM-XF, FIRE, fitCons, GERP++, PhastCons, PhyloP, REMM, also LINSIGHT, FunSeq2). Clustering, GO enrichment reveals distinct functionality. hg19 single nucleotide scores. [Download](https://figshare.com/articles/ConsHMM_100-state_Segmentation_of_hg19_Human_Genome/8162036/1). Java, Python. Conda. <details>
    <summary>Paper</summary>
    Arneson, Adriana, and Jason Ernst. “Systematic Discovery of Conservation States for Single-Nucleotide Annotation of the Human Genome.” Communications Biology 2, no. 1 (December 2019): 248. https://doi.org/10.1038/s42003-019-0488-1.
</details>

- [DoCM](http://www.docm.info) - database of manually curated cancer SNPs from 876 publications and databases. [API](http://www.docm.info/api) to get the data, [downloads](http://www.docm.info/downloads) in other formats. <details>
    <summary>Paper</summary>
    Ainscough, Benjamin J., Malachi Griffith, Adam C. Coffman, Alex H. Wagner, Jason Kunisaki, Mayank Nk Choudhary, Joshua F. McMichael, et al. “DoCM: A Database of Curated Mutations in Cancer.” Nature Methods 13, no. 10 (September 29, 2016): 806–7. https://doi.org/10.1038/nmeth.4000.
</details>

- [echtvar](https://github.com/brentp/echtvar) - fast simultaneous annotation and filtering genetic variants using 32-bit integer encoding of population variants and annotation fields (gnomAD, CADD) into a compressed archive. Input: VCF/BCF files, chromosome, and a JSON5 configuration. Faster/more efficient than slivar, vcfanno, bcftools, snpSift. Examples of filtering somatic variants with dbNSFP, whole genome recessive germline variants and somatic variants. Rust. <details>
    <summary>Paper</summary>
    Pedersen, Brent S, and Jeroen de Ridder. “Echtvar: Compressed Variant Representation for Rapid Annotation and FIltering of SNPs and Indels,” n.d.
</details>

- [FENRIR](https://fenrir.flatironinstitute.org/) - Tissue-specific enhancer functional networks for associating SNPs in distal regulatory regions to disease. Integrates thousands of disparate epigenetic and functional genomics datasets to infer tissue-specific functional relationships between en- hancers for 140 diverse human tissues and cell types, providing a regulatory-region-centric approach to systematically identify disease-associated enhancers. A Bayesian integration model that creates tissue-specific enhancer functional networks and, second, a network-based machine-learning model (elastic net) that prioritizes disease-enhancer association in a tissue-specific manner. Outperforms EMERGE, applied to autism variants. [GitHub](https://github.com/xichensf/fenrir). <details>
    <summary>Paper</summary>
    Chen, Xi, Jian Zhou, Ran Zhang, Aaron K. Wong, Christopher Y. Park, Chandra L. Theesfeld, and Olga G. Troyanskaya. “Tissue-Specific Enhancer Functional Networks for Associating Distal Regulatory Regions to Disease.” Cell Systems 12, no. 4 (April 2021): 353-362.e6. https://doi.org/10.1016/j.cels.2021.02.002.
</details>

- [FUMA](https://fuma.ctglab.nl/) - Functional Mapping and Annotation of GWAS using 43 scRNA-seq datasets (human and mouse). MAGMA Cell type-specific enrichment. Applied to 26 GWAS disorders. [Download](https://github.com/Kyoko-wtnb/FUMA_scRNA_data) processed data. <details>
    <summary>Paper</summary>
    Watanabe, Kyoko, Maša Umićević Mirkov, Christiaan A. de Leeuw, Martijn P. van den Heuvel, and Danielle Posthuma. “Genetic Mapping of Cell Type Specificity for Complex Traits.” Nature Communications 10, no. 1 (December 2019): 3222. https://doi.org/10.1038/s41467-019-11181-1.
</details>

- [ONCOTATOR](https://www.broadinstitute.org/oncotator/) - Cancer-oriented annotation of SNPs. 14 data sources. Goal - to identify somatic SNPs, eliminate germline. <details>
    <summary>Paper</summary>
    Ramos, Alex H., Lee Lichtenstein, Manaswi Gupta, Michael S. Lawrence, Trevor J. Pugh, Gordon Saksena, Matthew Meyerson, and Gad Getz. “Oncotator: Cancer Variant Annotation Tool.” Human Mutation 36, no. 4 (April 2015): E2423-2429. https://doi.org/10.1002/humu.22771.
</details>

- [OpenCRAVAT](https://opencravat.org/index.html) - Open Custom Ranked Analysis of Variants Toolkit, extension of the Cancer-Related Analysis of Variants Toolkit (CRAVAT). Command-line and GUI interface, extensive resource catalog. Input: VCF, annotated VCF, basic tabular file format, dbSNP identifiers, 23andMe, and Ancestry.com. Output: viewer with 4 tabs: Summary, Variant, Gene, and Filter. Python. <details>
    <summary>Paper</summary>
    Pagel, Kymberleigh A., Rick Kim, Kyle Moad, Ben Busby, Lily Zheng, Collin Tokheim, Michael Ryan, and Rachel Karchin. "Integrated informatics analysis of cancer-related variants." JCO clinical cancer informatics 4 (March 4, 2020): 310-317. https://doi.org/10.1200/cci.19.00132
</details>

- [RegTools](https://regtools.readthedocs.io/en/latest/) - tools o integrate somatic variants from genomic data with splice junctions from bulk or single cell transcriptomic data to identify cis-acting variants that may cause aberrant splicing. Three sub-modules: a _variants_ module to annotate genomic variant calls for their potential splicing relevance, a _junctions_ module to analyze aligned RNA-seq data to extract and annotate splice junctions, and a _cis-splice-effects_ module that associates these variants and junctions to identify potential splice-associated variants. Compared with SpliceAI, MiSplice, Veridical, SAVNet. Applied to over 9,000 tumor samples (TCGA, MGI), used the experimentally validated 10 splice-site-creating variants that Jayasinghe et al. Input - BAM, GTF, VCF, configurable splice-region window definition. C++ implementation, [GitHub](https://github.com/griffithlab/regtools/), [Docker](https://hub.docker.com/r/griffithlab/regtools/).<details>
    <summary>Paper</summary>
    Cotto, Kelsy C, Yang-Yang Feng, Avinash Ramu, Megan Richters, Sharon Freshour, Zachary L Skidmore, Huiming Xia, et al. “RegTools: Integrated Analysis of Genomic and Transcriptomic Data for the Discovery of Splice-Associated Variants in Cancer.” Preprint. Bioinformatics, October 5, 2018. https://doi.org/10.1101/436634.
</details>

- [The Personal Cancer Genome Reporter (PCGR)](https://github.com/sigven/pcgr) is a stand-alone software package for functional annotation and translation of individual cancer genomes for precision cancer medicine. Currently, it interprets both somatic SNVs/InDels and copy number aberrations. The software extends basic gene and variant annotations from the [Ensembl’s Variant Effect Predictor (VEP)](http://www.ensembl.org/info/docs/tools/vep/index.html) with oncology-relevant, up-to-date annotations retrieved flexibly through vcfanno, and produces interactive HTML reports intended for clinical interpretation. NOTE: If you also want to interrogate the clinical impact of germline variants in the same individual, try the accompanying tool [Cancer Predisposition Sequencing Reporter (CPSR)](https://github.com/sigven/cpsr). R, Python. <details>
    <summary>Paper</summary>
    Sigve Nakken, Ghislain Fournous, Daniel Vodák, Lars Birger Aaasheim, Ola Myklebost, and Eivind Hovig. [Personal Cancer Genome Reporter: variant interpretation report for precision oncology](https://doi.org/10.1093/bioinformatics/btx817) (2017). Bioinformatics
</details>

- [PORI](https://bcgsc.github.io/pori/) (Platform for Oncogenomic Reporting and Interpretation) - analytical framework for reporting and interpreting somatic variants in cancer, prioritizes therapeutically informative alterations. Two components - GraphKB (a curated (resolving aliases, generalisms, etc., into controlled vocabulary) graph-based knowledge base integrating multiple [data sources](https://bcgsc.github.io/pori/graphkb/loading_data/#popular-resources), [API](https://github.com/bcgsc/pori_graphkb_api), [web](https://github.com/bcgsc/pori_graphkb_client), [Python](https://github.com/bcgsc/pori_graphkb_python)) and [IPR](https://github.com/bcgsc/pori_ipr_client) (Integrated Pipeline Reports, [API](https://github.com/bcgsc/pori_ipr_api), [web](https://github.com/bcgsc/pori_ipr_client), [Python](https://github.com/bcgsc/pori_ipr_python)). Flexible input variant call formats, patient data. [Demo](https://bcgsc.github.io/pori/demo/). <details>
    <summary>Paper</summary>
    Reisle, Caralyn, Laura M. Williamson, Erin Pleasance, Anna Davies, Brayden Pellegrini, Dustin W. Bleile, Karen L. Mungall, et al. “A Platform for Oncogenomic Reporting and Interpretation.” Nature Communications 13, no. 1 (December 2022): 756. https://doi.org/10.1038/s41467-022-28348-y.
</details>

- [Seattle](https://snp.gs.washington.edu/SeattleSeqAnnotation154/) - The SeattleSeq Annotation server provides annotation of SNVs (single-nucleotide variations) and small indels, both known and novel. This annotation includes dbSNP rs IDs, gene names and accession numbers, variation functions (e.g. missense), protein positions and amino-acid changes, conservation scores, HapMap frequencies, PolyPhen predictions, and clinical association

- [snpEff](https://pcingola.github.io/SnpEff/se_introduction/) - a variant annotation and effect prediction tool

- [The Variant Interpretation for Cancer Consortium Meta-Knowledgebase](https://search.cancervariants.org/). Aggregate interpretations covering 3,437 unique variants in 415 genes, 357 diseases, and 791 drugs. Validation using GENIE database. <details>
    <summary>Paper</summary>
    Wagner, Alex Handler, Brian Walsh, Georgia Mayfield, David Tamborero, Dmitriy Sonkin, Kilannin Krysiak, Jordi Deu Pons, et al. “A Harmonized Meta-Knowledgebase of Clinical Interpretations of Cancer Genomic Variants,” July 11, 2018. https://doi.org/10.1101/366856.
</details>

- [VARAdb](http://www.licpathway.net/VARAdb) - a comprehensive human variation annotation database, focus on noncoding regulatory . A total of 577,283,813 variations, five annotation sections including ‘Variation information’, ‘Regulatory information’, ‘Related genes’, ‘Chromatin accessibility’ and ‘Chromatin interaction’. The information section includes motif changes, risk SNPs, LD SNPs, eQTLs, clinical variant-drug-gene pairs, sequence conservation, somatic mutations, enhancers, super enhancers, promoters, TFs, ChromHMM states, histone modifications, ATAC accessible regions and chromatin interactions from Hi-C and ChIA-PET. Two types of variation related genes: 1) variation that sets in enhancer may associate with enhancer target genes predicted by Lasso method; 2) variation related genes based on distance. Can prioritize variations based on score, annotate novel variants and perform pathway downstream analysis. VCF annotation, one section at a time. Chromosome-specific [downloads](http://www.licpathway.net/VARAdb/download.php). <details>
    <summary>Paper</summary>
    Pan, Qi, Yue-Juan Liu, Xue-Feng Bai, Xiao-Le Han, Yong Jiang, Bo Ai, Shan-Shan Shi, et al. “[VARAdb: A Comprehensive Variation Annotation Database for Human](https://doi.org/10.1093/nar/gkaa922),” Nucleic Acids Research, 8 January 2021
</details>

- [VAtools](https://github.com/griffithlab/VAtools) - A set of tools to annotate VCF files with expression and readcount data. Python. [vatools.org](https://vatools.readthedocs.io) - codumentation.

- [VCF-plotein](https://vcfplotein.liigh.unam.mx/#/) - graphical, interactive interptetation of exome sequencing data in VCF format. Includes major databases, from gnomAD, COSMIC to Human Phenotype Ontology and GO terms. [Shiny app](https://vcfplotein.liigh.unam.mx/#/), [GitHub](https://github.com/raulossio/VCF-plotein). <details>
    <summary>Paper</summary>
    Ossio, Raul, O Isaac Garcia-Salinas, Diego Said Anaya-Mancilla, Jair S Garcia-Sotelo, Luis A Aguilar, David J Adams, and Carla Daniela Robles-Espinoza. “[VCF/Plotein: Visualisation and Prioritisation of Genomic Variants from Human Exome Sequencing Projects](https://doi.org/10.1093/bioinformatics/btz458).” Bioinformatics, June 4, 2019
</details>

- [vcfanno](https://github.com/brentp/vcfanno) - annotate a VCF with other VCFs/BEDs/tabixed files. Go, Python, binary. <details>
    <summary>Paper</summary>
    Pedersen, Brent S., Ryan M. Layer, and Aaron R. Quinlan. “Vcfanno: Fast, Flexible Annotation of Genetic Variants.” Genome Biology 17, no. 1 (December 2016). https://doi.org/10.1186/s13059-016-0973-5.
</details>

## SNP signatures

- [deconstructSigs](https://github.com/raerose01/deconstructSigs) - deconstructSigs aims to determine the contribution of known mutational processes to a tumor sample. By using deconstructSigs, one can: Determine the weights of each mutational signature contributing to an individual tumor sample; Plot the reconstructed mutational profile (using the calculated weights) and compare to the original input sample. Data from Alexandrov, COSMIC, others. R implementation. <details>
    <summary>Paper</summary>
    Rosenthal, Rachel, Nicholas McGranahan, Javier Herrero, Barry S. Taylor, and Charles Swanton. “DeconstructSigs: Delineating Mutational Processes in Single Tumors Distinguishes DNA Repair Deficiencies and Patterns of Carcinoma Evolution.” Genome Biology 17 (February 22, 2016): 31. https://doi.org/10.1186/s13059-016-0893-4.
</details>

- [MutationalPatterns](https://bioconductor.org/packages/MutationalPatterns/) - R package to evaluate and visualize a multitude of mutational patterns in base substitution catalogues of e.g. healthy samples, tumour samples, or DNA-repair deficient cells. The package covers a wide range of patterns including: mutational signatures, transcriptional and replicative strand bias, lesion segregation, genomic distribution and association with genomic features, which are collectively meaningful for studying the activity of mutational processes. The package works with single nucleotide variants (SNVs), insertions and deletions (Indels), double base substitutions (DBSs) and larger multi base substitutions (MBSs). R.

- [sigminer](https://shixiangwang.github.io/sigminer-doc/) - Genomic Alteration Signature Analysis, by Shixiang Wang, novel and known signature extraction, visualization. [GitHub](https://github.com/ShixiangWang/sigminer). R. [sigminer.prediction](https://github.com/ShixiangWang/sigminer.prediction) - Train and Predict Cancer Subtype with Keras Model based on Mutational Signatures.

- [SignatureAnalyzer](https://software.broadinstitute.org/cancer/cga/msp) - finding mutation patterns in multiple samples, NMF. Part of [CGA](https://software.broadinstitute.org/cancer/cga/) Cancer Genome Analysis tools.

- [SparseSignatures](https://github.com/danro9685/SparseSignatures), [R package](https://bioconductor.org/packages/release/bioc/html/SparseSignatures.html) - De Novo Mutational Signature Discovery in Tumor Genomes using SparseSignatures. <details>
    <summary>Paper</summary>
    Lal, Avantika, Keli Liu, Robert Tibshirani, Arend Sidow, and Daniele Ramazzotti. "De novo mutational signature discovery in tumor genomes using SparseSignatures." PLoS computational biology 17, no. 6 (2021): e1009119. https://doi.org/10.1371/journal.pcbi.1009119
</details>

- [YAPSA](https://bioconductor.org/packages/YAPSA/) - Yet Another Package for Signature Analysis, R package. Functionality used in [L. Alexandrov et al., Nature 2013](https://doi.org/10.1038/nature12477)


## SNP pathogenicity scores

- [ClinPred](https://sites.google.com/site/clinpred/home) - pathogenicity prediction for all nonsynonymous SNPs. Trained on ClinVar, validated on nine other databases. Random forest and gradient boosted decision tree, comparison with other machine learning algorithms. Downloadable scores for all nonsynonymous SNPs. <details>
    <summary>Paper</summary>
    Alirezaie, Najmeh, Kristin D. Kernohan, Taila Hartley, Jacek Majewski, and Toby Dylan Hocking. “ClinPred: Prediction Tool to Identify Disease-Relevant Nonsynonymous Single-Nucleotide Variants.” The American Journal of Human Genetics 103, no. 4 (October 2018): 474–83. https://doi.org/10.1016/j.ajhg.2018.08.005.
</details>

- [regBase](https://github.com/mulinlab/regBase) - Prediction of regulatory impact of variants outside of protein-coding regions, human. Trained on prediction scores from 23 tools, Gradient Tree Boosting, thorough training and evaluation. hg19 predictions are available for download. Python. <details>
    <summary>Paper</summary>
    Zhang, Shijie, Yukun He, Huanhuan Liu, Haoyu Zhai, Dandan Huang, Xianfu Yi, Xiaobao Dong, et al. “RegBase: Whole Genome Base-Wise Aggregation and Functional Prediction for Human Non-Coding Regulatory Variants.” Nucleic Acids Research, September 12, 2019, gkz774. https://doi.org/10.1093/nar/gkz774.
</details>

- [SEMpl](https://github.com/Boyle-Lab/SEM_CPP) - predict the impact of SNPs on TF binding. Uses ChIP-seq data, DNAse-seq, and PWMs. Simulate all possible SNPs in PWMs, estimate effect on ChIP signal. C, C++. <details>
    <summary>Paper</summary>
    Nishizaki, Sierra S, Natalie Ng, Shengcheng Dong, Robert S Porter, Cody Morterud, Colten Williams, Courtney Asman, Jessica A Switzenberg, and Alan P Boyle. “Predicting the Effects of SNPs on Transcription Factor Binding Affinity.” Edited by John Hancock. Bioinformatics, August 2, 2019, btz612. https://doi.org/10.1093/bioinformatics/btz612.
</details>

## SNP visualization, clustering

- [CONQUER](https://github.com/roderickslieker/CONQUER) - an R package for visualizing individual and multiple SNPs in epigenomic context, omics expression values, analysis of QTL genes for pathway enrichment. Interactive D3 graphics, circos plots. Similar functionality - [DEPICT](https://data.broadinstitute.org/mpg/depict/index.html) Input - rsIDs. <details>
    <summary>Paper</summary>
    Bouland, Gerard A, Joline W J Beulens, Joey Nap, and Arnaud Zaldumbide. “[CONQUER: An Interactive Toolbox to Understand Functional Consequences of GWAS Hits](https://doi.org/10.1093/nargab/lqaa085),” NAR Genomics and Bioinformatics, October 27, 2021, 7.
</details>

- [CMplot](https://github.com/YinLiLin/CMplot) - circular and regular (horizontal) Manhattan Plots, chromosome density plots, highly customizable.

- [gwasTools](https://github.com/bnwolford/gwasTools) - A collection of R scripts for exploring and plotting GWAS results. Visualization examples

- [gw](https://github.com/kcleal/gw) - A fast browser and annotation tool for genomic sequencing data

- [gpart](https://bioconductor.org/packages/gpart/) - R package for defining LD blocks (Big-LD algorithm), and visualizing them. <details>
    <summary>Paper</summary>
    - Ah Kim, Sun, Myriam Brossard, Delnaz Roshandel, Andrew D Paterson, Shelley B Bull, and Yun Joo Yoo. “Gpart: Human Genome Partitioning and Visualization of High-Density SNP Data by Identifying Haplotype Blocks.” Edited by Alfonso Valencia. Bioinformatics, May 9, 2019, btz308. https://doi.org/10.1093/bioinformatics/btz308.
</details>

- [maftools](https://bioconductor.org/packages/maftools/) - Summarize, Analyze and Visualize MAF files from TCGA or in house studies. [GitHub](https://github.com/PoisonAlien/maftools). R. 

- [manhattanly](https://cran.r-project.org/web/packages/manhattanly/) - Interactive Manhattan plots. R.

- [mutcraft](https://github.com/EmilieT/mutcraft) - R tools to mine & craft somatic mutations from cancer genomes. R.

- [MutScan](https://github.com/OpenGene/MutScan) - Detect and visualize target mutations by scanning FastQ files directly. C, C++.

- [samplot](https://github.com/ryanlayer/samplot) - Plot structural variant signals from many BAMs and CRAMs. Python.

- [ttplot](https://github.com/YTLogos/ttplot) - Tao Yan's Plot Toolkit, plots LD Heatmap, Manhattan plot. R.

- [VarClust](https://github.com/fasterius/VarClust) - A Python package for clustering of single nucleotide variants from high-through seqencing data. Works on single-sample VCF files.

- [VIVA](https://github.com/compbiocore/VariantVisualization.jl) - VCF visualization tool, written in Julia. Competing tools - vcfR, IGVZ, Genome Browser, Genome Savant, svviz, jvarkit - JfxNgs. Input - VCF file and, optionally, variant list, sample list, sample metadata. Filtering. Heatmap visualization. Julia. <details>
    <summary>Paper</summary>
    Tollefson, George A, Jessica Schuster, Fernando Gelin, Ashok Ragavendran, Isabel Restrepo, Paul Stey, James Padbury, and Alper Uzun. “VIVA (VIsualization of VAriants): A VCF File Visualization Tool.” BioRxiv, March 28, 2019. https://doi.org/10.1101/589879.
</details>

## SNP, GWAS databases

- [CAUSALdb](http://mulinlab.tmu.edu.cn/causaldb/index.html) - curated summary statistics for GWASs, mapped to MeSH terms, Manhattan plot visualization. Download available. <details>
    <summary>Paper</summary>
    Wang, Jianhua, Dandan Huang, Yao Zhou, Hongcheng Yao, Huanhuan Liu, Sinan Zhai, Chengwei Wu, et al. “CAUSALdb: A Database for Disease/Trait Causal Variants Identified Using Summary Statistics of Genome-Wide Association Studies.” Nucleic Acids Research, November 6, 2019, gkz1026. https://doi.org/10.1093/nar/gkz1026.
</details>

- [GWASatlas](https://atlas.ctglab.nl/) resource, analysis of pleiotropy, genetic architecture of complex traits. <details>
    <summary>Paper</summary>
    Watanabe, Kyoko, Sven Stringer, Oleksandr Frei, Masa Umićević Mirkov, Tinca J.C. Polderman, Sophie van der Sluis, Ole A. Andreassen, Benjamin M. Neale, and Danielle Posthuma. “A Global View of Pleiotropy and Genetic Architecture in Complex Traits.” BioRxiv, January 1, 2018, 500090. https://doi.org/10.1101/500090.
</details>

- [GWAScentral](https://www.gwascentral.org) - central GWAS repository. Browser, download. <details>
    <summary>Paper</summary>
    Beck, Tim, Tom Shorter, and Anthony J Brookes. “GWAS Central: A Comprehensive Resource for the Discovery and Comparison of Genotype and Phenotype Data from Genome-Wide Association Studies.” Nucleic Acids Research, October 15, 2019, gkz895. https://doi.org/10.1093/nar/gkz895.
</details>

- [gnomAD data](https://gnomad.broadinstitute.org/downloads) - [Open access to gnomAD data on multiple cloud providers](https://gnomad.broadinstitute.org/blog/2020-10-open-access-to-gnomad-data-on-multiple-cloud-providers/)

- High-coverage (30X) WGS and analyses of the original 2,504 1kGP samples, and 698 additional samples including 602 trios. More SNPs, InDels. [Data portal](https://www.internationalgenome.org/data-portal/data-collection/30x-grch38), Processed data, hg38: [GVCFs](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20190425_NYGC_GATK/raw_calls_updated/), [annotated SNP/INDEL VCFs](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/), [Phased SNV/INDEL/SV VCFs](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/), [SV VCF](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20210124.SV_Illumina_Integration/), [Sample metadata file with pedigree and sex information](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/1kGP.3202_samples.pedigree_info.txt). Detailed methods: QC (Picard Tools), SNV/InDel calling using GATK, BCFtools, VCFtools. Easy and difficult-to-sequence regions, chrX analysis. SV discovery using svtools (smoove, Manta, CNVnator tools in the wdl-based [workflow](https://github.com/hall-lab/sv-pipeline)), phasing. Comparison with phase 3 1000 genomes data. [List of tools](https://www.cell.com/action/showFullTableHTML?isHtml=true&tableId=undtbl1&pii=S0092-8674%2822%2900991-6). <details>
    <summary>Paper</summary>
    Byrska-Bishop, Marta, Uday S. Evani, Xuefang Zhao, Anna O. Basile, Haley J. Abel, Allison A. Regier, André Corvelo, et al. “High-Coverage Whole-Genome Sequencing of the Expanded 1000 Genomes Project Cohort Including 602 Trios.” Cell 185, no. 18 (September 2022): 3426-3440.e19. https://doi.org/10.1016/j.cell.2022.08.004.
</details>

- [kgp](https://github.com/stephenturner/kgp) - 1000 Genomes Project Metadata R Package, including high-coverage samples. By [Stephen Turner]

- Havrilla, James M., Brent S. Pedersen, Ryan M. Layer, and Aaron R. Quinlan. “[A Map of Constrained Coding Regions in the Human Genome](https://doi.org/10.1038/s41588-018-0294-6).” Nature Genetics, (2019) - Constrained coding regions (CCRs), analysis of gnomAD. [GitHub](https://github.com/quinlan-lab/ccr). [Data availability](https://www.nature.com/articles/s41588-018-0294-6#data-availability) section contain links to many genomics datasets. [Supplementary data](https://www.nature.com/articles/s41588-018-0294-6#Sec28):
    - [Supplementary Table 1](https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-018-0294-6/MediaObjects/41588_2018_294_MOESM3_ESM.txt) - Genes with CCRs in the 99th percentile or higher
    - [Supplementary Table 2](https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-018-0294-6/MediaObjects/41588_2018_294_MOESM4_ESM.txt) - CCRs under purifying selection specifically in humans
    - [Supplementary Table 3](https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-018-0294-6/MediaObjects/41588_2018_294_MOESM5_ESM.txt) - CCR enrichment in Pfam domains
    - [Supplementary Table 4](https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-018-0294-6/MediaObjects/41588_2018_294_MOESM6_ESM.txt) - Highly constrained CCRs not covered by missense depletion
    - [CCR Browser](https://s3.us-east-2.amazonaws.com/ccrs/ccr.html), [CCR BED files, autosomes](https://s3.us-east-2.amazonaws.com/ccrs/ccrs/ccrs.autosomes.v2.20180420.bed.gz), [CCR BED file, X chromosome](https://s3.us-east-2.amazonaws.com/ccrs/ccrs/ccrs.xchrom.v2.20180420.bed.gz)

- [PhenoScanner V2](http://www.phenoscanner.medschl.cam.ac.uk/) - a curated database from large GWASs. Disease/trait association with genetic variants, gene expression, metabolite and protein levels, and epigenetic markers. Search by any feature. Download, API access. [GitHub](https://github.com/phenoscanner/phenoscanner). R. <details>
    <summary>Paper</summary>
    Kamat, Mihir A, James A Blackshaw, Robin Young, Praveen Surendran, Stephen Burgess, John Danesh, Adam S Butterworth, and James R Staley. “PhenoScanner V2: An Expanded Tool for Searching Human Genotype-Phenotype Associations.” Edited by Janet Kelso. Bioinformatics, June 24, 2019, btz469. https://doi.org/10.1093/bioinformatics/btz469.
</details>

- [Public cancer GWAS](https://docs.google.com/document/d/13jmOtCNKiB8GPN0oTMJJOSVCQ_MUuO5bIG3zDKFHU0c/edit) by Peter Kraft. [Tweet](https://twitter.com/GENES_PK/status/1418543289968611330?s=20)

### GWAS pipelines

- [Eagle](https://www.hsph.harvard.edu/alkes-price/software/) - long-range phasing, improves the following imputation. Benchmarked against Beagle, HAPI-UR, SHAPEIT2, faster, similar or better phasing accuracy. <details>
    <summary>Paper</summary>
    Loh, Po-Ru, Pier Francesco Palamara, and Alkes L Price. “Fast and Accurate Long-Range Phasing in a UK Biobank Cohort.” Nature Genetics 48, no. 7 (July 2016): 811–16. https://doi.org/10.1038/ng.3571.
</details>

- [IMPUTE2](https://mathgen.stats.ox.ac.uk/impute/impute_v2.html) - prephasing framework, improves running time and accuracy of imputation step. In contrast to IMPUTE1 analytical integration strategy, IMPUTE2 uses a haplotype sampling strategy. <details>
    <summary>Paper</summary>
    Howie, Bryan, Christian Fuchsberger, Matthew Stephens, Jonathan Marchini, and Gonçalo R Abecasis. “Fast and Accurate Genotype Imputation in Genome-Wide Association Studies through Pre-Phasing.” Nature Genetics 44, no. 8 (August 2012): 955–59. https://doi.org/10.1038/ng.2354.
</details>

- [GCTA](https://yanglab.westlake.edu.cn/software/gcta/#Overview) - genome-wide complex trait analysis testing, estimates the variance explained by all SNPs on a chromosome or the whole genome for a complex trait, instead of testing the association of any particular SNP. Five function types: data management (PLINK, MACH formats), estimation of the genetic relationships from SNPs, mixed linear model analysis of variance explained, estimation of linkage disequilibrium, GWAS simulation. <details>
    <summary>Paper</summary>
    Yang, Jian, S. Hong Lee, Michael E. Goddard, and Peter M. Visscher. “GCTA: A Tool for Genome-Wide Complex Trait Analysis.” The American Journal of Human Genetics 88, no. 1 (January 2011): 76–82. https://doi.org/10.1016/j.ajhg.2010.11.011.
</details>

- [GLIMPSE](https://odelaneau.github.io/GLIMPSE/) - genotype likelihoods imputation and phasing method. Low WGS (1X) outperforms SNP arrays and is appropriate for GWASs and ancestry inference. Input - a matrix of genotype likelihoods from low-coverage sequencing data. Iterative genotype imputation and haplotype phasing with a Gibbs sampling procedure. Outperforms Beagle, LOIMPUTE, STITCH, GeneImp. [GitHub](https://github.com/odelaneau/GLIMPSE). [Tutorial 1](https://odelaneau.github.io/GLIMPSE/tutorial_b38.html), [tutorial 2](https://odelaneau.github.io/GLIMPSE/tutorial_hg19.html), [tutorial 3](https://odelaneau.github.io/GLIMPSE/tutorial_chrX.html). <details>
    <summary>Paper</summary>
    Rubinacci, Simone, Diogo M. Ribeiro, Robin J. Hofmeister, and Olivier Delaneau. “Efficient Phasing and Imputation of Low-Coverage Sequencing Data Using Large Reference Panels.” Nature Genetics 53, no. 1 (January 2021): 120–26. https://doi.org/10.1038/s41588-020-00756-0.
    Chat, Vylyny, Robert Ferguson, Leah Morales, and Tomas Kirchhoff. “Ultra Low-Coverage Whole-Genome Sequencing as an Alternative to Genotyping Arrays in Genome-Wide Association Studies.” Frontiers in Genetics 12 (February 15, 2022): 790445. https://doi.org/10.3389/fgene.2021.790445.
</details>

- [Michigan Imputation Server](https://imputationserver.sph.umich.edu/index.html#!) - phasing an imputation using [Eagle2](https://alkesgroup.broadinstitute.org/Eagle/) and [minimac3](https://genome.sph.umich.edu/wiki/Minimac3), outperforms Beagle, IMPUTE2. <details>
    <summary>Paper</summary>
    Das, Sayantan, Lukas Forer, Sebastian Schönherr, Carlo Sidore, Adam E Locke, Alan Kwong, Scott I Vrieze, et al. “Next-Generation Genotype Imputation Service and Methods.” Nature Genetics 48, no. 10 (October 2016): 1284–87. https://doi.org/10.1038/ng.3656.
</details>

- [SHAPEIT4](https://odelaneau.github.io/shapeit4/) - Segmented HAPlotype Estimation & Imputation Tool. Haplotype phasing, allows integrating external phasing information (large reference panels, pre-phased variants, long sequencing reads). SNP arrays of sequencing data, fast running time, scalable. Novel methods include (1) the Li and Stephens model to capture long-range haplotype sharing between individuals, (2) the Positional Burrows-Wheeler Transform (PBWT), (3) the compact representation of the solution space. Outperforms SHAPEIT3, Eagle2, Beagle5, tested on UK Biobank, GIAB data. [Website](http:://www.shapeit.fr) <details>
    <summary>Paper</summary>
    Delaneau, Olivier, Jean-François Zagury, Matthew R. Robinson, Jonathan L. Marchini, and Emmanouil T. Dermitzakis. “Accurate, Scalable and Integrative Haplotype Estimation.” Nature Communications 10, no. 1 (December 2019): 5436. https://doi.org/10.1038/s41467-019-13225-y.
</details>

- [Tractor](https://github.com/Atkinson-Lab/Tractor) - a statistical framework for GWAS analysis considering admixed individuals (local ancestry inference, LAI, estimated using [RFMix2](https://github.com/slowkoni/rfmix)). Tested on simulated and experimental two-way admixed African-European cohorts. Generates accurate ancestry-specific effect-size estimates and p-values, improves power. It is generally not necessary to include local ancestry in a GWAS model when there is no difference in the estimated effect size between ancestry groups. [Post-QC](https://github.com/Atkinson-Lab/Post-QC) - an automated pipeline for post-genotyping QC, harmonization, phasing, and LA inference. [ancestry_pipeline](https://github.com/armartin/ancestry_pipeline) - helper scripts for inferring local ancestry, performing ancestry-specific PCA, make painted karyogram plots. [Scripts](https://github.com/eatkinson/Tractor_ms_results) to reproduce the paper. The Tractor pipeline is implemented as scrips and Jupyter notebools on [GitHub](https://github.com/Atkinson-Lab/Tractor). [snowcat](https://github.com/andreyshabalin/snowcat) - Code Alternative for TRACTOR, by [Andrey Shabalin](https://github.com/andreyshabalin). <details>
    <summary>Paper</summary>
    Atkinson, Elizabeth G., Adam X. Maihofer, Masahiro Kanai, Alicia R. Martin, Konrad J. Karczewski, Marcos L. Santoro, Jacob C. Ulirsch, et al. “Tractor Uses Local Ancestry to Enable the Inclusion of Admixed Individuals in GWAS and to Boost Power.” Nature Genetics 53, no. 2 (February 2021): 195–204. https://doi.org/10.1038/s41588-020-00766-y.
</details>

### Ancestry

- [ADMIXTURE](https://dalexander.github.io/admixture/) - ancestry estimation using a model-based approach. Similar model as in STRUCTURE, FRAPPE, detailed description of model similarities and differences, fast numerical optimization algorithm. Contrasted with PCA-based EIGENSTRAT, STRUCTURE, nearly the same speed (as EIGENSTRAT) and performance on experimental and simulated data. Binaries for Linux and Mac. [Manual](https://dalexander.github.io/admixture/admixture-manual.pdf). <details>
    <summary>Paper</summary>
    Alexander, David H., John Novembre, and Kenneth Lange. “Fast Model-Based Estimation of Ancestry in Unrelated Individuals.” Genome Research 19, no. 9 (September 2009): 1655–64. https://doi.org/10.1101/gr.094052.109.
</details>

- [AKT](https://github.com/Illumina/akt) - ancestry and kinship toolkit by Illumina team. Detects related samples, ancestry, calculate correlation between variants, clustering, Mendel consistency. Input: multi-sample BCF file, no thinning required. Pre-built hg19 and hg38 1000 genomes PCA coordinates, projecton of samples on them, visualization. Algorithms and results are similar to PLINK, KING, faster. [Documentation](https://illumina.github.io/akt/). C++, command line. <details>
    <summary>Paper</summary>
    Arthur, Rudy, Ole Schulz-Trieglaff, Anthony J. Cox, and Jared O’Connell. “AKT: Ancestry and Kinship Toolkit.” Bioinformatics 33, no. 1 (January 1, 2017): 142–44. https://doi.org/10.1093/bioinformatics/btw576.
</details>

- [AIPS](https://github.com/biomedicaldatascience/AIPS) - Ancestry Inference using Principal Component analysis and Spatial analysis, R package. Distance-based, incorporates an Inverse Distance Weighted (IDW) interpolation method from spatial analysis. Compared with EIGENSTRAT, STRUCTURE, fastSTRUCTURE, and ADMIXTURE (each briefly described). AIPS and EIGENSTRAT principal components and eigenvalues are highly correlated (close to 1). R. <details>
    <summary>Paper</summary>
    Byun, Jinyoung, Younghun Han, Ivan P. Gorlov, Jonathan A. Busam, Michael F. Seldin, and Christopher I. Amos. “Ancestry Inference Using Principal Component Analysis and Spatial Analysis: A Distance-Based Analysis to Account for Population Substructure.” BMC Genomics 18, no. 1 (December 2017): 789. https://doi.org/10.1186/s12864-017-4166-8.
</details>

- [ancestry_pipeline](https://github.com/armartin/ancestry_pipeline) - Ancestry pipeline tutorial, by Alicia Martin. Provides helper scripts for inferring local ancestry, performing ancestry-specific PCA, etc. <details>
    <summary>More info</summary>
    Carrot-Zhang, Jian, Seunghun Han, Wanding Zhou, Jeffrey S. Damrauer, Anab Kemal, Andrew D. Cherniack, Rameen Beroukhim, et al. “Analytical Protocol to Identify Local Ancestry-Associated Molecular Features in Cancer.” STAR Protocols 2, no. 4 (December 2021): 100766. https://doi.org/10.1016/j.xpro.2021.100766 - Implementation of Alicia Martin's pipeline https://github.com/armartin/ancestry_pipeline. 
    <summary>TCGA ancestry</summary>
    Carrot-Zhang, Jian, Nyasha Chambwe, Jeffrey S. Damrauer, Theo A. Knijnenburg, A. Gordon Robertson, Christina Yau, Wanding Zhou, et al. “Comprehensive Analysis of Genetic Ancestry and Its Molecular Correlates in Cancer.” Cancer Cell 37, no. 5 (May 2020): 639-654.e6. https://doi.org/10.1016/j.ccell.2020.04.012. - Ancestry-centric analysis of somatic alterations, methylation, mRNA and miRNA expression, and DNA methylation across TCGA tumor types. Ancestry effects are tissue specific. Methods, tools, [GitHub](https://github.com/jcarrotzhang/ancestry-from-panel/tree/master/GDAN_AIM). Supplementary material https://www.cell.com/cancer-cell/fulltext/S1535-6108(20)30211-7#supplementaryMaterial: Table S1 - Admixture and Ethnicity Calls of all TCGA samples. Table S2. Mutation Associations with FBXW7, Alexandrov signatures. Table S3. Methylation Associations, genes associated with ancestry, models corrected for various confounders. Table S4. mRNA Associations. Table S5. miRNA Associations.
</details>

- [EIGENSTRAT](https://www.hsph.harvard.edu/alkes-price/software/) method, [FAQ](https://www.hsph.harvard.edu/alkes-price/eigensoft-frequently-asked-questions/). [GitHub](https://github.com/DReichLab/EIG) source code, Linux support, C.

- [EthSEQ](https://github.com/cibiobcg/EthSEQ) - Ethnicity Annotation from Whole-Exome and Targeted Sequencing Data. Input - reference genotype data of individuals with known enthnicity (e.g., 1000 genomes)  and either BAM files or VCF data. PCA using SNPRelate R package. Compared with fastSTRUCTURE. Multi-step refinement improves performance. 2000 SNPs are sufficient to each more than 98% precision. R package. <details>
    <summary>Paper</summary>
    Romanel, Alessandro, Tuo Zhang, Olivier Elemento, and Francesca Demichelis. “EthSEQ: Ethnicity Annotation from Whole Exome Sequencing Data.” Edited by Oliver Stegle. Bioinformatics 33, no. 15 (August 1, 2017): 2402–4. https://doi.org/10.1093/bioinformatics/btx165.
</details>

- [FastPCA](https://github.com/gabraham/flashpca) - population stratification by genotype, identifying SNPs differentiating subpopulations. Uses random matrix theory ([blanczos](http://tygert.com/blanczos.pdf) method) to approximate top PCs, reduced time and memory to linear scaling with the number of indivisuals. Compared with smartpca, PLINK2-pca, flashpca. Accuracy assessment on simulated data, nearly identical to PLINK2-pca. Combined witn SNPweights annotations. Applied to a [GERA](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000674.v1.p1) cohort of over 60K EAs genotyped on over 670K SNP array, filtered by QC (no sex, close relatedness, non european, non-autosomal, rare variants) and Hardy-Weinberg Equilibrium test and several rounds of LD pruning. Identified 5 subpopulations. C++, R. <details>
    <summary>Paper</summary>
    Galinsky, Kevin J., Gaurav Bhatia, Po-Ru Loh, Stoyan Georgiev, Sayan Mukherjee, Nick J. Patterson, and Alkes L. Price. “Fast Principal-Component Analysis Reveals Convergent Evolution of ADH1B in Europe and East Asia.” The American Journal of Human Genetics 98, no. 3 (March 2016): 456–72. https://doi.org/10.1016/j.ajhg.2015.12.022.
</details>

- [Gnomix](https://github.com/AI-sandbox/gnomix) - local ancestry inference (LAI) method, whole genome, higher accuracy, speed. Two stages: a set of classifiers to predict ancestry within genomic windows (Covariance reduction string kernel CovRSK performs best, others - tree-based (RF, XGBoost, etc.), linear (logistic regression, LDA, etc.), SVM, kNN) and Gnofix algorithm for smoothing, authomatic phasing. Outperforms RXMix, ELA, Loter. Applied to reconstruct LAI of dogs. Pretrained models, including human, are available. <details>
    <summary>Paper</summary>
    Hilmarsson, Helgi, Arvind S. Kumar, Richa Rastogi, Carlos D. Bustamante, Daniel Mas Montserrat, and Alexander G. Ioannidis. “High Resolution Ancestry Deconvolution for Next Generation Genomic Data.” Preprint. Genomics, September 21, 2021. https://doi.org/10.1101/2021.09.19.460980.
</details>

- [Neural ADMIXTURE](https://github.com/AI-sandbox/neural-admixture) - fast autoencoder network implementing ADMIXTURE principles to decompose genomes into fractional cluster assignments with each cluster representing a vector of DNA marker frequencies. Analog - NMF implemented as autoencoder. SNPs are encoded as 0, 1. Multi-head Neural ADMIXTURE combines multiple decoders to simultaneously obtain different clusters (aka ADMIXTURE runs with different priors for numbers of clusters). Encoders can be nonlinear, decoder - linear to help interpretability. Binary cross entropy loss function to minimize. Compared with ADMIXTURE, AIStructure, TeraStructure, HaploNet. PyTorch implementation. <details>
    <summary>Paper</summary>
    Mantes, Albert Dominguez, Daniel Mas Montserrat, Carlos D. Bustamante, Xavier Giró-i-Nieto, and Alexander G. Ioannidis. “Neural ADMIXTURE: Rapid Population Clustering with Autoencoders.” Preprint. Genomics, June 28, 2021. https://doi.org/10.1101/2021.06.27.450081.
</details>

- [RFMix](https://github.com/slowkoni/rfmix) - local ancestry inference (LAI), a discriminative approach that models ancestry along an admixed chromosome given observed haplotype sequences of reference ancestry (known or inferred). Phased chromosomes from reference populations (HapMap3) are divided into windows and, within each window, a random forest is trained to distinguish ancestry by using reference data and generates a fractional vote for each ancestry that are summed producing posterior ancestry probabilities within each window. These posterior probabilities are used for determining the most likely sequence of ancestry across windows via MAP inference. Phasing using Beagle. Outperforms LAMP-HAP, SupportMix. [Manual](https://github.com/slowkoni/rfmix/blob/master/MANUAL.md). <details>
    <summary>Paper</summary>
    Maples, Brian K., Simon Gravel, Eimear E. Kenny, and Carlos D. Bustamante. “RFMix: A Discriminative Modeling Approach for Rapid and Robust Local-Ancestry Inference.” The American Journal of Human Genetics 93, no. 2 (August 2013): 278–88. https://doi.org/10.1016/j.ajhg.2013.06.020.
</details>

- [somalier](https://github.com/brentp/somalier) - fast sample-swap and relatedness checks on BAMs/CRAMs/VCFs/GVCFs, by Brent Pedersen. Commands: extract, relate, ancestry, find-sites. [Ancestry prediction wiki](https://github.com/brentp/somalier/wiki/ancestry)

- [SNPweights](https://www.hsph.harvard.edu/alkes-price/software/) - ancestry inference using SNP weights from external reference panels. Greater accuracy than ancestry-informative markers (AIMs are 10-fold more informative than random SNPs for ancestry inference). PCA on the reference panel, prediction of PCs for the unseen samples. SNPs excluded: chrX, A/T, C/G. Predicted PCs using genome-wide SNPs have the highest accuracy. SNP weights for various populations are available. <details>
    <summary>Paper</summary>
    Chen, Chia-Yen, Samuela Pollack, David J. Hunter, Joel N. Hirschhorn, Peter Kraft, and Alkes L. Price. “Improved Ancestry Inference Using Weights from External Reference Panels.” Bioinformatics 29, no. 11 (June 1, 2013): 1399–1406. https://doi.org/10.1093/bioinformatics/btt144.
</details>

- [STRUCTURE](https://web.stanford.edu/group/pritchardlab/structure.html) - a model-based clustering method for inferring population structure from genotypes and assign individuals to populations. Assumptions - K populations (may be unknown), each is characterized by a set of allele frequencies at each locus, loci are unlinked and in linkage disequilibrium and Hardy-Weinberg equilibrium in populations. Model-based clustering, observations from each cluster are random draws from some parametric model. Inference using a Bayesian approach. MCMC and Gibbs sampling, overview in the [appendix](https://academic.oup.com/genetics/article/155/2/945/6048111#325587632). <details>
    <summary>Paper</summary>
    Pritchard, Jonathan K, Matthew Stephens, and Peter Donnelly. “Inference of Population Structure Using Multilocus Genotype Data.” Genetics 155, no. 2 (June 1, 2000): 945–59. https://doi.org/10.1093/genetics/155.2.945.
</details>

- Ancestry inference from any genomics data. Data synthesis framework - repeatedly adding population reference genomics data with known ancestry to the data of interest with unknown ancestry and assessing performance in inferring ancestry (PCA, five dimensions, K nearest neighborhood clustering). Idea - if ancestry backgrounds match, addition of the reference genomics data will improve performance. Code upon request. <details>
    <summary>Paper</summary>
    Belleau, Pascal, Astrid Deschênes, David A. Tuveson, and Alexander Krasnitz. “Accurate and Robust Inference of Genetic Ancestry from Cancer-Derived Molecular Data across Genomic Platforms.” Preprint. Bioinformatics, February 4, 2022. https://doi.org/10.1101/2022.02.01.478737.
</details>

- [ancestryinformativemarkers](https://github.com/leipzig/ancestryinformativemarkers) - Public AIMs, by [Jeremy Leipzig](https://github.com/leipzig)

- AIM panels have been developed for European Americans (Paschou et al., 2008; Price et al., 2008; Seldin and Price, 2008; Tian et al., 2008) and Latino Americans (Galanter et al., 2012). [Source](https://doi.org/10.1093/bioinformatics/btt144)

### eQTLs

- [The eQTL catalog](https://www.ebi.ac.uk/eqtl/) providing uniformly processed gene expression, imputed genotypes, and splicing QTLs from [29 human public studies](https://www.ebi.ac.uk/eqtl/Studies/) (GTeX, BLUEPRINT, GEUVADIS etc.). 69 distinct cell types and tissues, harmonized using [Zooma](https://www.ebi.ac.uk/spot/zooma/) ontology annotation. High reproducibility across datasets. [Methods, pipelines, and code](https://www.ebi.ac.uk/eqtl/Methods/) for RNA-seq quantification, QC and normalization, genotype QC and imputation, association testing. [FTP](ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL) with tab-delimited and HDF5 data summarized by chromosome, SNP, study etc. [Links to all tab-delimited files](https://github.com/kauralasoo/eQTL-Catalogue-resources/blob/master/tabix/). [REST API](https://www.ebi.ac.uk/eqtl/api-docs/), [Tabix](https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources/blob/master/tutorials/tabix_use_case.md), and [Python](https://github.com/eQTL-Catalogue/eQTL-SumStats/blob/master/querying_hdf5_basics.ipynb) access options are available. <details>
    <summary>Paper</summary>
    Kerimov, Nurlan. et al. “[A Compendium of Uniformly Processed Human Gene Expression and Splicing Quantitative Trait Loci](https://doi.org/10.1038/s41588-021-00924-w).” Nature Genetics, 2021.
</details>


### Polygenic risk score

- [Polygenic risk score (PRS) analysis guidelines and tutorial](https://choishingwan.github.io/PRS-Tutorial/). PRS is an estimate of an individual’s genetic liability to a trait or disease, calculated according to their genotype profile and relevant genome-wide association study (GWAS) data. Box 1 - definitions.  Several options how PRSs are calculated. Data and R/command line examples using PLINK, PRSice-2, LDpred-2 and lassosum for QC of base and target data, PRS calculations, visualizing results. <details>
    <summary>Paper</summary>
    Choi, S.W., Mak, T.S. & O’Reilly, P.F. Tutorial: a guide to performing polygenic risk score analyses. Nat Protoc (2020). https://doi.org/10.1038/s41596-020-0353-1
</details>

- [pgsc_calc](https://github.com/PGScatalog/pgsc_calc) - The Polygenic Score Catalog Calculator. Nextflow.

### Regulatory

- [Activity-by-Contact (ABC) Model](https://www.engreitzlab.org/resources/) predicting SNP-enhancer regularoty pairs. [Predictions](ftp://ftp.broadinstitute.org/outgoing/lincRNA/ABC/AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz) in 131 cell types and tissues (all element-gene connections with ABC scores >= 0.015. [GitHub](https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction).

- [3DGeOD](https://www.csuligroup.com/3DGeOD/home) - The 3D Genome Organization in Disease curated database of disease-associated genetic variants (inter-chromosomal translocations ICTs from Mitelman, COSMIC, ChimerPub databases, SVs from dbSNP, ClinVar, SNPs from eQTL catalogue) disrupting 3D structure (Hi-C data from 4D Nucleoume, territory, compartment, TAD, loop boundaries) and the associated gene expression in 35 tissue/cell types. Search for curated associations, predicted with [3DFunc](https://github.com/CSUBioGroup/3DFunc), a two-phase scoring algorithm to prioritize gene expression changes associated with 3D changes. [GitHub](https://github.com/CSUBioGroup/3DFunc). <details>
    <summary>Paper</summary>
    Tang, Li, Matthew C. Hill, Mingxing He, Junhao Chen, Patrick T. Ellinor, and Min Li. “A 3D Genome Atlas of Genetic Variants and Their Pathological Effects.” Preprint. Bioinformatics, November 27, 2022. https://doi.org/10.1101/2022.11.27.518071.
</details>

- [SNP2TFBS](https://ccg.epfl.ch/snp2tfbs/) is a Web interface aimed at studying variations (SNPs/indels) that affect transcription factor binding (TFB) in the Human genome.

## InDels

- Kosugi, Shunichi, Yukihide Momozawa, Xiaoxi Liu, Chikashi Terao, Michiaki Kubo, and Yoichiro Kamatani. “[Comprehensive Evaluation of Structural Variation Detection Algorithms for Whole Genome Sequencing](https://doi.org/10.1186/s13059-019-1720-5).” Genome Biology, (December 2019) - Benchmarking of structural variant detection tools. Introduction to types of structural variants. No tool detects all. Table 1 prioritizes best tools for deletion, duplication, insertion, invertion detection.

- [Breakfast](https://github.com/annalam/breakfast) - a software for detecting genomic structural variants from DNA sequencing data.

- [MindTheGap](https://gatb.inria.fr/software/mind-the-gap/) - detection and assembly of DNA insertion variants. [GitHub](https://github.com/GATB/MindTheGap). C++.

- [Pindel](http://gmt.genome.wustl.edu/packages/pindel/) - breakpoints of large deletions, medium sized insertions, inversions, tandem duplications and other structural variants. C++.

- [SVAFotate](https://github.com/fakedrtom/SVAFotate) - Annotate a (lumpy) structual variant (SV) VCF with allele frequencies (AFs) from large population SV cohorts. Python.


## CNV, SV

- [Awesome papers and projects about CNV and SV using NGS data](https://github.com/geocarvalho/sv-cnv-studies) - Relevant studies with Structual Variants and Copy Number Variants in NGS (Genome, Exome and Amplicon Sequencing) pipelines

- Benchmarking of 15 WGS-based structural variant callers (focusing on deletions), selected out of 55 (Supplementary Table 1). Gold standard - homozygous deletions in inbred mouse chromosomes. Lumpy, Manta give most optimal performance under various test strategies. [GitHub](https://github.com/Mangul-Lab-USC/benchmarking_SV), [Tweet](https://twitter.com/brent_p/status/1265803854785855490?s=20). <details>
    <summary>Paper</summary>
    Sarwal, Varuni, Sebastian Niehus, Ram Ayyala, Sei Chang, Angela Lu, Nicholas Darci-Maher, Russell Littman, et al. “[A Comprehensive Benchmarking of WGS-Based Structural Variant Callers](https://doi.org/10.1101/2020.04.16.045120).” Preprint. April 18, 2020. 
</details>

- Benchmark of 10 CNV callers. LUMPY performs best overall, Canvas is good for high specificity, CNVnator and RDXplorer are good for high sensitivity. Table 1 summarizes functionality of each caller. Used the Database of Genomic Variants as a gold standard, call CNVs from NA12878 genome. <details>
    <summary>Paper</summary>
    - Zhang, Le, Wanyu Bai, Na Yuan, and Zhenglin Du. “[**Comprehensively Benchmarking Applications for Detecting Copy Number Variation.**](https://doi.org/10.1371/journal.pcbi.1007069)” PLOS Computational Biology 15, no. 5 (May 28, 2019)
</details>

- Review of structural variant callers. De novo-based approaches (graph- or scaffold-based), short-read DNA-seq and RNA-seq (gene fusion) mapping, long-read (PacBio, Oxford Nanopore) mapping, multimethods approaches. SV calling from newer technologies, such as optical mapping, strand-seq, 10X Genomics linked reads, Hi-C. Brief description of tools, their performance, references to reviews. [Table 1](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1828-7/tables/1) - categorized list of tools, brief description and links. <details>
    <summary>Paper</summary>
    Mahmoud, Medhat, Nastassia Gobet, Diana Ivette Cruz-Dávalos, Ninon Mounier, Christophe Dessimoz, and Fritz J. Sedlazeck. “[Structural Variant Calling: The Long and the Short of It](https://doi.org/10.1186/s13059-019-1828-7).” Genome Biology 20, no. 1 (December 2019): 246. 
</details>

- Benchmarking of seven SV callers (BreakDancer, Pindel, Delly, Lumpy, Manta, GRIDSS, SvABA) to detect different SV types, sizes, the effect of SV abundance and sequencing coverage, sequence similarity, biases (GC content and homopolymers), and mapping quality. Overview of read-pair, read-depth, split-read, and local assembly. SV types and definitions. Manta, Lumpy, of GRIDSS perform well. [Supplementary material](https://academic.oup.com/bib/article/22/3/bbaa056/5831479#supplementary-data) - code examples. <details>
    <summary>Paper</summary>
    Gong, Tingting, Vanessa M Hayes, and Eva K F Chan. “[Detection of Somatic Structural Variants from Short-Read next-Generation Sequencing Data](https://doi.org/10.1093/bib/bbaa056),” Briefings in Bioinformatics, May 2021
</details>

- [ClassifyCNV](https://github.com/Genotek/ClassifyCNV) - clinical annotation of copy-number variants. Input - BED file with DEL/DUP. Output - numeric pathogenicity scores. Python3

- [CopywriteR](https://github.com/PeeperLab/CopywriteR) - R package for CNV detection from WES off-target sequences, discarding reads covering capture bait regions and adjacent. Calling peaks with MACS 1.4 and removing reads overlapping peaks. Compared with onTarget CBS (circular binary segmentation), propSeq, EXCAVATOR, higher sensitivity, better than low-coverage WGS, more robust to noise, improves downstream segmentation analysis. Works without reference. Performs well on ChIP-seq data (similarly removing reads in MACS peaks). 20kb bin resolution, bins are filtered, corrected for GC content, mappability, excluding blacklist-overlapping reads. Parallelized.  preCopywriteR(), CopywriteR() and plotCNA() functions. BAM files as an input. <details>
    <summary>Paper</summary>
    Kuilman, Thomas, Arno Velds, Kristel Kemper, Marco Ranzani, Lorenzo Bombardelli, Marlous Hoogstraat, Ekaterina Nevedomskaya, et al. “CopywriteR: DNA Copy Number Detection from off-Target Sequence Data.” Genome Biology 16, no. 1 (December 2015): 49. https://doi.org/10.1186/s13059-015-0617-1.
</details>

- [Control-FREEC](http://boevalab.com/FREEC/) - assess copy number and genotype information in whole genome and exome sequencing data. Corrects for contamination by normal cells and variable sample ploidy. With a matched normal sample, distinguishes somatic from germline events. R. <details>
    <summary>Paper</summary>
    Boeva, Valentina, Tatiana Popova, Kevin Bleakley, Pierre Chiche, Julie Cappo, Gudrun Schleiermacher, Isabelle Janoueix-Lerosey, Olivier Delattre, and Emmanuel Barillot. “[Control-FREEC: A Tool for Assessing Copy Number and Allelic Content Using next-Generation Sequencing Data](https://doi.org/10.1093/bioinformatics/btr670).” Bioinformatics, (February 1, 2012)
</details>

- [copyCat](https://github.com/chrisamiller/copyCat) - a parallel R package for detecting copy-number alterations from short sequencing reads. To call high-coverage regions. Has downloadable annotations (mappability, GC-content). R.

- [copynumber](https://bioconductor.org/packages/copynumber/) - Segmentation of single- and multi-track copy number data by penalized least squares regression. Similar functionality to [DNAcopy](https://bioconductor.org/packages/DNAcopy/). R.

- [CREST](https://www.stjude.org/research/labs/zhang-lab/crest.html) - mapping somatic structural variations in cancer genomes ad base-pair resolution. realignment of soft-clipped subsequences.Compared with BreakDancer, Pindel. Tested on experimental and simulated data, hg18, some SVs experimentally validated. <details>
    <summary>Paper</summary>
    Wang, Jianmin, Charles G Mullighan, John Easton, Stefan Roberts, Sue L Heatley, Jing Ma, Michael C Rusch, et al. “[CREST Maps Somatic Structural Variation in Cancer Genomes with Base-Pair Resolution](https://doi.org/10.1038/nmeth.1628).” Nature Methods 8, no. 8 (August 2011)
</details>

- [CNVkit](https://github.com/etal/cnvkit) - capturing CNVs in on-target and off-target genomic regions. Existing tools (CNVer, ExomeCNV, exomeCopy, CONTRA, CoNIFER, ExomeDepth, VarScan2, XHMM, ngCGH, EXCAVATOR, CANOES, PatternCNV, CODEX, Control-FREEC, cn.MOPS, cnvOffSeq, CopyWriteR). Account for GC content, mappability. [Examples](https://github.com/etal/cnvkit-examples). Python. <details>
    <summary>Paper</summary>
    Talevich, Eric, A. Hunter Shain, Thomas Botton, and Boris C. Bastian. “[CNVkit: Genome-Wide Copy Number Detection and Visualization from Targeted DNA Sequencing](https://doi.org/10.1371/journal.pcbi.1004873).” PLOS Computational Biology, (April 21, 2016)
</details>

- [CNVpytor](https://github.com/abyzovlab/CNVpytor) - CNV/CNA detection from read depth and SNPs (allelic imbalance). Python implementation of CNVnator, faster parsing (pysam) smaller .pytor files (H5). Works on a cloud, in Jupyter notebooks. Visualization using JBrowse. Python. <details>
    <summary>Paper</summary>
    Suvakov, Milovan, Arijit Panda, Colin Diesh, Ian Holmes, and Alexej Abyzov. “[CNVpytor: A Tool for CNV/CNA Detection and Analysis from Read Depth and Allele Imbalance in Whole Genome Sequencing](https://doi.org/10.1101/2021.01.27.428472
),” biorXiv, January 27, 2021
</details>

- [CNVnator](https://github.com/abyzovlab/CNVnator) - CNV using read depth. Bin the genome, mean-shift technique to quantify CNVs. Poor agreement in methods detection. Misses retrotransposon-CNVs but detects >50% other regions missed by other methods. C++, Python. <details>
    <summary>Paper</summary>
    Abyzov, Alexej, Alexander E. Urban, Michael Snyder, and Mark Gerstein. “[**CNVnator: An Approach to Discover, Genotype, and Characterize Typical and Atypical CNVs from Family and Population Genome Sequencing.**](https://doi.org/10.1101/gr.114876.110)” Genome Research 21, no. 6 (June 2011)
</details>

- [DELLY](https://tobiasrausch.com/delly/) - detection of structural variants, such as CNVs, duplications, inversions, translocations. Paired-end and split-read analysis. C, C++. <details>
    <summary>Paper</summary>
    Rausch, Tobias, Thomas Zichner, Andreas Schlattl, Adrian M. Stütz, Vladimir Benes, and Jan O. Korbel. “[DELLY: Structural Variant Discovery by Integrated Paired-End and Split-Read Analysis](https://doi.org/10.1093/bioinformatics/bts378).” Bioinformatics, (September 15, 2012)
</details>

- [DNAcopy](https://bioconductor.org/packages/DNAcopy/) - DNA copy number data analysis using Circular Binary Segmentation (CBS) algorithm. Input - CHG-like data, regions with comparative copy number values. Smoothing, outlier detection, segmentation using changepoint detection, visualization. R.

- [dysgu](https://github.com/kcleal/dysgu) - a collection of tools for calling structural variants using short paired-end or long reads (ONT, PacBio). Detects SVs from alignment gaps, discordant and supplementary mappings (a fast consensus sequence algorithm, inspired by the positional de Brujin graph, followed by remapping of anomalous sequences to discover additional small SVs), generates consensus contigs, classifies events using macnine learning (pre-trained on manually curated gold standard), generates confidence score. Outperforms [svbench](https://github.com/kcleal/svbench) library for benchmarking) manta, delly, lumpy, strelka, gatk on GIAB, simulated datasets, fast. Cython, Python. <details>
    <summary>Paper</summary>
    Cleal, Kez, and Duncan Baird. "Dysgu: efficient structural variant calling using short or long reads." Nucleic Acids Research (31 January 2022). https://doi.org/10.1093/nar/gkac039
</details>

- [instasv](https://github.com/venyao/intansv) - Integrative analysis of structural variations called by different tools. R.

- [gnomAD-SV](https://gnomad.broadinstitute.org/) - structural variants from deep WGS, added to gnomAD browser. 14,891 genomes, average statistics of SVs in general population. <details>
    <summary>Paper</summary>
    Collins, Ryan L., Harrison Brand, Konrad J. Karczewski, Xuefang Zhao, Jessica Alföldi, Amit V. Khera, Laurent C. Francioli, et al. “[An Open Resource of Structural Variation for Medical and Population Genetics](https://doi.org/10.1101/578674).” BioRxiv, March 14, 2019.
</details>

- [GRIDSS2](https://github.com/PapenfussLab/gridss) - structural variant caller detecting single breakpoints. Phasing of complex rearrangements. Assembles reads supporting a structural variant using a positional de Bruijn graph. Single breakends are called on soft-clipped reads, reads and assemblies with unmapped or ambiguously mapping mates (Methods). Benchmarked against GRIDSS1, Manta, svaba, novobreak on manually curated experimental and simulated data. Input - one or more aligned SAM/BAM. 5 phases: (i) preprocessing, (ii) assembly, (iii) variant calling, (iv) annotation, and (v) somatic filtering. Annotations go into INFO and FORMAT VCF fields. [Quickstart](https://github.com/PapenfussLab/gridss/blob/master/QuickStart.md). Java. <details>
    <summary>Paper</summary>
    Cameron, Daniel L., Jonathan Baber, Charles Shale, Jose Espejo Valle-Inclan, Nicolle Besselink, Arne van Hoeck, Roel Janssen, Edwin Cuppen, Peter Priestley, and Anthony T. Papenfuss. “GRIDSS2: Comprehensive Characterisation of Somatic Structural Variation Using Single Breakend Variants and Structural Variant Phasing.” Genome Biology 22, no. 1 (December 2021): 202. https://doi.org/10.1186/s13059-021-02423-x.
</details>

- [HATCHet](https://github.com/raphael-group/hatchet) (Holistic Allele-specific Tumor Copy-number Heterogeneity) is an algorithm that infers allele and clone-specific CNAs and WGDs jointly across multiple tumor samples from the same patient, and that leverages the relationships between clones in these samples. Roff, Conda. <details>
    <summary>Paper</summary>
    - Zaccaria, Simone, and Benjamin J. Raphael. “[Accurate Quantification of Copy-Number Aberrations and Whole-Genome Duplications in Multi-Sample Tumor Sequencing Data](https://doi.org/10.1101/496174).” BioRxiv, January 1, 2018
</details>

- [LUMPY-SV](https://github.com/arq5x/lumpy-sv/) - a general probabilistic framework for structural variant discovery. Integrates multiple signals - read-pair, split-read, read-depth and prior knowledge. Operates on multiple samples. C, C++. <details>
    <summary>Paper</summary>
    Layer, Ryan M, Colby Chiang, Aaron R Quinlan, and Ira M Hall. “[**LUMPY: A Probabilistic Framework for Structural Variant Discovery.**](https://doi.org/10.1186/gb-2014-15-6-r84)” Genome Biology 15, no. 6 (2014): R84. 
</details>

- [Manta](https://github.com/Illumina/manta) - SV detection in single- and tumor-normal samples. parallelized for within-sample performance. Fast, detects more variants of different types. C++, Python. <details>
    <summary>Paper</summary>
    Chen, Xiaoyu, Ole Schulz-Trieglaff, Richard Shaw, Bret Barnes, Felix Schlesinger, Morten Källberg, Anthony J. Cox, Semyon Kruglyak, and Christopher T. Saunders. “[Manta: Rapid Detection of Structural Variants and Indels for Germline and Cancer Sequencing Applications](https://doi.org/10.1093/bioinformatics/btv710).” Bioinformatics, (April 15, 2016)
</details>

- [QDNASeq](https://bioconductor.org/packages/release/bioc/html/QDNAseq.html) - Quantitative DNA sequencing for chromosomal aberrations. R.

- [samplot](https://github.com/ryanlayer/samplot) - Plot structural variant signals from many BAMs and CRAMs. 

- [smoove](https://github.com/brentp/smoove) - structural variant calling and genotyping with existing tools, but, smoothly. Go. Docker.

- [SV-Bay](https://github.com/BoevaLab/SV-Bay) - structural variant detection, Bayesian, corrects for GC-content and mappability. Uses both normal and abnormal reads, paired-end and depth information. Somatic variants if a normal sample is available. Detailed methods, stats. Compared with GASVPro, Lumpy, BreakDancer, DELLY on simulated (TGsim) and experimental neuroblastoma datasets. Improve sensitivity and specificity of SV detection, less false positives. Reasonably fast. Python, R. <details>
    <summary>Paper</summary>
    Iakovishina, Daria, Isabelle Janoueix-Lerosey, Emmanuel Barillot, Mireille Regnier, and Valentina Boeva. “[SV-Bay: Structural Variant Detection in Cancer Genomes Using a Bayesian Approach with Correction for GC-Content and Read Mappability](https://doi.org/10.1093/bioinformatics/btv751).” Bioinformatics, (April 1, 2016)
</details>


## Power

- [GeneticsDesign](https://www.bioconductor.org/packages/release/bioc/html/GeneticsDesign.html) - GWAS power analysis, functions for designing genetics studies. R.

- [Genetic Association Study (GAS) Power Calculator](http://csg.sph.umich.edu/abecasis/gas_power_calculator/index.html) - online tool to compute statistical power for large one-stage genetic association studies. The underlying method is derived from the CaTS power calculator for two-stage association studies (2006). 

- [powerEQTL](https://bwhbioinfo.shinyapps.io/powerEQTL/) - power analysis for eQTL studies. Includes two models, one-way unbalanced ANOVA (categorical genotypes) and linear regression (additive counting genotypes). Applicable to bulk and scRNA-seq. [Shiny app](https://bwhbioinfo.shinyapps.io/powerEQTL/), [CRAN](https://CRAN.R-project.org/package=powerEQTL). <details>
    <summary>Paper</summary>
    Dong, Xianjun, Xiaoqi Li, Tzuu-Wang Chang, Scott T Weiss, and Weiliang Qiu. “[PowerEQTL: An R Package and Shiny Application for Sample Size and Power Calculation of Bulk Tissue and Single-Cell EQTL Analysis](https://doi.org/10.1101/2020.12.15.422954),” bioRxiv, December 16, 2020.
</details>

- [SEQPower](http://bioinformatics.org/spower/start) - GWAS power analysis for case/control and quantitative studies, rare variants. Shell, binary. <details>
    <summary>Paper</summary>
    Wang, Gao T., Biao Li, Regie P. Lyn Santos-Cortez, Bo Peng, and Suzanne M. Leal. “Power Analysis and Sample Size Estimation for Sequence-Based Association Studies.” Bioinformatics 30, no. 16 (August 15, 2014): 2377–78. https://doi.org/10.1093/bioinformatics/btu296.
</details>

## Miscellaneous

- [Consensus ICGC SNV](https://www.synapse.org/#!Synapse:syn7357330) - The set of somatically acquired SNVs and indels across PCAWG tumour samples. Variant calls were generated by three pipelines run independently on each sample, with subsequent merging into a consensus set of high-quality calls.

- [Consensus ICGC CNA](https://www.synapse.org/#!Synapse:syn8042988) - Somatically acquired copy number alterations across PCAWG tumour samples. Variant calls were generated by three pipelines run independently on each sample, with subsequent merging into a consensus set of high-quality calls. 

- [Awesome papers and projects about CNV and SV using NGS data](https://github.com/geocarvalho/sv-cnv-studies)

- [Atlas of Variant Age](https://human.genome.dating/) - age estimate of \~45M SNPs. Method - Genealogical Estimation of Variant Age (GEVA), performs similar or better to PSMC. https://human.genome.dating/
    - Albers, Patrick K., and Gil McVean. “Dating Genomic Variants and Shared Ancestry in Population-Scale Sequencing Data.” Edited by Nick H. Barton. PLOS Biology 18, no. 1 (January 17, 2020): e3000586. https://doi.org/10.1371/journal.pbio.3000586.

- [chromeister](https://github.com/estebanpw/chromeister) - An ultra fast, heuristic approach to detect conserved signals in extremely large pairwise genome comparisons (dot-plot). C, R, shell.

- [DNA-seq-analysis](https://github.com/crazyhottommy/) DNA sequencing analysis notes from Ming Tang. 

- [gggenomes](https://github.com/thackl/gggenomes) - A grammar of graphics for comparative genomics. R. [Website](https://thackl.github.io/gggenomes/)

- [genomepy](https://github.com/vanheeringen-lab/genomepy) - Download genomes the easy way. Python.

- [GWAS tutorial](https://github.com/MareesAT/GWA_tutorial/) Quality control with PLINK, population stratification (MDS), association tests (binary and quantitative using PLINK), polygenic risk scores. PLINK file formats. Box 1 - GWAS definitions. <details>
    <summary>Paper</summary>
    Marees, Andries T., Hilde de Kluiver, Sven Stringer, Florence Vorspan, Emmanuel Curis, Cynthia Marie-Claire, and Eske M. Derks. “A Tutorial on Conducting Genome-Wide Association Studies: Quality Control and Statistical Analysis.” International Journal of Methods in Psychiatric Research 27, no. 2 (2018): e1608. https://doi.org/10.1002/mpr.1608.
</details>

- Long-range sequencing and software for data processing. SMRT by Pacific Bioscience, nanopore-based by Oxford Nanopore, genome partitioning and barcoding by 10X Genomics, Hi-C based, BioNano optical mapping (Table 1). Table 2 - Bioinformatics tools for de novo genome assembly, SNP/CNV etc. variant detection, phasing, RNA-seq, methylation. Text describes each tool. <details>
    <summary>Paper</summary>
    Sedlazeck, Fritz J., Hayan Lee, Charlotte A. Darby, and Michael C. Schatz. “[Piercing the Dark Matter: Bioinformatics of Long-Range Sequencing and Mapping](https://doi.org/10.1038/s41576-018-0003-4).” Nature Reviews Genetics 19, no. 6 (June 2018)
</details>

- [popgen-notes](https://github.com/cooplab/popgen-notes) - Population genetics notes, by Graham Coop. [PDF](https://github.com/cooplab/popgen-notes/releases), [R demo code](https://github.com/cooplab/Popgen_teaching_code/)

- [refGenie](http://refgenie.databio.org/en/latest/) - reference genome manager. Python.

- [somailer](https://github.com/brentp/somalier) [NGSCheckMate](https://github.com/parklab/NGSCheckMate) - QC, Sample swap check.

- [SNPhylo](https://github.com/thlee/SNPhylo) - A pipeline to generate a phylogenetic tree from huge SNP data. R, shell.

