# Notes on SNP-related tools and genome variation analysis

Issues with suggestions and pull requests are welcome!

# Table of content

* [SNP tools](#snp-tools)
* [SNP callers](#snp-callers)
  * [Deep learning SNP callers](#deep-learning-snp-callers)


## SNP tools

- `PyVCF` - A Variant Call Format Parser for Python. https://pyvcf.readthedocs.io/en/latest/. Has `vcf_melt` tool to reformat a VCF into long format.


## SNP callers

- Xu, Chang. “A Review of Somatic Single Nucleotide Variant Calling Algorithms for Next-Generation Sequencing Data.” Computational and Structural Biotechnology Journal 16 (2018): 15–24. https://doi.org/10.1016/j.csbj.2018.01.003. - Overview of 46 somatic Single Nucleotide Variant (SNV) caller tools. Pre-processing, variant evaluation, and post-filtering steps. Four categories of algorithms, description of each, and the corresponding tools: matched tumor-normal (position-, haplotype-, mathine learning-based methods, [Table 1](https://www.sciencedirect.com/science/article/pii/S2001037017300946?via%3Dihub#t0005)), single-sample ([Table 2](https://www.sciencedirect.com/science/article/pii/S2001037017300946?via%3Dihub#t0010), some offer somatic-germline classification), UMI-based (UMI technology, [Figure 1](https://www.sciencedirect.com/science/article/pii/S2001037017300946?via%3Dihub#f0005), [Table 3](https://www.sciencedirect.com/science/article/pii/S2001037017300946?via%3Dihub#t0015)), and RNA-seq (Technology, issues, [Table 4](https://www.sciencedirect.com/science/article/pii/S2001037017300946?via%3Dihub#t0020)) variant calling. Benchmarking using tools for generating synthetic reads, spike-ins, GiAB, melanoma-normal samples, performance evaluation metrics. Issues in representing complex variants and tools for variant normalization. Deep neural network-based algorithms perform best.

- `MutationSeq` - somatic SNV detection from tumor-normal pairs. http://compbio.bccrc.ca/software/mutationSeq/




### Deep learning SNP callers

- `NeuSomatic` - Convolutional neural network (9 layers) for somatic mutations calling. Reads within 7bp window around candidate SNP are extracted, realigned, summarized into matrices for tumor-normal samples,  used for classifying mutation type, length, position. Tested on GiB samples and on DREAM datasets. Comparison with other SNP callers. https://github.com/bioinform/neusomatic
    - Sahraeian, Sayed Mohammad Ebrahim, Ruolin Liu, Bayo Lau, Marghoob Mohiyuddin, and Hugo Y. K. Lam. “Deep Convolutional Neural Networks for Accurate Somatic Mutation Detection,” September 4, 2018. https://doi.org/10.1101/393801.

