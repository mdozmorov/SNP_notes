# Notes on SNP-related tools and genome variation analysis

Issues with suggestions and pull requests are welcome!

# Table of content

* [SNP callers](#snp-callers)
  * [Deep learning SNP callers](#deep-learning-snp-callers)

## SNP callers

### Deep learning SNP callers

- `NeuSomatic` - Convolutional neural network (9 layers) for somatic mutations calling. Reads within 7bp window around candidate SNP are extracted, realigned, summarized into matrices for tumor-normal samples,  used for classifying mutation type, length, position. Tested on GiB samples and on DREAM datasets. Comparison with other SNP callers. https://github.com/bioinform/neusomatic
    - Sahraeian, Sayed Mohammad Ebrahim, Ruolin Liu, Bayo Lau, Marghoob Mohiyuddin, and Hugo Y. K. Lam. “Deep Convolutional Neural Networks for Accurate Somatic Mutation Detection,” September 4, 2018. https://doi.org/10.1101/393801.

