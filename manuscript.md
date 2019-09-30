---
author-meta:
- Jake Crawford
- Jane Roe
date-meta: '2019-09-30'
keywords:
- machine-learning
- deep-learning
- gene-networks
- pathway-databases
- literature-review
lang: en-US
title: Incorporating biological structure into machine learning models in biomedicine
...






<small><em>
This manuscript
([permalink](https://greenelab.github.io/biopriors-review/v/f72bcb4019d50b6a1584a2f0bf0cf3bef96eb8b6/))
was automatically generated
from [greenelab/biopriors-review@f72bcb4](https://github.com/greenelab/biopriors-review/tree/f72bcb4019d50b6a1584a2f0bf0cf3bef96eb8b6)
on September 30, 2019.
</em></small>

## Authors



+ **Jake Crawford**<br>
    ![ORCID icon](images/orcid.svg){.inline_icon}
    [0000-0001-6207-0782](https://orcid.org/0000-0001-6207-0782)
    · ![GitHub icon](images/github.svg){.inline_icon}
    [jjc2718](https://github.com/jjc2718)
    · ![Twitter icon](images/twitter.svg){.inline_icon}
    [jjc2718](https://twitter.com/jjc2718)<br>
  <small>
     Graduate Group in Genomics and Computational Biology, Perelman School of Medicine, University of Pennsylvania, Philadelphia, PA; Department of Systems Pharmacology and Translational Therapeutics, Perelman School of Medicine, University of Pennsylvania, Philadelphia, PA
     · Funded by Grant XXXXXXXX
  </small>

+ **Jane Roe**<br>
    ![ORCID icon](images/orcid.svg){.inline_icon}
    [XXXX-XXXX-XXXX-XXXX](https://orcid.org/XXXX-XXXX-XXXX-XXXX)
    · ![GitHub icon](images/github.svg){.inline_icon}
    [janeroe](https://github.com/janeroe)<br>
  <small>
     Department of Something, University of Whatever; Department of Whatever, University of Something
  </small>



## Abstract {.page_break_before}




## Introduction

When applying machine learning techniques to biomedical datasets, it can be challenging to distinguish signal from noise, particularly in the presence of limited amounts of data.
Biological knowledge can take many forms, including genomic sequences, pathway databases, gene interaction networks, and knowledge hierarchies such as the Gene Ontology [@eH3LaU5K].
Incorporating these resources in machine learning models can be helpful in identifying patterns in noisy data [@Zm1iJl59] and in interpreting model predictions [@1HLd4LXUB].
However, there is often no canonical way to encode these structures as real-valued predictors.
This means modelers must be creative when deciding how to encode biological knowledge that they expect will be relevant to the task in models.

Biomedical datasets often contain more input predictors than data samples [@1B2Ue9gnv; @JPGgWnoo].
For example, a genetic study may genotype millions of single nucleotide polymorphisms (SNPs) in hundreds of patients, or a gene expression study may profile the expression of thousands of genes in only a handful of samples.
Thus, it can be useful to include prior information describing the relationships between the predictors to inform the representation learned by the model.
This stands in contrast to non-biological applications of machine learning, where one might fit a model on millions of images [@lt4BNUoG] or tens of thousands of documents [@hIUpKMFp], making inclusion of prior information unnecessary.

In this review, we survey approaches to learning models from biomedical data that incorporate external information about the structure of desirable solutions.
One class of commonly used approaches involves using raw sequence data to learn a representation that considers the context of each base pair.
For models that operate on gene input, such as gene expression data or genetic variants, it can be useful to incorporate networks or pathways describing relationships between genes.
We also consider other examples in this review, such as neural network architectures that are constrained based on biological knowledge.



## Sequence models

Early neural network models primarily used hand-engineered sequence features as input to a fully connected neural network [@12aqvAgz6; @17sgPdcMT].
As convolutional neural network (CNN) approaches matured for image processing and computer vision, researchers were able to use similar ideas to leverage biological sequence proximity in modeling.
CNNs are a neural network variant in which input data are grouped by spatial context to extract features for prediction.
The definition of "spatial context" is specific to the input.
For example, for images pixels that are nearby in 2D space might be grouped, or for genomic sequences base pairs that are nearby in the linear genome might be grouped.

These approaches work by first encoding input into a numeric matrix (for DNA "one-hot" encoding is often used: A=[1,0,0,0], C=[0,1,0,0], G=[0,0,1,0], T=[0,0,0,1]).
They then apply spatial filters to the encoded input, which are adjusted based on the data during model training.
In this way, CNNs are able to consider context without making strong assumptions about exactly how much context is needed or how it should be encoded; the data informs the encoding.
A detailed description of how CNNs are applied to sequence data can be found in [@irSe12Sm].

### Applications in regulatory biology

Many of the first applications of deep learning to biological sequence data were in regulatory biology.
Example early applications of CNNs on sequence data include prediction of binding protein sequence specificity from DNA or RNA sequence [@jJHZHWrl], prediction of variant effects from noncoding DNA sequence [@2UI1BZuD], and prediction of chromatin accessibility from DNA sequence [@2CbHXoFn].

Recent sequence models take advantage of hardware advances and methodological innovation to incorporate more sequence context and rely on fewer modeling assumptions.
BPNet, a CNN used to predict transcription factor binding profiles from raw DNA sequences, was able to accurately map known locations of binding motifs in mouse embryonic stem cells [@miFkgNt5].
The BPNet model considers 1000 base pairs of context around each position when predicting binding probabilities, using a technique called dilated convolutions [@18MdGd1xW] to make a large input field feasible.
This context is particularly important because motif spacing and periodicity can influence protein function.
cDeepbind [@NyVCOHer] combines RNA sequences with information about secondary structure to predict RNA binding protein affinities.
Its convolutional model acts on a feature vector combining sequence and structural information, simultaneously using context for both to inform predictions.
APARENT [@17KNm3K0B] is a CNN used to predict alternative polyadenylation (APA) from a training set of over 3 million synthetic APA reporter sequences.
These diverse applications underscore the power of modern deep learning models to synthesize large sequence datasets.

Models that consider sequence context have also been applied to impute and make predictions from epigenetic data.
DeepSignal [@phJJvCFv] is a CNN that uses contextual electrical signals from Oxford Nanopore single-molecule sequencing data to predict 5mC or 6mA DNA methylation status.
MRCNN [@NzYX9e9i] uses sequences of length 400, centered at CpG sites, to predict 5mC methylation status.
Deep learning models have also been used to predict gene expression from histone modifications [@126y5dSh0; @rxMdCSQm].
Here, a neural network model consisting of long short-term memory (LSTM) units was used to encode the long-distance interactions of histone marks in both the 3' and 5' genomic directions.
In each of these cases, proximity in the linear genome proved useful for modeling the complex interactions between DNA sequence and epigenome.

### Applications in variant calling and mutation detection

Identification of genetic variants can also benefit from models that take into account sequence context.
DeepVariant [@YqAWSEkm] applies a CNN to images of sequence read pileups, using read data around each candidate variant to accurately distinguish true variants from sequencing errors.
<!-- could mention GATK4 here which uses a CNN, although nothing has been published? -->
<!-- https://gatkforums.broadinstitute.org/gatk/discussion/10996/deep-learning-in-gatk4 -->
CNNs have also been applied to single molecule (PacBio and Oxford Nanopore) sequencing data [@UsuUETZK], using a different sequence encoding that results in better performance than DeepVariant on single molecule data.
However, many variant calling models still use hand-engineered sequence features as input to a classifier, including current state-of-the-art approaches to insertion/deletion calling [@KDXezqfK; @1ENge146H].
Detection of somatic mutations is a distinct but related challenge to detection of germline variants, and has also recently benefitted from use of CNN models [@t17iSYjA].

### Applications in CRISPR guide selection

With the recent rise in popularity of CRISPR gene editing technology, sequence models have proven useful in improving design of single-guide RNAs (sgRNAs).
To select the best sgRNA from multiple possibilities, one might be interested in balancing on-target efficiency with likelihood of off-target effects, both of which depend on the sgRNA sequence and its genomic context.
Early models for prediction of on-target cleavage efficiency [@UBa1ZsWW] and off-target effects [@YL6D0uxU] used hand-engineered sequence features, in addition to other RNA-specific information such as thermodynamic sequence properties.
More recently, CNN-based models have demonstrated improved performance using data-derived sequence features.
CNNs have been successfully applied to CRISPR-Cas9 on-target efficiency prediction [@12MXxroFF; @GwBnBz2w], CRISPR-Cas12a (Cpf1) on-target efficiency prediction [@65fi127i], and CRISPR-Cas9 off-target effect prediction [@11czNQNnW; @12MXxroFF].
In each case, the authors show that using a CNN to operate on raw sequence data improves sgRNA design relative to models that use hand-engineered sequence features as input.


## Network- and pathway-based models

Rather than operating on raw DNA sequence, many machine learning models in biomedicine operate on inputs without an intrinsic order.
For instance, models may make use of gene expression data matrices (e.g. from RNA sequencing or microarray experiments), where each row represents a sample and each column a gene.
When modeling data that is indexed by gene, one might incorporate knowledge that describes the relationships or correlations between genes, in an effort to take these relationships into account when making predictions or generating a low-dimensional representation of the data.
This is comparable to the manner in which sequence context encourages models to consider nearby base pairs similarly.

### Applications in transcriptomics

Models that take gene expression data as input can benefit from incorporating gene-level relationships.
One form that this knowledge commonly takes is a database of gene sets, which may represent biological pathways or gene signatures related to a biological state of interest.
PLIER [@Ki2ij7zE] uses gene set information from MSigDB [@15p5LWIVP] and cell type markers to extract a representation of whole blood gene expression data that corresponds to biological processes and reduces technical noise.
The resulting gene set-aligned representation was used to perform accurate cell type mixture decomposition.
MultiPLIER [@14rnBunuZ] applied the PLIER framework to the recount2 gene expression compendium [@6SPTvFXq] to develop a model that shares information across multiple tissues and diseases, including rare diseases with limited sample sizes.
PASNet [@1Bb8CyeTY] uses data from MSigDB to inform the structure of a neural network for predicting patient outcomes in glioblastoma multiforme (GBM) from gene expression data.
This approach has the added benefit of straightforward interpretation, as pathway nodes in the network having high weights can be inferred to correspond to important pathways in GBM outcome prediction.

Alternatively, gene-level relationships can take the form of a network.
Nodes in these networks typically represent genes, and real-valued edges in these networks may represent interactions or correlations between genes, often in a tissue or cell type context of interest.
netNMF-sc [@17fvHtbrH] incorporates coexpression networks [@3VYPTgXw] as a smoothing term to perform dimension reduction and impute dropouts in single-cell gene expression data.
The authors show that using a coexpression network to extract a low-dimensional representation increases performance for cell type identification and identification of cell cycle marker genes, both as compared to using raw gene expression data and as compared to other single-cell dimension reduction methods.
Combining gene expression data with a network-derived smoothing term has also been shown to improve performance at predicting patient drug response in acute myeloid leukemia [@LLInUBEI] and at detecting mutated cancer genes [@1BbuXJuIl].
PIMKL [@12cJO5Pse] combines network and pathway data to predict disease-free survival from breast cancer cohorts.
This method takes as input both RNA-seq gene expression data and copy number alteration data, but can be applied to gene expression data alone as well.

Gene regulatory networks can also augment models for gene expression data.
These networks describe how the expression of genes is modulated by biological regulators such as transcription factors, microRNAs, or small molecules.
creNET [@g8OoyIPj] integrates a gene regulatory network, derived from STRING [@q9Fhy8eq], with a sparse logistic regression model to predict phenotypic response in clinical trials for ulcerative colitis and acute kidney rejection based on gene expression data.
The gene regulatory information allows the model to identify the biological regulators that are associated with the response, potentially giving mechanistic insight into differential clinical trial response.
GRRANN [@19wuAzYvo] uses a gene regulatory network to inform the structure of a neural network, applying it to the same clinical trial data as creNET.
Several other methods [@owp8L957; @e4tSAJkK] have also used gene regulatory network structure to constrain the structure of a neural network, reducing the number of parameters to be fit by the network and facilitating interpretation of network predictions.

### Applications in genetics

* 10.1093/bioinformatics/bty429 (HINT PPI)

1. Genotype-phenotype associations (GWAS/eQTL)
   * 10.1093/bioinformatics/btu293
   * 10.1093/bioinformatics/btx677 (multivariate version of ^)
   * 10.1111/biom.13072 (phenotype prediction using variant and GE data)
2. Using pathways as a prior to study SNP-SNP or gene-gene interactions
 (not sure if this is truly a constraint on the model, more just a
  strategy to reduce the number of hypothesis tests)
   * 10.1101/182741, 10.1371/journal.pgen.1006973 (pathways + GWAS data -> GIs)
   * 10.1371/journal.pgen.1006516 (pathways + GWAS data -> GIs)
3. Using SNP-SNP or gene-gene networks to reduce hypothesis testing
  burden for detecting GIs
   * can find some examples if we decide this is a direction worth going


## Other constrained models

### Ontology-constrained models

1. Phenotype prediction
    * 10.1038/nmeth.4627 (gene ontology -> NN structure, to predict effects of mutations on growth in yeast)
2. Function prediction
    * 10.1093/bioinformatics/btx624 (GO hierarchy relationships -> NN structure, for function prediction)
    * 10.1093/bioinformatics/btx252 (PPI network + tissue ontology -> function prediction)

### Otherwise constrained models

1. Cell cycle information
    * 10.1038/nbt.3102
    * 10.1101/526848
   (in principle the denoising method could be generalized to other gene sets, but here they used cell cycle-relevant gene sets and emphasized the utility of this)
2. Circadian rhythms: 10.1073/pnas.1619320114 (circular node autoencoder for modeling periodic gene expression)
3. TAD/3D chromatin structure information? Can probably find some examples of this


## Conclusions

1. What is outside of the scope of this review?
   (this can also go in introduction?)
   a. Biological "constraints" vs. feature selection or feature extraction
      from heterogeneous biological data (e.g. network embedding approaches)
       * Example: one could use a network-based feature extraction method
         (e.g. Node2Vec) to convert each gene in a PPI network into a set of
         real-valued features, then use those + gene expression as input to
         a model
       * For purposes of keeping this review short enough, I'm trying to
         stay away from papers like ^, but still unclear to me where exactly
         the line should be drawn. Almost any ML model that operates on
         sequence data can be viewed as having a feature extraction component,
         for example.
   b. Could kind-of consider many single-cell dimension reduction methods
      biologically constrained (e.g. dropout/zero inflation modeling
      approaches, etc), but this is way too broad for this review - maybe
      refer the reader to other recent reviews of these methods.
   c. Could also consider omics integration methods (combining, for example,
      gene expression and epigenetic data) to be biologically constrained, but
      we refer the reader to 10.1016/j.inffus.2018.09.012 for further detail
      on these methods.



## References {.page_break_before}

<!-- Explicitly insert bibliography here -->
<div id="refs"></div>
