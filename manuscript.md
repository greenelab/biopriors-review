---
author-meta:
- Jake Crawford
- Jane Roe
date-meta: '2019-10-09'
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
([permalink](https://greenelab.github.io/biopriors-review/v/7e8a3a2c7f58bf8fee3c1547870f54f98141fdb1/))
was automatically generated
from [greenelab/biopriors-review@7e8a3a2](https://github.com/greenelab/biopriors-review/tree/7e8a3a2c7f58bf8fee3c1547870f54f98141fdb1)
on October 9, 2019.
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
This stands in contrast to non-biological applications of machine learning, where one might fit a model on millions of images [@lt4BNUoG] or tens of thousands of documents [@xPLjeqyD], making inclusion of prior information unnecessary.

In this review, we survey approaches to learning models from biomedical data that incorporate external information about the structure of desirable solutions.
One class of commonly used approaches involves using raw sequence data to learn a representation that considers the context of each base pair.
For models that operate on gene expression data or genetic variants, it can be useful to incorporate networks or pathways describing relationships between genes.
We also consider other examples in this review, such as neural network architectures that are constrained based on biological knowledge.

In addition to the methods surveyed here, there are many complementary options for utilizing heterogeneous sources of biomedical data.
These include feature extraction or representation learning prior to modeling, and/or other data integration methods that do not necessarily involve customizing the model itself.
Many of these methods have been covered extensively elsewhere [@15hcJCu89; @NfKILi8i].



## Sequence models

Early neural network models primarily used hand-engineered sequence features as input to a fully connected neural network [@12aqvAgz6; @17sgPdcMT].
As convolutional neural network (CNN) approaches matured for image processing and computer vision, researchers were able to use similar ideas to leverage biological sequence proximity in modeling.
CNNs are a neural network variant in which input data are grouped by spatial context to extract features for prediction.
The definition of "spatial context" is specific to the input.
For example, for images pixels that are nearby in 2D space might be grouped, or for genomic sequences base pairs that are nearby in the linear genome might be grouped.
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

Models that consider sequence context have also been applied to epigenetic data.
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
The authors show that using a coexpression network to extract a low-dimensional representation improves performance for cell type identification and identification of cell cycle marker genes, as compared to using raw gene expression or other single-cell dimension reduction methods.
Combining gene expression data with a network-derived smoothing term has also been shown to improve prediction of patient drug response in acute myeloid leukemia [@LLInUBEI] and detection of mutated cancer genes [@1BbuXJuIl].
PIMKL [@12cJO5Pse] combines network and pathway data to predict disease-free survival from breast cancer cohorts.
This method takes as input both RNA-seq gene expression data and copy number alteration data, but can be applied to gene expression data alone as well.

Gene regulatory networks can also augment models for gene expression data.
These networks describe how the expression of genes is modulated by biological regulators such as transcription factors, microRNAs, or small molecules.
creNET [@g8OoyIPj] integrates a gene regulatory network, derived from STRING [@q9Fhy8eq], with a sparse logistic regression model to predict phenotypic response in clinical trials for ulcerative colitis and acute kidney rejection based on gene expression data.
The gene regulatory information allows the model to identify the biological regulators that are associated with the response, potentially giving mechanistic insight into differential clinical trial response.
GRRANN [@19wuAzYvo] uses a gene regulatory network to inform the structure of a neural network, applying it to the same clinical trial data as creNET.
Several other methods [@owp8L957; @e4tSAJkK] have also used gene regulatory network structure to constrain the structure of a neural network, reducing the number of parameters to be fit by the network and facilitating interpretation of network predictions.

### Applications in genetics

Approaches to incorporating gene set or network structure into genetic studies have a long history (see, e.g. [@PuUYD4zV; @MPgqEKuv]).
Recent applications of these methods include expression quantitative trait loci (eQTL) mapping studies, which aim to identify associations between genetic variants and gene expression.
netReg [@1AuHIuXiR] implements the graph-regularized dual LASSO algorithm for eQTL mapping described in [@9SBhyy2x] in a publicly available R package, based on an efficient C++ backend.
This model smooths regression coefficients simultaneously based on networks describing associations between genes (target variables in the eQTL regression model) and between variants (predictors in the eQTL regression model).
eQTL information can also be used in conjunction with genetic variant information to predict phenotypes, in an approach known as Mendelian randomization (MR).
In [@m7CdzEUG], a smoothing term derived from a gene regulatory network is used in an MR model.
The model with the network smoothing term, applied to a human liver data set, more robustly identifies genes that influence enzyme activity than a network-agnostic model.
As genetic datasets become larger, it is likely that researchers will continue to develop models that leverage gene set and network databases.


## Other models incorporating biological structure

Knowledge about biological entities is often organized in an ontology, which is a directed graph that encodes the relationships between entities.
The Gene Ontology (GO) [@eH3LaU5K] describes the relationships between cellular subsystems and other attributes describing proteins or genes.
DCell [@qQP20moO] uses GO to inform the connectivity of a neural network predicting the effects of gene deletions on yeast growth.
The authors show that DCell performs comparably to an unconstrained neural network for this task.
Additionally, it has the advantage of interpretability: a cellular subsystem with high neuron output values under a particular gene deletion can be inferred to be strongly affected by the gene deletion, providing a putative genotype-phenotype association.
DeepGO [@TIQTmEOG] uses a similar approach to predict protein function from amino acid sequence, using a neural network constrained by the dependencies of GO.
However, a follow-up paper by the same authors [@Cf5duPBD] showed that this hierarchy-aware approach can be outperformed by a hierarchy-naive CNN, which uses only raw amino acid sequence and similarity to labeled proteins in the training set.
This suggests a potential tradeoff between interpretability and predictive accuracy in the context of protein function prediction.

Phylogenetic trees, or hierarchies describing the evolutionary relationships between species, can be useful for a similar purpose.
glmmTree [@uQ5z1fAc] uses a phylogenetic tree describing the relationship between microorganisms to improve predictions of age based on gut microbiome sequencing data.
The same authors combine a similar phylogeny smoothing strategy with a sparse regression model in [@VvllhwW1], showing its utility to model caffeine intake and smoking status based on microbiome data.
Phylogenetic trees can also be useful in describing the relationships between subclones of a tumor, which are fundamental to understanding cancer evolution and development.
In [@VzWLIpJn], the authors use a tumor phylogeny inferred from copy number aberration (CNA) sequencing data as a smoothing term in an algorithm for deconvolving tumor subclones, providing more robust predictions than a phylogeny-free model.
The tree structure of the phylogeny and the subclone mixture model are fit jointly to the CNA data.

Depending on the application, incorporating other forms of structure or prior knowledge can help to inform predictions and interpretation of the model's output.
CYCLOPS [@1Dk2FpWk6] uses a circular node autoencoder [@15jrSlkXD] to order periodic gene expression data and measure circadian rhythms.
The authors validated the method by correctly ordering samples without temporal labels and identifying genes with known circadian expression.
They then applied it to compare gene expression in normal and cancerous liver biopsies, identifying drug targets with circadian expression as candidates for chronotherapy.
NetBiTE [@kKiwlzZq] uses drug-gene interaction information from GDSC [@lJFuND9F], in addition to protein interaction data, to build a tree ensemble model with splits that are biased toward high-confidence drug-gene interactions.
The model is used to predict sensitivity to drugs that inhibit critical signaling pathways in cancer, showing improved predictive performance compared to random forests, another commonly used tree ensemble model.


## Conclusions and future directions

As the quantity and richness of biomedical data has increased, resources such as sequence repositories and interaction databases have expanded and become more robust.
This has created unique opportunities for integrating these resources into machine learning models in a way that considers their structure.
Going forward, there is an outstanding need for benchmarks comparing these approaches across diverse datasets and prediction problems, along the lines of the evaluation in [@LL5fLwtS] but updated and expanded to include recent methods and applications.
Ideally, improved benchmarking will lead to a better understanding of which datasets can benefit from which approaches, guiding application of similar models to new datasets.
Many of the methods described in this review have open-source implementations available; however, increased availability of performant and extensible implementations of the models and algorithms described in this review would also facilitate further use and development.
In the future, we foresee that incorporating structured biomedical data will become commonplace for improving model interpretability and boosting performance when sample size is limited.



## References {.page_break_before}

<!-- Explicitly insert bibliography here -->
<div id="refs"></div>
