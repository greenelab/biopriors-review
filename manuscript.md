---
author-meta:
- Jake Crawford
- Casey S. Greene
date-meta: '2019-12-10'
keywords:
- machine learning
- deep learning
- gene networks
- pathway databases
- literature review
lang: en-US
title: Incorporating biological structure into machine learning models in biomedicine
...






<small><em>
This manuscript
([permalink](https://greenelab.github.io/biopriors-review/v/b40606a323ed9e9ede40c79947895186f9717c14/))
was automatically generated
from [greenelab/biopriors-review@b40606a](https://github.com/greenelab/biopriors-review/tree/b40606a323ed9e9ede40c79947895186f9717c14)
on December 10, 2019.
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
     · Funded by National Human Genome Research Institute (R01 HG010067)
  </small>

+ **Casey S. Greene**<br>
    ![ORCID icon](images/orcid.svg){.inline_icon}
    [0000-0001-8713-9213](https://orcid.org/0000-0001-8713-9213)
    · ![GitHub icon](images/github.svg){.inline_icon}
    [cgreene](https://github.com/cgreene)
    · ![Twitter icon](images/twitter.svg){.inline_icon}
    [greenescientist](https://twitter.com/greenescientist)<br>
  <small>
     Department of Systems Pharmacology and Translational Therapeutics, Perelman School of Medicine, University of Pennsylvania, Philadelphia, PA; Childhood Cancer Data Lab, Alex's Lemonade Stand Foundation, Philadelphia, PA
     · Funded by The Gordon and Betty Moore Foundation (GBMF 4552), National Human Genome Research Institute (R01 HG010067), National Cancer Institute (R01 CA237170), Alex's Lemonade Stand Foundation (CCDL)
  </small>



## Abstract {.page_break_before}
In biomedical applications of machine learning, relevant information often has a rich structure that is not easily encoded as real-valued predictors.
Examples of such data include DNA or RNA sequences, gene sets or pathways, gene interaction or coexpression networks, ontologies, and phylogenetic trees.
We highlight recent examples of machine learning models that use structure to constrain model architecture or incorporate structured data into model training.
For machine learning in biomedicine, where sample size is limited and model interpretability is critical, incorporating prior knowledge in the form of structured data can be particularly useful.
The area of research would benefit from performant open source implementations and independent benchmarking efforts.



## Graphical abstract

![Schematic showing the main categories of models incorporating structured biological data covered in this review. The first panel shows an example of a model operating on raw sequence data, the second panel shows a model in which dimension reduction is constrained by the connections in a gene network, and the third panel shows a neural network with structure constrained by a phylogeny or ontology.](images/all_models.svg){.white}


## Introduction

It can be challenging to distinguish signal from noise in biomedical datasets, and machine learning methods are particularly hampered when the amount of available training data is small.
Incorporating biomedical knowledge into machine learning models can reveal patterns in noisy data [@Zm1iJl59] and aid model interpretation [@1HLd4LXUB].
Biological knowledge can take many forms, including genomic sequences, pathway databases, gene interaction networks, and knowledge hierarchies such as the Gene Ontology [@eH3LaU5K].
However, there is often no canonical way to encode these structures as real-valued predictors.
Modelers must creatively decide how to encode biological knowledge that they expect will be relevant to the task.

Biomedical datasets often contain more input predictors than data samples [@1B2Ue9gnv; @JPGgWnoo].
A genetic study may genotype millions of single nucleotide polymorphisms (SNPs) in thousands of individuals, or a gene expression study may profile the expression of thousands of genes in tens of samples.
Thus, it can be useful to include prior information describing relationships between predictors to inform the representation learned by the model.
This contrasts with non-biological applications of machine learning, where one might fit a model on millions of images [@lt4BNUoG] or tens of thousands of documents [@xPLjeqyD], making inclusion of prior information unnecessary.

We review approaches that incorporate external information about the structure of desirable solutions to learn from biomedical data.
One class of commonly used approaches learns a representation that considers the context of each base pair from raw sequence data.
For models that operate on gene expression data or genetic variants, it can be useful to incorporate networks or pathways describing relationships between genes.
We also consider other examples, such as neural network architectures that are constrained based on biological knowledge.

There are many complementary ways to incorporate heterogeneous sources of biomedical data into the learning process, which have been covered elsewhere [@15hcJCu89; @NfKILi8i].
These include feature extraction or representation learning prior to modeling and/or other data integration methods that do not necessarily involve customizing the model itself.


## Sequence models

Early neural network models primarily used hand-engineered sequence features as input to a fully connected neural network [@12aqvAgz6; @17sgPdcMT] (Figure {@fig:sequence_features}).
As convolutional neural network (CNN) approaches matured for image processing and computer vision, researchers leveraged biological sequence proximity similarly.
CNNs are a neural network variant that groups input data by spatial context to extract features for prediction.

The definition of "spatial context" is specific to the input: one might group image pixels that are nearby in 2D space, or genomic base pairs that are nearby in the linear genome.
In this way, CNNs consider context without making strong assumptions about exactly how much context is needed or how it should be encoded; the data informs the encoding.
A detailed description of how CNNs are applied to sequences can be found in Angermueller et al. [@irSe12Sm].

![Contrasting approaches to extracting features from DNA or RNA sequence data. Early models defined features of interest by hand based on prior knowledge about the prediction or clustering problem of interest, such as GC content or sequence melting point, as depicted in the left branch in the figure. Convolutional models, depicted in the right branch, use sequence convolutions to derive features directly from sequence proximity, without requiring quantities of interest to be identified before the model is trained. Red or blue emphasis denotes inputs to the predictive model (either the hand-defined numeric features on the left or the outputs of convolutional filters on the right).](images/sequence_features_revised.svg){#fig:sequence_features .white}

### Applications in regulatory biology

Many early applications of deep learning to biological sequences were in regulatory biology.
Early CNNs for sequence data predicted binding protein sequence specificity from DNA or RNA sequence [@jJHZHWrl], variant effects from noncoding DNA sequence [@2UI1BZuD], and chromatin accessibility from DNA sequence [@2CbHXoFn].

Recent sequence models take advantage of hardware advances and methodological innovation to incorporate more sequence context and rely on fewer modeling assumptions.
BPNet, a CNN that predicts transcription factor binding profiles from DNA sequences, accurately mapped known locations of binding motifs in mouse embryonic stem cells [@miFkgNt5].
BPNet considers 1000 base pairs of context around each position when predicting binding probabilities with a technique called dilated convolutions [@18MdGd1xW], which is particularly important because motif spacing and periodicity can influence binding.
cDeepbind [@NyVCOHer] combines RNA sequences with information about secondary structure to predict RNA binding protein affinities.
Its convolutional model acts on a feature vector combining sequence and structural information, using context for both to inform predictions.
APARENT [@17KNm3K0B] is a CNN that predicts alternative polyadenylation (APA) from a training set of over 3 million synthetic APA reporter sequences.
These diverse applications underscore the power of modern deep learning models to synthesize large sequence datasets.

Models that consider sequence context have also been applied to epigenetic data.
DeepSignal [@phJJvCFv] is a CNN that uses contextual electrical signals from Oxford Nanopore single-molecule sequencing data to predict 5mC or 6mA DNA methylation status.
MRCNN [@NzYX9e9i] uses sequences of length 400, centered at CpG sites, to predict 5mC methylation status.
Deep learning models have also been used to predict gene expression from histone modifications [@126y5dSh0; @rxMdCSQm].
Here, a neural network model consisting of long short-term memory (LSTM) units was used to encode the long-distance interactions of histone marks in both the 3' and 5' genomic directions.
In each of these cases, proximity in the linear genome helped model the complex interactions between DNA sequence and epigenome.

### Applications in variant calling and mutation detection

Identification of genetic variants also benefits from models that include sequence context.
DeepVariant [@YqAWSEkm] applies a CNN to images of sequence read pileups, using read data around each candidate variant to accurately distinguish true variants from sequencing errors.
<!-- could mention GATK4 here which uses a CNN, although nothing has been published? -->
<!-- https://gatkforums.broadinstitute.org/gatk/discussion/10996/deep-learning-in-gatk4 -->
CNNs have also been applied to single molecule (PacBio and Oxford Nanopore) sequencing data [@UsuUETZK], using a different sequence encoding that results in better performance than DeepVariant on single molecule data.
However, many variant calling models still use hand-engineered sequence features as input to a classifier, including current state-of-the-art approaches to insertion/deletion calling [@KDXezqfK; @1ENge146H].
Detection of somatic mutations is a distinct but related challenge to detection of germline variants, and has also recently benefitted from use of CNNs [@t17iSYjA].



## Network- and pathway-based models

Rather than operating on sequences, many machine learning models in biomedicine operate on inputs that lack intrinsic order.
Models may make use of gene expression data matrices from RNA sequencing or microarray experiments in which rows represent samples and columns represent genes.
To account for relationships between genes, one might incorporate known interactions or correlations when making predictions or generating a low-dimensional representation of the data (Figure {@fig:network_models}).
This is comparable to the manner in which sequence context pushes models to consider nearby base pairs similarly.

![The relationships between genes provide structure that can be incorporated into machine learning models. One common approach is to use a network or collection of gene sets to embed the data in a lower-dimensional space, in which genes that are in the same gene sets or that are well-connected in the network have a similar representation in the lower-dimensional space. The embedded data can then be used for classification or clustering tasks. The "x" values in the data table represent gene expression measurements.](images/network_models_revised.svg){#fig:network_models .white}

### Applications in transcriptomics

Models built from gene expression data can benefit from incorporating gene-level relationships.
One form that this knowledge commonly takes is a database of gene sets, which may represent biological pathways or gene signatures for a biological state of interest.
PLIER [@Ki2ij7zE] uses gene set information from MSigDB [@15p5LWIVP] and cell type markers to extract a representation of gene expression data that corresponds to biological processes and reduces technical noise.
The resulting gene set-aligned representation accurately decomposed cell type mixtures.
MultiPLIER [@14rnBunuZ] applied PLIER to the recount2 gene expression compendium [@6SPTvFXq] to develop a model that shares information across multiple tissues and diseases, including rare diseases with limited sample sizes.
PASNet [@1Bb8CyeTY] uses MSigDB to inform the structure of a neural network for predicting patient outcomes in glioblastoma multiforme (GBM) from gene expression data.
This approach aids interpretation, as pathway nodes in the network with high weights can be inferred to correspond to certain pathways in GBM outcome prediction.

Gene-level relationships can also be represented with networks.
Network nodes typically represent genes and real-valued edges may represent interactions or correlations between genes, often in a tissue or cell type context of interest.
netNMF-sc [@17fvHtbrH] incorporates coexpression networks [@3VYPTgXw] as a smoothing term for dimension reduction and dropout imputation in single-cell gene expression data.
The coexpression network improves performance for identifying cell types and cell cycle marker genes, as compared to using raw gene expression or other single-cell dimension reduction methods.
Combining gene expression data with a network-derived smoothing term also improved prediction of patient drug response in acute myeloid leukemia [@LLInUBEI] and detection of mutated cancer genes [@1BbuXJuIl].
PIMKL [@12cJO5Pse] combines network and pathway data to predict disease-free survival from breast cancer cohorts.
This method takes as input both RNA-seq gene expression data and copy number alteration data, but can also be applied to gene expression data alone.

Gene regulatory networks can also augment models for gene expression data.
These networks describe how the expression of genes is modulated by biological regulators such as transcription factors, microRNAs, or small molecules.
creNET [@g8OoyIPj] integrates a gene regulatory network, derived from STRING [@q9Fhy8eq], with a sparse logistic regression model to predict phenotypic response in clinical trials for ulcerative colitis and acute kidney rejection.
The gene regulatory information allows the model to identify the biological regulators associated with the response, potentially giving mechanistic insight into differential clinical trial response.
GRRANN [@19wuAzYvo], which was applied to the same data as creNET, uses a gene regulatory network to inform the structure of a neural network.
Several other methods [@owp8L957; @e4tSAJkK] have also used gene regulatory network structure to constrain the structure of a neural network, reducing the number of parameters to be fit and facilitating interpretation.

### Applications in genetics

Approaches that incorporate gene set or network structure into genetic studies have a long history [@PuUYD4zV; @MPgqEKuv].
Recent applications include expression quantitative trait loci (eQTL) mapping studies, which aim to identify associations between genetic variants and gene expression.
netReg [@1AuHIuXiR] implements a graph-regularized dual LASSO algorithm for eQTL mapping [@9SBhyy2x] in a publicly available R package.
This model smooths regression coefficients simultaneously based on networks describing associations between genes (target variables in the eQTL regression model) and between variants (predictors in the eQTL regression model).
eQTL information is also used in conjunction with genetic variant information to predict phenotypes, in an approach known as Mendelian randomization (MR).
In [@m7CdzEUG], a smoothing term derived from a gene regulatory network is used in an MR model.
The model with the network smoothing term, applied to a human liver dataset, more robustly identifies genes that influence enzyme activity than a network-agnostic model.
As genetic datasets grow, we expect that researchers will continue to develop models that leverage gene set and network databases.


## Other models incorporating biological structure

Knowledge about biological entities is often organized in an ontology, which is a directed graph that encodes relationships between entities (see Figure {@fig:ontology_models} for a visual example).
The Gene Ontology (GO) [@eH3LaU5K] describes the relationships between cellular subsystems and other attributes describing proteins or genes.
DCell [@qQP20moO] uses GO to inform the connectivity of a neural network predicting the effects of gene deletions on yeast growth.
DCell performs comparably to an unconstrained neural network for this task.
Additionally, it is easier to interpret: a cellular subsystem with high neuron outputs under a particular gene deletion can be inferred to be strongly affected by the gene deletion, providing a putative genotype-phenotype association.
DeepGO [@TIQTmEOG] uses a similar approach to predict protein function from amino acid sequence with a neural network constrained by the dependencies of GO.
However, a follow-up paper by the same authors [@Cf5duPBD] showed that this hierarchy-aware approach can be outperformed by a hierarchy-naive CNN, which uses only amino acid sequence and similarity to labeled training set proteins.
This suggests a tradeoff between interpretability and predictive accuracy for protein function prediction.

![Directed graph-structured data, such as an ontology or phylogenetic tree, can be incorporated into machine learning models. Here, the connections in the neural network used to predict a set of labels parallel those in the tree graph. This type of constraint can also be useful in model interpretation: for example, if the nodes in the right tree branch have high neuron outputs for a given sample, then the subsystem encoded in the right branch of the tree graph is most likely important in making predictions for that sample. The "x" values in the data table represent gene expression measurements.](images/ontology_models_revised.svg){#fig:ontology_models .white}

Phylogenetic trees, or hierarchies describing the evolutionary relationships between species, can be useful for a similar purpose.
glmmTree [@uQ5z1fAc] uses a phylogenetic tree describing the relationship between microorganisms to improve predictions of age based on gut microbiome data.
The same authors combine a similar phylogeny smoothing strategy with sparse regression to model caffeine intake and smoking status based on microbiome data [@VvllhwW1].
Phylogenetic trees can also describe the relationships between subclones of a tumor, which are fundamental to understanding cancer evolution and development.
Using a tumor phylogeny inferred from copy number aberration (CNA) sequencing data as a smoothing term for deconvolving tumor subclones provided more robust predictions than a phylogeny-free model [@VzWLIpJn].
The tree structure of the phylogeny and the subclone mixture model are fit jointly to the CNA data.

Depending on the application, other forms of structure or prior knowledge can inform predictions and interpretation of the model's output.
CYCLOPS [@1Dk2FpWk6] uses a circular node autoencoder [@15jrSlkXD] to order periodic gene expression data and estimate circadian rhythms.
The authors validated the method by correctly ordering samples without temporal labels and identifying genes with known circadian expression.
They then applied it to compare gene expression in normal and cancerous liver biopsies, identifying drug targets with circadian expression as candidates for chronotherapy.
NetBiTE [@kKiwlzZq] uses drug-gene interaction information from GDSC [@lJFuND9F], in addition to protein interaction data, to build a tree ensemble model with splits that are biased toward high-confidence drug-gene interactions.
The model predicts sensitivity to drugs that inhibit critical signaling pathways in cancer, showing improved predictive performance compared to random forests, another commonly used tree ensemble model.


## Conclusions and future directions

As the quantity and richness of biomedical data has increased, sequence repositories and interaction databases have expanded and become more robust.
This raises opportunities to integrate these resources into the structure of machine learning models.
Going forward, there is an outstanding need for benchmarks comparing these approaches across diverse datasets and prediction problems, along the lines of the evaluation in [@LL5fLwtS] but updated and expanded to include recent methods and applications.
Improved benchmarking should lead to a better understanding of which dataset characteristics align with which approaches.

Many methods described in this review have open-source implementations available; however, increased availability of performant and extensible implementations of these models and algorithms would facilitate further use and development.
In the future, we foresee that incorporating structured biomedical data will become commonplace for improving model interpretability and boosting performance when sample size is limited.



## Acknowledgements

The authors would like to thank Daniel Himmelstein for a critical reading of the manuscript and helpful discussion.


## Reference Annotations

[**] Annotation for BPNet [@miFkgNt5]:

This paper describes BPNet, a neural network for predicting transcription
factor (TF) binding profiles from raw DNA sequence. The model is able to
accurately infer the spacing and periodicity of pluripotency-related TFs in
mouse embryonic stem cells, leading to an improved understanding of the motif
syntax of combinatorial TF binding in cell development.

[*] Annotation for cDeepbind [@NyVCOHer]:

cDeepbind is a neural network model for predicting RNA binding protein (RBP)
specificity from RNA sequence and secondary structure information. The authors
show that this combined approach provides an improvement over previous models
that use only sequence information.

[*] Annotation for DeepDiff [@rxMdCSQm]:

DeepDiff uses a long short-term memory neural network to predict differential
gene expression from the spatial structure of histone modification
measurements. The network has a multi-task objective, enabling gene expression
predictions to be made simultaneously in multiple cell types.

[**] Annotation for DeepVariant [@YqAWSEkm]:

This paper describes DeepVariant, a neural network model for distinguishing
true genetic variants from errors in next-generation DNA sequencing data. The
model adapts techniques from the image processing community to fit a model on
images of read pileups around candidate variants, using information about the
sequence around the candidate variant site to make predictions about the true
genotype at the site.

[**] Annotation for PLIER [@Ki2ij7zE]:

This paper describes a "pathway-level information extractor" (PLIER), a method
for reducing the dimension of gene expression data in a manner that aligns with
known biological pathways or informative gene sets. The method can also reduce
the effects of technical noise. The authors show that PLIER can be used to
improve cell type inference and as a component in eQTL studies.

[**] Annotation for netNMF-sc [@17fvHtbrH]:

netNMF-sc is a dimension reduction method that uses network information to
"smooth" a matrix factorization of single-cell gene expression data, such that
genes that are connected in the network have a similar low-dimensional
representation. Inclusion of network information is particularly useful when
analyzing single-cell expression data, due to its ability to mitigate "dropouts"
and other sources of variability that are present at the single cell level.

[*] Annotation for Attribution Priors [@LLInUBEI]:

This paper describes "model attribution priors", or a method for constraining
a machine learning model's behavior during training with prior beliefs or
expectations about the data or problem structure. As an example of this concept,
the authors show that incorporation of network data improves the performance of
a model for drug response prediction in acute myeloid leukemia.

[*] Annotation for PIMKL [@12cJO5Pse]:

In this paper, the authors present an algorithm for combining gene expression
and copy number data with prior information, such as gene networks and pathways
or gene set annotations, to predict survival in breast cancer. The weights
learned by the model are also interpretable, providing a putative set of
explanatory features for the prediction task.

[**] Annotation for creNET [@g8OoyIPj]:

This work describes creNET, a regression model for gene expression data that
uses information about gene regulation to differentially weight or penalize
gene sets that are co-regulated. The authors show that the model can be used to
predict phenotype from gene expression data in clinical trials. The model also
provides interpretable weights for each gene regulator.

[**] Annotation for DCell [@qQP20moO]:

This paper presents DCell, a neural network model for prediction of yeast
growth phenotype from gene deletions. The structure of the neural network is
constrained by the relationships encoded in the Gene Ontology (GO), enabling
predictions for a given input to be interpreted based on the subsystems of GO
that they activate. Thus, the neural network can be seen as connecting genotype
to phenotype.

[*] Annotation for DeepGO [@TIQTmEOG]:

Here, the authors describe a method for predicting protein function from amino
acid sequence, incorporating the dependency structure of the Gene Ontology (GO)
into their neural network used for prediction. Using the GO information
provides a performance improvement over similar models that do not incorporate
this information.

[**] Annotation for NetBITE [@kKiwlzZq]:

This paper describes a method for using prior knowledge about drug targets to
inform the structure of a tree ensemble model, used for predicting IC50 drug
sensitivity from gene expression data. The model also uses a protein
interaction network to "smooth" the gene weights, such that genes that are
related in the network will have a similar influence on predictions.



## References {.page_break_before}

<!-- Explicitly insert bibliography here -->
<div id="refs"></div>
