# Thesis outline
1. Preamble
2. Introduction
3. Statistical model of methylC-seq data
4. Standard analysis of methylC-seq data
5. m-tuples
6. Correlation of $\beta$-values
7. Simulation model: Part I
8. Simulation model: Part II
9. Data analysis - EPISCOPE
10. Concluding Remarks
11. Bibliography
12. Appendices

## Chapter outlines
### Preamble
* Abstract
* Preface
	* Data sets: Who made them and which chapters are they used in
	* Analyses: Who performed them
	* Code: Who wrote it
	* Publications: Any arising from chapters of the thesis
* Acknowledgements
	* Supervisors
	* Specific
		* Collaborators
		* WEHI Bioinformatics, especially Melanie Bahlo
		* Felix Krueger, Ryan Lister, Kasper Hansen, Sue Clarke, Rafael Irizarry, Mark Robinson
	* IT, especially Keith and Jakub
	* Other WEHI: Communications, etc.
	* General
		* Family
		* Friends
	* Financial
		* APA, VLSCI, EMBL-PhD, Melanie Bahlo, conference subsidies
* Table of Contents
* List of Figures
* List of Tables
* List of Abbreviations

### Introduction
* Chapter Overview
* DNA and DNA methylation
	* DNA
	* DNA methlyation and its function
	* Types of DNA methylation
* Assays for studying DNA methylation
	* Bisulfite-sequencing
	* Overview of assays (chronological?)
		* Sanger bisulfite-sequencing
			* Can't distinguish mC and hmC
		* Microarrays
		* Enrichment assays
			* methylC-seq, BS-seq, RRBS, PBAT
	* Why use methylC-seq
* Description of methylC-seq
	* History
	* Protocol
	* Data
	* Design of methylC-seq experiments
* Standard analyses
	* Brief, non-mathematical (see also Chapter 4)
* Other biology
	* Primary tisue
	* Cell lines
		* IPS cells
		* ES cells  
* Thesis outline

### Statistical model of methylC-seq data
* Chapter Overview
* Notation
	* See `notation.Rmd`
* Estimating $Z$, $M$, $U$ and $\beta$
	* See `notation.Rmd`
	* Current approaches
* Statistical properties of $z$, $m$, $u$ and $\beta$
	* See `notation.Rmd`
	* Marginal distributions
	* Dependence/correlation over $i$ (see also __Chapter X__)
		* Li _et al._ compute $\beta$ values and compute correlations within a sample as a funciton of genomic distance.
		* Arabidopsis compute $\beta$ values and compute correlations within a sample as a funciton of genomic distance.
		* Lister _et al._ compute both within-read and across-read correlations or dependencies of DNA methylation.
		* Akulenko, R., _et al._ compute gene-level $\beta$ values and compute the correlations between pairs of genes across samples as a function of genomic distance.
	* Dependence/correlation over $j$ (?)
		* Akulenko, R., _et al._ compute gene-level $\beta$ values and compute the correlations between pairs of genes across samples as a function of genomic distance.

### Standard analysis of methylC-seq data
* Chapter Overview
* QC, mapping, post-processing
* Methylation calling
	* Reference-based calling
* Biases
	* M-bias
	* Cell-composition, age, global hypo- or hyper-methylation, SNPs, alignment artefacts
* Univariate tests of $\beta_i$

### m-tuples
* Chapter Overview
* Previous work
	* Amos Tanay
	* Lister _et al._
	* _Arabidopsis_
	* Others
* Within-fragment co-methylation
	* Concept
	* Pros
	* Cons
* Methods
	* `comethylation.py`
	* `cometh`
* Results
	* Lister data
	* Seisenberg data
	* Other?
* Discussion

### Correlation of $\beta$-values
* Chapter Overview
* Previous work
	* _Arabidopsis_
	* _PMBC_
	* Other
	* What are they actually computing?
* Methods
	* `comethylation.py`
	* `cometh`
	* $NIC = 0$ vs. $NIC >= 0$
* Results
	* Lister data
	* Seisenberg data
	* Other
* Discussion

### Simulation model: Part I
* Chapter Overview
* Previous work
	* Simulators for comparing aligners
		* SHERMAN, LAST, etc.
	* Simulators for comparing DMR finders
		* Michelle Lacey, `bsseq`, `BiSeq`, etc.
* Limitations of previous simulation models
* First simple model

### Simulation model: Part II
* Chapter Overview
* A more complex model
* Methods
	* `methsim`
* Results
	* Comparison to real data
	* Limitations
* Discussion

### Data analysis - EPISCOPE
* Chapter overview
* Introduction
* Methods
	* Pre-processing
		* `Bismark`
		* `MarkDuplicates`
	* `comethylation.py`
	* `cometh`
* Results
* Discussion

### Concluding Remarks
### Bibliography
### Appendices

### Others(?)
#### Differential methylation calling
#### Simulation-based inference
* See Brian Ripley's lecture notes (email from Terry)
	* "_The basic idea is quite simple – simulate data from one or more plausible models (or for a parametric model, at a range of plausible parameter values), apply the same (or similar) procedure to the simulated datasets as was applied to the original data, and then analyse the results._"
#### Nick Wong's data

## Proposed chapters in 2-year PhD review
1. Introduction
2. Public datasets used in thesis
3. Pre-processing methylC-seq data and technical biases
4. Statistical model of methylC-seq data
5. Within haplotype dependence structure
6. Across haplotype dependence structure
7. Simulation of methlC-seq data
8. Data analysis
9. Conclusions

## Questions
* Formal literature review chapter or just intersperse lit review at beginning of relevant chapter(s)?
* Terminology (see Comethylation GitHub issue #48)
