# Notes from outline.md
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
		
# A statistical framework for analysing bisulfite sequencing data
## Chapter overview
In this chapter I set out a statistical framework for analysing bisulfite sequencing data. I begin by explaining the various levels of stochasticity in a methylC-seq experiment. Following this, I set out the mathematical notation that I use throughout my thesis and, using this notation, formalise some of the concepts introduced in the previous section. 

I also use this framework to formulate in statistical terminology the key variables and common questions in studies of DNA methylation and compare previous estimators of these key variables.

## A stochastic model of methylC-seq


## Notation
To begin, I consider the simple experiment shown in __FIGURE__. This experiment has $n$ samples, for example, blood samples from $n$ humans. This experiment has a hierarchical structure that I exploit by starting at the level of a sample and working down to the level of a single methylation loci on a single fragment of DNA within that sample.

I separate this description into two stages, as illustrated in __CARTOON FIGURE OF METHYLC-SEQ EXPERIMENT__:

1. Prior to sequencing. At this point we have a pool containing a large (~$10^10$)) number of bisulfite-treated short (~100-500nt) DNA fragments.
2. Following sequencing.

This description momentarily ignore some complications but I address these in __SECTION__.


### Prior to sequencing
The $j^{th}$ sample's genome contains $N_{j}$ methylation loci. A methylation locus is a single cytosine, that is, a CpG, CHG or CHH. However, the methylation state of a locus can vary within a sample due to the fact that DNA for each sample is extracted from hundreds or thousands of cells and each cell can have a slightly different methylation profile. Furthermore, within a diploid cell there are two copies of each chromosome, and therefore two copies of each methylation locus, and these two copies can have different methylation states. Therefore, it is also necessary to consider the next level down in the hierarchy; the DNA fragments within the sample.

I suppose that in the pool of DNA fragments for sample $j$ that there are $H_{i, j}$ fragments containing the $i^{th}$ methylation locus. In general, $H_{i, j}$ is unknown and will vary from locus to locus within a sample[^H_ij]. 

[^H_ij]: Knowing $H_{i, j}$ would require knowing: (1) the number of cells used as input (which is generally only known to within an order of magnitude), (2) the ploidy of each cell (easy) and (3) the number of PCR cycles (easy). But the real problem is that none of the steps in creating the pool of DNA fragments is perfect, in particular, PCR introduces biases -- some molecules are preferrentially amplified while others "drop out". So even if we knew (1), (2) and (3) we cannot simply multiply these together to compute $H_{i, j}$, although this might at least give us a rough estimate.

Although we do not know the number of fragments in the pool, we can define (and measure) the methylation state of a locus on a single DNA fragment. I denote by $Z_{h, i, j}$ the methylation state on the $h^{th}$ DNA fragment containing the $i^{th}$ methylation locus on the for the $j^{th}$ sample.

\begin{equation}
Z_{h, i, j} = \left\{ 
  \begin{array}{l l}
    1 & \quad \text{if methylated}\\
    0 & \quad \text{if unmethylated}
  \end{array} \right.
\end{equation}

__TODO__: Define $M, U, B$, etc.

This defines the three levels of the hierarchy -- individual molecules ($h$), individual methylation loci ($j$) and individual samples ($j$).

### Following sequencing
 
* Sequencing as sampling
* The $i^{th}$ individual has covariates $X_{i}$, which may be a matrix.
* Sequencing error as measurement error
* Aggregate measurements $Y_{i, j}$

### Some complications
 
 __Things I've ignored__

 * Each individual has their own set of methylation loci
 * A sample can be heterozygous for a CpG 



