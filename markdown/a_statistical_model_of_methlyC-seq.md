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

## One sample
To begin, I consider the simple experiment of performing methylC-seq on a single sample.
Even with only a single sample, this experiment has a hierarchical structure. I separate this description into two stages, as illustrated in __CARTOON FIGURE OF METHYLC-SEQ EXPERIMENT__:

1. Pre-sequencing.
2. Post-sequencing.

Once I have described the framework for a single sample, I extend it to allow for multiple samples. Complications and limitations of this framework are addressed after the general description.

### Pre-sequencing
There are $N$ methylation loci in our sample's genome. A methylation locus is a single cytosine, that is, a CpG, CHG or CHH. Generally, $N$ is not known exactly, although estimates can be made based on a reference genome, but this is no great concern. 

The set of these loci is labelled $\mathcal{I} = \{pos_{i}: i = 1, \ldots, N \}$, where $pos_{i}$ is the genomic co-ordinates of the $i^{th}$ locus with respect to the forward strand, e.g. chr1:723461-723461. I frequently refer to loci by the subscript $i$ rather than by $pos_{i}$. This means that the distance between the $i^{th}$ and $(i + 1)^{th}$ methylation loci varies along the genome and, for a small number of instances, that the $i^{th}$ and $(i + 1)^{th}$ methylation loci are on separate chromosomes.

The methylation state of a locus can vary within a sample due to the fact that DNA for each sample is extracted from hundreds or thousands of cells and each cell can have a slightly different methylation profile. Furthermore, within a diploid cell there are two copies of each chromosome, and therefore two copies of each methylation locus, and these two copies can have different methylation states. Therefore, it is also necessary to consider the next level down in the hierarchy; the DNA fragments within the sample.

I suppose that in the pool of DNA fragments for the sample that there are $H_{i}$ fragments containing the $i^{th}$ methylation locus. In general, $H_{i}$ is unknown and will vary from locus to locus within a sample[^H_i]. __$H_{i}$ is post-PCR; therefore, it can give a grossly distorted picture of the true representation of the cells__

[^H_i]: Knowing $H_{i}$ would require knowing: (1) the number of cells used as input (which is generally only known to within an order of magnitude), (2) the ploidy of each cell (easy) and (3) the number of PCR cycles (easy). But the real problem is that none of the steps in creating the pool of DNA fragments is perfect, in particular, PCR introduces biases -- some molecules are preferrentially amplified while others "drop out". So even if we knew (1), (2) and (3) we cannot simply multiply these together to compute $H_{i}$, although this might at least give us a rough estimate.

Although we do not know the number of fragments in the pool, we can define (and measure) the methylation state of a locus on a single DNA fragment. I denote by the indicator random variable, $Z_{h, i}$, the methylation state of $i^{th}$ methylation locus on the $h^{th}$ DNA fragment containing the $i^{th}$ methylation locus:

\begin{equation}
Z_{h, i, j} = \left\{ 
  \begin{array}{l l}
    1 & \quad \text{if methylated}\\
    0 & \quad \text{if unmethylated}
  \end{array} \right.
\end{equation}

By summing over the number of fragments containing the $i^{th}$ locus, we obtain the number of fragments that are methylated at the $i^{th}$ locus ($M_{i}$) and unmethylated at the $i^{th}$ locus ($U_{i}$):

\begin{align}
	M_{i} &= \sum_{h = 1}^{H = H_{i}} Z_{h, i} \\
	U_{i} &= \sum_{h = 1}^{H = H_{i}} ()1 - Z_{h, i})
\end{align}

Related to these is the proportion of fragments that methylated at the $i^{th}$ locus:
\begin{equation}
	B_{i} = \frac{M_{i}}{M_{i} + U_{i}}
\end{equation}

Again, I emphasise that $H_{i, j}, Z_{h, i, j}, M_{i, j}, U_{i, j} \text{ and } $B_{i, j}$ are unobservable. However, by sequencing the pools of DNA fragments we aim to estimate these variables.

### Post-sequencing
We do not sequence every fragment in the pool. Rather, sequencing can be thought of as sampling without replacement from the pool of DNA fragments. We have a large number(~$10^10$) of fragments in the pool and each methylation locus is only present on a small number of those fragments. Therefore, we can approximate this sampling by Poisson sampling (__CITE__), where the rate parameter for locus $i$ is proportional to the number of fragments in the pool and inversely proportional to $H_{i}$.

At this point I make three simplifying assumptions:

1. We perform single-end sequencing
2. Sequencing is performed without error
3. Read mapping is perfect

I also ignore reads containing no methylation loci as these are not salient to this discussion.

The effect of the first assumption is minor. When using single-end sequencing, the methylation loci from a single read will always form a positively ordered sequence without "gaps", i.e., $(i, i + 1, i + 2)$ and not $(i, i - 1, i - 2)$ nor $(i, i + 1, i + 3)$. However, when using paired-end sequencing the methylation loci from a paired-end read will still be an ordered sequence but one of the following may occur (__FIGURE__):

1. There may be gaps due to the insert size being longer than the sum of the read lengths, e.g. $(i, i + 1, i + 3, i + 4))$. In effect, we have missing data for any intervening methylation loci; the $(i + 2)^{th}$ loci in this example.
2. Loci may be measured twice if the insert size is less than the sum of the read lengths, e.g. read_1 gives us $(i, i + 1)$ and read_2 gives us $(i + 1, i + 2, i + 3)$. In this example we must use only one of read_1 or read_2 as the measurement of the $(i + 1)^{th}$ locus because we are "double-counting" (__SOURCE__)

I discuss the implications of the latter two assumptions in __SECTION__. 

Each read measures the methylation state of one or more loci from a single DNA fragment. I denote by $\mathcal{R}_i$ the set of all reads containing the $i^{th}$ locus. Therefore, the number of reads containing the $i^{th}$ locus is $|\mathcal{R}_{i}|$, where $|\mathcal{R}_{i}| \leq H_{i}$ with strict inequality for almost all $i$.

A single read containing the $i^{th}$ locus is denoted by $z: z \in \mathcal{R}_{i}$; the observed methylation state is given by:

\begin{equation}
z: z \in \mathcal{R}_{i} = \left\{ 
  \begin{array}{l l}
    1 & \quad \text{if methylated at the } i^{th} \text{ locus}\\
    0 & \quad \text{if unmethylated at the } i^{th} \text{ locus}
  \end{array} \right.
\end{equation}

By summing over the number of reads containing the $i^{th}$ locus we obtain the number of reads that are methylated at the $i^{th}$ locus ($m_{i}$) and unmethylated at the $i^{th}$ locus ($u_{i}$):

\begin{equation}
	m_{i} = \sum_{z: z \in \mathcal{R}_{i}} z
	u_{i} = \sum_{z: z \in \mathcal{R}_{i}} (1 - z)
\end{equation}

Similiarly, we obtain the proportion of reads that are methylated at the $i^{th}$ locus as:
\begin{equation}
	\beta_{i} = \frac{m_{i}}{m_{i} + u_{i}}
\end{equation}

These are the so-called $\beta$-values, which are commonly interpreted as an estimate of the proportion of cells in the sample that are methylated at the $i^{th}$ locus. We will discuss this interpretation, and other estimators of the "sample-average" methylation, in __SECTION__.

### Some complications for a single sample
I now discuss some complications and how this framework might accommodate this issues. I also discuss how these issues are dealt with in practice.

#### What is $\mathcal{I}$? {-}

As mentioned in __SECTION__, studies using bisulfite-conversion assays rely on either a reference genome or, less commonly, separate DNA sequencing of the sample that is assayed. Different analysis strategies lead to different definitions of $\mathcal{I}$, which are approximations to the "true" set of methylation loci in the sample, $\mathcal{I}^{truth}$. Listed here are definitions of $\mathcal{I}$ from least closely matching to most closely matching $\mathcal{I}^{truth}$:

1. $\mathcal{I}^{ref}$: Defined by the set of methylation loci in the reference genome. Ignores all genetic variation between the sample and the reference.
2. $\mathcal{I}^{refFilter}$: Defined by filtering out problematic loci from $\mathcal{I}^{ref}$. A conservative approach that removes many sites of genetic variation between the sample and the reference as well as sites that do not display genetic variation between the sample and the reference. Misses sample-specific methylation loci.
3. $\mathcal{I}^{BisSNP}$: Defined by calling genetic variants from the bisulfite-sequencing data using `Bis-SNP` ([@Liu:2012ge]). Identifies sample-specific methylation loci and removes reference-specific methylation loci. The best that can be done if only bisulfite-sequencing data is available. Genetic variation detection is not as good as that from DNA-seq data.
4. $\mathcal{I}^{WGS}$: Defined by identifying all methylation loci from _whole-genome sequencing_ (WGS) of the sample's genome. The gold standard. All methylation loci are defined with respect to the sample's genome. The only differences between $\mathcal{I}^{WGS}$ and the "truth" are due to sequencing errors and variant calling errors of the WGS data.

#### A sample can be heterozygous for a methylation locus  {-}

In a diploid cell there are are sites where, for example, the maternal chromosome is a CpG whereas the paternal chromosome is an ApG. In effect, the maternal and paternal chromosomes within the sample have different $\mathcal{I}^{truth}$. The number of loci that differ between $\mathcal{I}^{maternal}$ and $\mathcal{I}^{paternal}$ is often small enough not to worry about (__I THINK__). However, in some studies, such as those of allele-specific methylation, these loci can be very important. 


* How do I deal with these in practice?
* How might you detect these loci? 
	* Lots of non-C/T (resp. non-G/A) base calls on forward (resp. reverse) strand.
	* WGS of same sample?
	* The trick that BisSNP does?


#### Issues list {-}

__Expand on each of these points:__
* Mapping isn't perfect
* Sequencing error
* More?

## $n$ samples
To move from a single sample for $n$ samples simply requires an additional subscript, $j = 1, \ldots, n$. For example, $\mathcal{I}_{j}$ is the set of methylation loci in the $j^{th}$ sample and $\beta_{i, j}$ is the beta value for the $i^{th}$ locus in the $j^{th}$ sample.

This defines the three levels in the hierarchy of a typical experiment -- individual molecules ($h$), individual methylation loci ($j$) and individual samples ($j$).

#### Some complications for $n$ samples {-}

* Each individual has their own set of methylation loci, i.e. $\mathcal{I}_{j}$ differs for all $j$
* The $i^{th}$ individual has covariates $X_{i}$, which may be a matrix.

[^chimeric_reads]: Might need to note the possibility of chimeric reads.

## m-tuples
* Extend $z: z \in \mathcal{G}_{i}$ to $z: z \in \mathcal{G}_{(i, i + 1)}$. And not that we do not know which $h$ each read came from, only that they came from the same DNA fragment[^chimeric_reads].