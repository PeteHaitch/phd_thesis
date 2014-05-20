# A statistical framework for analysing bisulfite sequencing data

## Notes from outline.md (plus other notes added along the way)
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
		* Li _et al._ compute $\beta$ values and compute correlations within a sample as a function of genomic distance.
		* Arabidopsis compute $\beta$ values and compute correlations within a sample as a function of genomic distance.
		* Lister (2009) Fig. 3 and Sup. Fig. 9
		* Lister _et al._ compute both within-read and across-read correlations or dependencies of DNA methylation.
		* Akulenko, R., _et al._ compute gene-level $\beta$ values and compute the correlations between pairs of genes across samples as a function of genomic distance.
	* Dependence/correlation over $j$ (?)
		* Akulenko, R., _et al._ compute gene-level $\beta$ values and compute the correlations between pairs of genes across samples as a function of genomic distance.

## Chapter overview
In this chapter I set out a statistical framework for analysing bisulfite sequencing data. I begin by explaining the various levels of stochasticity in a methylC-seq experiment. Following this, I define the mathematical notation that I use throughout my thesis and formalise some of the concepts introduced in the previous section. 

I also use this framework to formulate common questions in studies of DNA methylation. In particular, I define key variables and describe the statistical properties of common estimators of these variables.

## One sample

To begin, I consider the simple experiment of performing methylC-seq on a single sample.
Even with only a single sample, this experiment has a hierarchical structure. I separate this structure into two stages, as illustrated in __CARTOON FIGURE OF METHYLC-SEQ EXPERIMENT__:

1. Pre-sequencing
2. Post-sequencing

Once I have described the framework for a single sample, I extend it to allow for multiple samples. The ideas here are simple, although when extended to their full generality the notation becomes messy. I address complications and limitations of this framework at the end of this sub-section. 

### Single locus analysis

As described in __SECTION__, most analyses to date of WGBS data have focused on the "sample-average" level methylation at individual loci. For this reason, I initially describe the framework in terms of single loci. However, bisulfite sequencing data contain much more information than provided by these univariate, marginal summaries of the data. In __SECTION__ I introduce _m-tuples_ that can summarise some of this extra information and which I use extensively in __CHAPTERS__.


#### Pre-sequencing {-}
A methylation locus is a single cytosine, that is, a CpG, CHG or CHH. The set of these loci is labelled $\mathcal{I} = \{pos_{i}: i = 1, \ldots, N_{loci} \}$, where $pos_{i}$ is the genomic co-ordinates of the $i^{th}$ locus with respect to the forward strand, e.g. chr1:723,461-723,461. I frequently refer to loci by the subscript $i$ rather than by $pos_{i}$. This means that the distance between the $i^{th}$ and $(i + 1)^{th}$ methylation loci varies along the genome and, for a small number of instances, that the $i^{th}$ and $(i + 1)^{th}$ methylation loci are on separate chromosomes. Generally, $N_{loci}$ is not known exactly, although estimates can be made based on a reference genome, but this is no great concern. 

The methylation state of a locus can vary within a sample due to the fact that DNA for each sample is extracted from hundreds or thousands of cells and each cell can have a slightly different methylation profile. Furthermore, within a diploid cell there are two copies of each chromosome, and therefore two copies of each methylation locus, and these two copies can have different methylation states. Therefore, it is also necessary to consider the next level down in the hierarchy; the DNA fragments within the sample.

I suppose that in the pool of DNA fragments for the sample that there are $H_{i}$ fragments containing the $i^{th}$ methylation locus. In general, $H_{i}$ is unknown and will vary from locus to locus within a sample[^H_i]. __$H_{i}$ is post-PCR; therefore, it can give a grossly distorted picture of the true representation of the cells__. I denote by $\mathcal{H}_{i}$ the set of all fragments containing the $i^{th}$ locus.

[^H_i]: Knowing $H_{i}$ would require knowing: (1) the number of cells used as input (which is generally only known to within an order of magnitude), (2) the ploidy of each cell (easy) and (3) the number of PCR cycles (easy). But the real problem is that none of the steps in creating the pool of DNA fragments is perfect, in particular, PCR introduces biases -- some molecules are preferrentially amplified while others "drop out". So even if we knew (1), (2) and (3) we cannot simply multiply these together to compute $H_{i}$, although this might at least give us a rough estimate.

Although we do not know the number of fragments in the pool, we can define (and measure) the methylation state of a locus on a single DNA fragment. I denote by the indicator random variable, $Z_{h, i}$, the methylation state of $i^{th}$ methylation locus on the $h^{th}$ DNA fragment containing the $i^{th}$ methylation locus:

\begin{equation*}
Z_{h, i} = \left\{ 
  \begin{array}{l l}
    1 & \quad \text{if methylated on the } h^{th} \text{ fragment}\\
    0 & \quad \text{if unmethylated on the } h^{th} \text{ fragment}
  \end{array} \right.
\end{equation*}

By summing over the number of fragments containing the $i^{th}$ locus, we obtain the number of fragments that are methylated at the $i^{th}$ locus ($M_{i}$) and unmethylated at the $i^{th}$ locus ($U_{i}$):

\begin{align*}
	M_{i} &= \sum_{h = 1}^{H = H_{i}} Z_{h, i} = |\{Z: Z \in \mathcal{H}_{i}, Z = 1 \}| \\
	U_{i} &= \sum_{h = 1}^{H = H_{i}} (1 - Z_{h, i}) = |\{Z: Z \in \mathcal{H}_{i}, Z = 0 \}|
\end{align*}

Related to these is the proportion of fragments that methylated at the $i^{th}$ locus:
\begin{equation*}
	B_{i} = \frac{M_{i}}{M_{i} + U_{i}}
\end{equation*}

Again, I emphasise that $H_{i, j}, Z_{h, i, j}, M_{i, j}, U_{i, j} \text{ and } B_{i, j}$ are unobservable. However, by sequencing the pools of DNA fragments we aim to estimate these variables.

#### Post-sequencing {-}

We do not sequence every fragment in the pool. Rather, sequencing can be thought of as sampling without replacement from the pool of DNA fragments. We have a large number(~$10^{10}$) of fragments in the pool and each methylation locus is only present on a small number of those fragments. Therefore, we can approximate this sampling by Poisson sampling (__CITE__), where the rate parameter for locus $i$ is proportional to the number of fragments in the pool and inversely proportional to $H_{i}$.

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

\begin{equation*}
z: z \in \mathcal{R}_{i} = \left\{ 
  \begin{array}{l l}
    1 & \quad \text{if methylated at the } i^{th} \text{ locus}\\
    0 & \quad \text{if unmethylated at the } i^{th} \text{ locus}
  \end{array} \right.
\end{equation*}

By summing over the number of reads containing the $i^{th}$ locus we obtain the number of reads that are methylated at the $i^{th}$ locus ($m_{i}$) and unmethylated at the $i^{th}$ locus ($u_{i}$):

\begin{align*}
	m_{i} &= \sum_{z: z \in \mathcal{R}_{i}} z \\
	      &= |\{z: z \in \mathcal{R}_{i}, z = 1 \}| \\
	u_{i} &= \sum_{z: z \in \mathcal{R}_{i}} (1 - z) \\
		  &= |\{z: z \in \mathcal{R}_{i}, z = 0 \}
\end{align*}

Similiarly, we obtain the proportion of reads that are methylated at the $i^{th}$ locus as:
\begin{equation*}
	\beta_{i} = \frac{m_{i}}{m_{i} + u_{i}}
\end{equation*}

These are the so-called $\beta$-values, which are commonly interpreted as an estimate of the proportion of cells in the sample that are methylated at the $i^{th}$ locus. We will discuss this interpretation, and other estimators of the "sample-average" methylation, in __SECTION__. I now move from single locus summaries to multi-loci summaries using m-tuples.

### m-tuples: Multiple loci analysis

I define an m-tuple to be a tuple of methylation loci, where m = $1, 2, \ldots$ is the size of the tuple. For example, a CpG 3-tuple is 3 CpGs. A CHH 1-tuple is a single CHH; 1-tuples are the basis of single locus analyses of bisulfite sequencing data. We can learn more from bilsufite sequencing data by using m-tuples with m $> 1$.

An observation on an m-tuple is the methylation pattern from a single read that overlaps the m-tuple. I require that the observation is from a single read as this ensures that each observation ultimately comes from a single cell[^chimeric_reads]. For example, suppose we have a read containing 3 CpGs -- the first two CpGs are methylated while the last one is unmethylated. From this read we can observe 3 $\times$ CpG 1-tuples, 2 $\times$ CpG 2-tuples and 1 $\times$ CpG 3-tuple. We can have multiple observations on an m-tuple by observing multiple reads, each containing the m-tuple. __FIGURE__ illustrates this example.

[^chimeric_reads]: Might need to note the possibility of chimeric reads. 

Note that in __FIGURE__ I haven't constructed a 2-tuple using the first CpG and the last CpG. In general, I focus on m-tuples where the m methylation loci are adjacent. That is, I focus on m-tuples where the _number of intervening loci_ is zero ($NIL = 0$). There are 3 reasons for this:

1. Quantity: From a sequence containing $M$ methylation loci there are $M - \text{m} + 1$ m-tuples, provided that we restrict ourselves to those m-tuples with $NIL = 0$. In contrast, if we allow $NIL \geq 0$ then there are $\binom{M}{\text{m}} ~ $ m-tuples. (__TODO__: Describe how these terms grow asymptotically).
2. Interpretability: Discussed in __SECTION__
3. Measurability: We cannot observe m-tuples where the methylation loci are far apart due to the read length limitations of the Illumina sequencing technology. This is true even when $NIL = 0$ but is especially the case if we allow $NIL \geq 0$.

When I refer to m-tuples I implicitly mean $NIL = 0$; I will explicitly use the notation $NIL \geq 0$ when I wish to make clear that there may be intervening methylation loci in the m-tuple.

For each m-tuple, I define the _intra-pair distance_ (IPD) as vector containing the $m - 1$ pair-wise distances (measured in bp) between adjacent methylation loci in the m-tuple. For example, the 2-tuple (chr7:145, 163) has $IPD = (163 - 145) = (18)$. The 5-tuple (chr2:560, 570, 572, 588, 612) has $IPD = (570 - 560, 572 - 570, 588 - 572, 612 - 588) = (10, 2, 16, 24)$. The IPD vector of a 1-tuple is not defined. In __CHAPTER__ I use the IPDs, along with other features such as the genomic context, to define "similar" m-tuples.

#### Pre-sequencing {-}

Mathemtically, an m-tuple is denoted by a sequence of methylation loci, $(i, i + 1, \ldots, i + m - 1)$. The remaining definitions are analogous to those for single methylation loci, that is, when m $= 1$.

I denote by the vector of indicator random variable, $Z_{h, (i, i + 1, \ldots, i + \text{m} - 1)}$, the methylation pattern of the m-tuple $(i, i + 1, \ldots, i + m - 1)$ on the $h^{th}$ DNA fragment containing the m-tuple $(i, i + 1, \ldots, i + \text{m} - 1)$:

\begin{equation*}
Z_{h, i} = \left\{ 
  \begin{array}{l l}
    (0, 0, \ldots, 0) & \quad \text{if unmethylated at the } i^{th}, (i + 1)^{th}, \ldots, (i + \text{m} - 1)^{th}) \text{ locus on the } h^{th} \text{ fragment}\\
    (0, 0, \ldots, 1) & \quad \text{if unmethylated at the } i^{th}, (i + 1)^{th}, \ldots, (i + \text{m} - 2)^{th}) \text{ locus and methlyated at the } (i + \text{m} - 1)^{th} \text{ locus on the } h^{th} \text{ fragment} \\
    \vdots \\
    (1, 1, \ldots, 1) & \quad \text{if methylated at the } i^{th}, (i + 1)^{th}, \ldots, (i + \text{m} - 1)^{th}) \text{ locus on the } h^{th} \text{ fragment}
   \end{array} \right.
\end{equation*}

$\mathcal{H}_{(i, i + 1, i + m - 1)}$ denotes the set of all fragments containing the m-tuple $(i, i + 1, i + m - 1)$ and $\mathcal{R}_{(i, i + 1, i + m - 1)}$ denotes the set of all reads containing the m-tuple $(i, i + 1, i + m - 1)$.

There are $2^{\text{m}}$ possible methylation patterns at an m-tuple. I also write these using $U$ and $M$ instead of $0$ and $1$; for example, the possible methylation patterns at a 2-tuple are $MM, MU, UM$ and $UU$. 

Analogously to the definition of $M_{i}$ and $U_{i}$ and $m_{i}$ and $u_{i}$ for the case where $\text{m} = 1$, we have when $\text{m} = 2$:

\begin{align*}
	MM_{(i, i + 1)} &= |\{Z: Z \in \mathcal{H}_{(i, i + 1)}, Z = (1, 1)\}| \\
	MU_{(i, i + 1)} &= |\{Z: Z \in \mathcal{H}_{(i, i + 1)}, Z = (1, 0)\}| \\
	UM_{(i, i + 1)} &= |\{Z: Z \in \mathcal{H}_{(i, i + 1)}, Z = (0, 1)\}| \\
	UU_{(i, i + 1)} &= |\{Z: Z \in \mathcal{H}_{(i, i + 1)}, Z = (0, 0)\}| \\
\end{align*}

We can extend the $B_{i}$ values to m-tuples, although the intuitive interpretation is somewhat lost. Here are the definitions for $\text{m} = 2$:
\begin{align*}
	B_{(i, i + 1)}^{MM} &= \frac{MM_{(i, i + 1)}}{MM_{(i, i + 1)} + MU_{(i, i + 1)} + UM_{(i, i + 1)} + UU_{(i, i + 1)}} \\
	B_{(i, i + 1)}^{MU} &= \frac{MU_{(i, i + 1)}}{MM_{(i, i + 1)} + MU_{(i, i + 1)} + UM_{(i, i + 1)} + UU_{(i, i + 1)}}  \\
	B_{(i, i + 1)}^{UM} &= \frac{UM_{(i, i + 1)}}{MM_{(i, i + 1)} + MU_{(i, i + 1)} + UM_{(i, i + 1)} + UU_{(i, i + 1)}}  \\
	B_{(i, i + 1)}^{UU} &= \frac{UU_{(i, i + 1)}}{MM_{(i, i + 1)} + MU_{(i, i + 1)} + UM_{(i, i + 1)} + UU_{(i, i + 1)}} 
\end{align*}
The definitions for $\text{m} > 3$ follow in the obvious manner.

Again, I emphasise that $H_{(i, i + 1, i + \text{m} - 1), j}, Z_{h, (i, i + 1, i + \text{m} - 1), j}, B_{(i, i + 1, i + \text{m} - 1), j}$ and the set of methylation patterns are unobservable. However, by sequencing the pools of DNA fragments we aim to estimate these variables.

#### Post-sequencing {-}

The the set of all reads containing the m-tuple $(i, i + 1, \ldots, i + \text{m} - 1)$ is denoted by $\mathcal{R}_i$. A single read containing the m-tuple $(i, i + 1, \ldots, i + \text{m} - 1)$ is denoted by $z: z \in \mathcal{R}_{(i, i + 1, \ldots, i + \text{m} - 1)}$; the observed methylation state is given by:

\begin{equation*}
z: z \in \mathcal{R}_{(i, i + 1, \ldots, i + \text{m} - 1)} = \left\{ 
  \begin{array}{l l}
    (0, 0, \ldots, 0) & \quad \text{if unmethylated at the } i^{th}, (i + 1)^{th}, \ldots, (i + \text{m} - 1)^{th}) \text{ locus}\\
    (0, 0, \ldots, 1) & \quad \text{if unmethylated at the } i^{th}, (i + 1)^{th}, \ldots, (i + \text{m} - 2)^{th}) \text{ locus and methlyated at the } (i + \text{m} - 1)^{th} \text{ locus} \\
    \vdots \\
    (1, 1, \ldots, 1) & \quad \text{if methylated at the } i^{th}, (i + 1)^{th}, \ldots, (i + \text{m} - 1)^{th}) \text{ locus}
      \end{array} \right.
\end{equation*}


Note that we do not know which $h$ each read came from, only that all methylation loci in the read came from the same DNA fragment.

By "summing" over the number of reads containing the m-tuple $(i, i + 1, \ldots, i + \text{m} - 1)$ we obtain the number of reads contain each methylation pattern. Here are the definitions for $\text{m} = 2$:

\begin{align*}
	mm_{(i, i + 1)} &= |\{z: z \in \mathcal{R}_{(i, i + 1)}, z = (1, 1)\}| \\
	mu_{(i, i + 1)} &= |\{z: z \in \mathcal{R}_{(i, i + 1)}, z = (1, 0)\}| \\ 
	um_{(i, i + 1)} &= |\{z: z \in \mathcal{R}_{(i, i + 1)}, z = (0, 1)\}| \\ 
	uu_{(i, i + 1)} &= |\{z: z \in \mathcal{R}_{(i, i + 1)}, z = (0, 0)\}| 
\end{align*}

The definitions for $\text{m} > 2$ follow in the obvious manner.

We can extend the $\beta_{i}$ values to m-tuples, although the intuitive interpretation is lost when m $> 1$. Here are the definitions for $\text{m} = 2$:
\begin{align*}
	\beta_{(i, i + 1)}^{mm} &= \frac{mm_{(i, i + 1)}}{mm_{(i, i + 1)} + mu_{(i, i + 1)} + um_{(i, i + 1)} + uu_{(i, i + 1)}} \\
	\beta_{(i, i + 1)}^{mu} &= \frac{MU_{(i, i + 1)}}{mm_{(i, i + 1)} + mu_{(i, i + 1)} + um_{(i, i + 1)} + uu_{(i, i + 1)}} \\
	\beta_{(i, i + 1)}^{um} &= \frac{UM_{(i, i + 1)}}{mm_{(i, i + 1)} + mu_{(i, i + 1)} + um_{(i, i + 1)} + uu_{(i, i + 1)}} \\
	\beta_{(i, i + 1)}^{uu} &= \frac{UU_{(i, i + 1)}}{mm_{(i, i + 1)} + mu_{(i, i + 1)} + um_{(i, i + 1)} + uu_{(i, i + 1)}}
\end{align*}
The definitions for $\text{m} > 3$ follow in the obvious manner.

__TODO__: Change definition of $M$, $m$, etc. when $m = 1$

\citet{Landan:2012kp} compute the average methylation of a read containing the m-tuple $(i, \ldots, i + m - 1)$. For each read, $z \in \mathcal{R}_{(i, i + 1, i + 2, i + 3)}$, the average methylation of the read, $\zeta_{z}$, is defined as the proportion of methylation loci in the read that are methylated. Thus, $\zeta_{z} = 0, \frac{1}{m}, \frac{2}{m}, \ldots, 1$.



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

Parent-specific methylation loci can be identified by calling heterozygous genetic variants using `Bis-SNP` or from WGS of the sample. In practice, the existance of parent-specific methylation loci is often ignored.

#### Sequencing error destroys methylation signal {-}

Sequencing errors from the Illumina technology are typically base substitutions, for example, a G being read as a T. These errors occur more frequently at the 3' end of the reads (__SOURCE__). Other sequencing technologies have different error profiles (__SOURCE__). In DNA sequencing, the base quality scores are used to weed out or down-weight likely sequencing errors. However, the base quality scores of methylC-seq data are not as reliable as those generated for DNA-seq because the base quality scoring algorithm is tuned for DNA sequencing (__SOURCE or EXAMPLE__).

There are two negative effects of sequencing errors. Firstly, if a methylated cytosine on the forward strand, which should be a C in the read, is erroneously read as another base -- A, G or T -- then the methylation signal is lost. If the erroneous basecall is an A or a G then at least we might be able to infer that a sequencing error occured. However, if the erroneous basecall is a T then we might instead mistakenly conclude that the cytosine was unmethylated.


The second effect of sequencing errors is on read mapping, which I discuss below.

#### Mapping against a reference genome produces errors {-}

Read mapping is not perfect and produces both false positive and false negative results. False positives are reads mapped to the wrong location and reads mapped to multiple locations with equal mapping scores. False negatives are reads that are not mapped to any location; these reads are effectively lost from any downstream analysis. The parameter settings used by the mapping software determine the false positive and false negative rates. 

There are biological and technical reasons why mapping against a reference genome can produce these errors. Biologically, if the sample contains sequences that are too genetically divergent from the reference genome then these sequences will be difficult, even impossible, to map. A particularly problematic class of sequences are those from repetitive regions of the genome. These repetitive sequences will map to multiple locations in the reference genome equally well. Furthermore, the number of times these repetitive sequences occur differs between the reference genome and sample's genome. 

Technically, reads from Illumina sequencing are often too short to "resolve" the mapping location of these repetitive sequences. Resolving the mapping location of repetitive sequences can be achieved byusing other sequencing technologies, such as Pacific Biosciences (__CITE__) and Oxford Nanopore (__CITE__), which produce longer reads.

Another source of technical error in read mapping is sequencing error. A sequencing error can transform a uniquely mapping read to one that maps equally well to multiple locations or, worse still, a read that maps uniquely to a new but incorrect location. Sequencing errors can also corrupt a read so badly that it no longer can be mapped, which leads to that read being a false negative.

In practice, most people try to mitigate these problems through their choice of parameters used by the mapping software. Many studies have been published that seek to provide the "optimal" parameters for a variety of common scenarios (__CITE__). 

In theory, reads might be down-weighted in downstream analyses based on the mapping quality score. Ideally, mapping software assigns the degree of confidence it has that the read is "correctly" mapped via a mapping quality score (`mapQ`). However, these mapping quality scores are often poorly calibrated, particularly for methylC-seq data, which makes them less useful. For this reason, the popular methylC-seq mapping software, Bismark (__VERSION__), does not even provide mapping quality scores.

All these problems are general challenges of read mapping and not specific to methylC-seq data, although the reduced complexity of methylC-seq reads exacerbates these issues. The difficulty of mapping to repetitive regions of the genome is a particularly frustrating one for methylC-seq data. Repetitive sequences, such as LINEs and SINEs (__CHECK__), are frequently methylated and so of interest to methylation researchers. Poor mapping efficiencies of these regions means that there is often limited data of these elements from methylC-seq data (__CHECK IN MY DATA__).

__Move some of this section to Introduction__

__Other issues?__


## $n$ samples

__DISCUSS WITH TERRY: Is it better to use $j$ to denote samples and encode group information in $X_j$, or to use $j$ to denote replicates and $k$ to denote groups? I think the former will be more general, e.g. a sample belonging to multiple groups.__

To move from a single sample to $n$ samples simply requires an additional subscript, $j = 1, \ldots, n$. $\mathcal{I}_{j}$ is the set of methylation loci in the $j^{th}$ sample and $\beta_{i, j}$ is the beta value for the $i^{th}$ locus in the $j^{th}$ sample. This defines the three levels in the hierarchy of a typical experiment -- individual molecules ($h$), individual methylation loci ($j$) and individual samples ($j$). A fourth level is how the samples relate in terms of an outcome of interest, such as phenotype. This fourth level might be defined up-front, such as in a designed experiment looking for differences in methylation between pre-defined groups of samples, or the aim of the experiment might be to discover this level.

A common experiment of the first type is the two-group design in which $n_{1}$ samples are from group $1$ and $n_{2}$ samples are from group $2$ ($n_{1} + n_{2} = n$). This can be represented by a design matrix $X = [X_{j}]$, where $X_{j} = 1$ if the sample is from group $1$ and $X_{j} = 0$ if the sample is from group $2$. 

More generally, we can study multi-group experiments via a suitable definition of the design matrix $X$. We can also include covariates by allowing $X_{j}$ to be a row vector, $X_{j} = (x_{1, j}, \ldots, x_{P, j})$, where $x_{p, j}$ encodes the information on the $p^{th}$ covariate for the $j^{th}$ sample.

The second type of experiment is one that aims to cluster individuals based on methylation patterns. This can be thought of as an experiment where the outcome of interest is unknown and we wish to discover it by lumping together samples with "similar" methylation patterns. We might also wish to then correlate these clusters with some external information such as disease-status. 

### Some complications for $n$ samples

In addition to the complications of the previous section, we now have sample-to-sample variation.

#### $\mathcal{I}$ {-}

Each individual has their own set of methylation loci, that is, $\mathcal{I}_{j}$ differs for all $j$. Furthermore, sequencing coverage varies from sample-to-sample. This means that even if the samples have exactly the same $\mathcal{I}_{j}$, e.g. the samples are genetically identical, each sample will have a different set of loci with "sufficient" sequencing coverage. Loci without sufficient coverage are effectively missing data.

In practice, we might choose to study $\mathcal{I}^{common} = \cap_{j} \mathcal{I}_{j}$ or some other combination of the $\mathcal{I}_{j}$, such as all methylation loci present in at least some fraction of the $n$ samples. 

A conservative analysis might only analyse those loci where at least some fraction of the $n$ samples have sufficient sequencing coverage. A less conservative analysis might try to impute the missing values based on methylation levels at neighbouring loci (__bsseq__).

#### Complications that aren't merely notational {-}

* $n$ is typically small while $N_{loci}$ is typically large
* Sequencing coverage may be low. Suggests smoothing.


## Parameter estimation
In this section I summarise current techniques for parameter estimation from WGBS data. I do not describe the processing of the raw data. When necessary, I have 'translated' the original work into my notation to make these methods more readily comparable.

All methods estimate parameters from files containing the aligned reads. 

### Estimating $M$ and $U$

The simplest and most commonly used estimators of $M_{i}$ and $U_{i}$ are $m_{i}$ and $u_{i}$, that is, the number of reads that are methylated and unmethylated at the $i^{th}$ locus (e.g. \citet{Cokus:2008fc, Lister:2008bh, Lister:2009hy, Lister:2011kg, Hansen:2011gu}). $M$ and $U$ may be estimated per-strand or aggregated across strands[^strand_aggregation]. CpG methylation is typically aggregated across strands.

[^strand_aggregation]: Estimates can only be aggregated across strands if the methylation marks is symmetric across strands, which CpG is and CHG and CHH methylation are not __TODO: CHECK SYMMETRY ARGUMENT__.

Not all reads, and positions within reads, are equally treated in this counting process. Essentially, each sequenced nucleotide is assigned a weight. By far the most used weighting scheme is a set of filters that assigns each sequenced nucleotide a weight of 0 (observation excluded) or 1 (observation included). 

When using a set of filters, at each step a read either "survives", and is subjected to the proceeding filter, or "dies", and is excluded from the estimation of $M$ and $U$[^read_filters]. While strictly speaking each sequenced nucleotide is assigned a weight, in practice weights are normally first assigned to reads and then to all nucleotides within "surviving" reads. 

[^read_filters]: A read that is not used to estimate $M$ and $U$ may still be used in other analyses, such as estimating copy number variation, and vice-versa.


__TABLE___ describes a standard set of filters for bisulfite-sequencing data, although the choice of filters and thresholds strongly depends on the assay used, design of the experiment and data quality. In the table I have also included the relevant argument for my `comethylation` software, which can estimate $M$ and $U$ using the simple counting procedure.

#### Read-level filters

A read survives if:

1. Read is mapped (SE, PE) and mapped in the expected orientation (PE only).
2. Read is not marked as a PCR duplicate. `--ignoreDuplicates`.
3. Read has a mapping quality score (`mapQ`) greater than some threshold. `minMapQ <int>`.
4. Read does not contain more than a certain level of non-CpG methylation. 
5. Read passes aligner-specific filters designed to remove known biases of the alignment software.

#### Base-level filters {-}

A sequenced base survives if:

1. Position of base in read means that it is unlikely to be effected by significant M-bias.
2. Base has Phred quality score greater than some threshold.
3. Base is a "methylation mismatch" (e.g. the sequenced base is a C or T at a C in the reference sequence) and not a "non-methylation mismatch" (e.g. the sequenced base is an A or a G at a C in the reference sequence).


#### Additional filters {-}

Paired-end sequencing data often contain reads where `read_1` and `read_2` overlap. It is important to ensure that any positions in the overlap are only counted once, otherwise evidence

#### More sophisticated weights {-}

An obvious extension is to weight a methylation call by its base quality or `mapQ` value of the read. However, as previously noted, these qualities are often not well calibrated for bisulfite sequencing data, which reduces their utility. Similarly, it might be possible to weight a methylation call by the expected level of M-bias. However, M-bias, as it is currently computed, also has problems with calibration (__SEE INTRODUCTION__).

### Estimating methylation patterns at m-tuples

`comethylation` is the only software that estimates methylation patterns at m-tuples, such as $MM, MU, UM$ and $UU$ for 2-tuples. `comethylation` uses the same "filter + count" method to estimate methylation patterns for all m-tuples (m $= 1, 2, 3, \ldots $).


### Estimating $B$

Most analyses of bisulfite sequencing data have focused on the on the average level of methylation at individual methylation loci, $B_{i}$, and use the simple estimator[^simple_beta] $\beta_{i} = \frac{m_{i}}{m_{i} + u_{i}}$ (e.g. \citet{Cokus:2008fc, Lister:2008bh, Lister:2009hy, Lister:2011kg}). Recently, more sophisticated estimators have also been proposed, such as the empirical Bayes frameworks separately proposed by \citet{Feng:2014iq} and \citet{Sun:2014fk} and the smoothing-methods proposed by \cite{Hansen:2011gu, Hansen:2012gr} and \citet{Hebestreit:2013ko}.

[^simple_beta]: Under the Binomial model, $M_{i} | (M_{i} + U_{i}) = Bin(M_{i} + U_{i}, B_{i})$, $\beta_{i} = \frac{m_{i}}{m_{i} + u_{i}}$ corresponds to the maximum likelihood estimator (__CHECK__). It would also have a Bayesian interpretation, I think (__CHECK__).

#### Empirical Bayes models of $\beta$-values {-}

Both \citet{Feng:2014iq} and \citet{Sun:2014fk} propose a Beta-Binomial empirical Bayes hierarchical model for the $B_{i}, i = 1, \ldots, N_{loci}. The actual models slightly differ, as do the algorithms for parameter estimation, but the main idea is the same. The R/Bioconductor package `DSS` implements the model of \citet{Feng:2014iq} and the `MOABS` software implements the model of \citet{Sun:2014fk}. I focus on \citet{Feng:2014iq} because the model is better described than that of \citet{Sun:2014fk}. In fact, the empirical Bayes method described in \citet{Sun:2014fk} is poorly written and I think is wrong in several places, which makes it very confusing. (__DISCUSS WITH TERRY___). 
  
\citet{Feng:2014iq} model the number of methylated reads at each locus by $M_{i, j, k} = Binomial(m_{i, j, k} + u{i, j, k}, B_{i, j, k})$, where $k = 1, 2$ is the group of each sample in a two-group experiment. They then assume that the $B_{i, j, k}$ follow a $Beta(\mu{i, k}, \theta_{i, k})$ distribution, which is the conjugate prior for the Binomial distribution, where $\mu$ is the mean and $\theta$ is the dispersion. \citet{Feng:2014iq} make the additional modelling assumption that $\theta_{i, k} = \text{log-Normal}(m_{0, k}, r_{0, k}^{2})$. Posterior estimates of $\mu_{i, k}$ are obtained using an empirical Bayes framework. 

__TODO: DISCUSS WITH TERRY - Feng perform tests based on  the parameters from the _prior_ distribution, not the posterior. That doesn't make sense, does it? Sun estimate the posterior parameters $\beta_{i, j, k}$.__

#### Smoothing $\beta$-values {-}

`BSmooth`, published in \citet{Hansen:2011gu, Hansen:2012gr} and available in the R/Bioconductor package `bsseq`, and `BiSeq`, published in \citet{Hebestreit:2013ko} and available in the R/Bioconductor package `BiSeq`, take a different approach to getting improved estimates of the $\beta$-values. Rather than developing an empirical Bayes model, `BSmooth` and `BiSeq` both use statistical smoothing of the "raw" $\beta_{i} = \frac{m_{i}}{m_{i} + u_{i}}$. Smoothing is motivated and justified by the fact that the $\beta$-values are spatially correlated. 

Smoothing is particularly powerful for loci with low sequencing coverage, where the denominimator $m_{i} + u_{i}$ is small and the corresponding standard error estimates of $\beta_{i}$ are large. Is smoothing of the raw $\beta$-values still useful when you have high-coverage sequencing data (__DISCUSS WITH TERRY__)? The smoothed $\beta$-values, and not the raw $\beta$-values, are then generally used in all downstream analyses.

Both `BSmooth` and `BiSeq` use a binomial local likelihood smoother. This smoother was chosen because `BSmooth` and `BiSeq` model the number of methylated reads at the $i^{th}$ locus in the $j^{th}$ sample by $M_{i, j} = Binom(m_{i, j} + u_{i, j}, B_{i})$ and it is "local" because "methylation levels are strongly correlated across the genome" \citep{Hansen:2012gr}. \citet{Hansen:2012gr} cite \citet{Eckhardt:2006gh} as evidence that DNA methylation levels are similar at proximal CpGs). 

The raw $\beta$-values are weighted according to the binomial likelihood and the kernel function. The binomial likelihood weights points inversely to their standard error, $se(\beta_{i, j})$, and the kernel gives greater weight to those $\beta_{i, j}$ near the centre of the window. \citet{Lacey:2013iy} note that loci with very high sequencing coverage will strongly influence the smoother, potentially biasing estimates at neighbouring loci with low coverage.

`BSmooth` assumes that for each sample that the underlying methylation level, $f_{j}(i)$, is a smoothly varying function of the position in the genome, $i$. In contrast, `BiSeq` first creates clusters of CpGs and only assumes that the underlying methylation level is smooth at positions within each cluster.

When smoothing, a key parameter is the bandwidth, which is the size of the window in which observations are included at each iteration of the smoother. `BSmooth` uses a much larger window size than `BiSeq`; the default window size in `BSmooth` is one that contains at least 70 CpGs and is at least 2000kb wide, whereas the default window size in `BiSeq` is 80bp, regardless of CpG-density. This is due to `BiSeq` being developed for RRBS data, which has a high CpG-density per window, whereas `BSmooth` was developed for WGBS data, which has a lower and more variable CpG density per window.

Another 'parameter' choice when smoothing is the choice of kernel, although this is less important than the choice of bandwidth, __correct?__. `BSmooth` uses a tricube kernel and `BiSeq` uses a triangular kernel.

\citet{Hebestreit:2013ko} and \citet{Lacey:2013iy} the smoothing results of `BiSeq` to `BSmooth`. Both \citet{Hebestreit:2013ko} and \citet{Lacey:2013iy} provide instances where they claim `BiSeq` gives more 'reasonable' smoothed values than `BSmooth`, but the comparison is selective and limited to a few regions. Furthermore, the comparisons are made using RRBS data, both real and simulated, and, `BSmooth` is designed from WGBS[^rrbs_sim].

[^rrbs_sim]: Both \citet{Hebestreit:2013ko} and \citet{Lacey:2013iy} altered the default `BSmooth` parameters to try to make them comparable to `BiSeq`. \citet{Hebestreit:2013ko} changed the default minimum window size to 80bp but still required at least 20 CpGs per window. \citet{Lacey:2013iy} kept the default minimum window size of 2000 bp but reduced the minimum number of CpGs per window to 50 from the default of 70.


## Statistical properties of $\beta$

Under the binomial model, $M_{i, j} | (M_{i, j} + U_{i, j}, B_{i, j}) = Bin(M_{i, j} + U_{i, j}, B_{i, j})$ ,$\beta_{i, j} = \frac{m_{i, j}}{m_{i, j} + u_{i, j}}$ is an unbiased estiamtor of $B_{i, j}$ with standard error $se(\beta_{i, j}) = \sqrt(\frac{\beta_{i, j}(1 - \beta_{i, j})}{M_{i, j} + U_{i, j}})$ \citep{Hansen:2012gr} (__CHECK WITH TERRY: the standard error should be defined in terms of estimates not parameters, i.e. $\beta$ instead of $B$, correct?__). The natural interpretation of $\beta_{i, j}$ is then an estimator of the average level of methylation at the $i^{th}$ locus in the $j^{th}$ sample. In this section I discuss this interpretation and statistical properties of this estimator.

### Marginal distributions

The average level of methylation $B_{i}$ varies widely across the genome. Much of this variation is associated with different genomic elements. For example, most CpGs in CGIs are unmethylated, which means that the corresponding $\beta$-values are close to zero. In contrast, CpGs outside of CGIs are typically methylated, which means that the corersponding $\beta$-values are closer to one. In humans, __WHAT PERCENTAGE__ of CpGs are in CGIs, which leads to a bimodal distribution of $\beta$-values for most human genomes (__CITE and SHOW PLOTS__) .

* __Show a bunch of plots of $\beta$-values for different samples__

This bimodality has led to several researchers modelling the distribution of the $\beta$-values by the Beta distribution. For example, \cite{Hebestreit:2013ko, Lacey:2013iy} simulate the "true" methylation level in each sample from a Beta distribution and both \citet{Feng:2014iq, Sun:2014fk} assume a Beta distribution as the prior distribution of the $\beta_{i, j}$ in their empirical Bayes models of bisulfite sequencing data. 

The Beta distribution is a flexible 2-parameter distribution on $[0, 1]$. It can be unimodal, "U"-shaped or "J"-shaped, depending on the choice of parameters. The Beta distribution also includes the Uniform and arcsine distributions as special cases (__Source: [http://en.wikipedia.org/wiki/Beta_distribution](http://en.wikipedia.org/wiki/Beta_distribution)__).

In the above examples, the researchers are modelling the distribution of $B_{i, j}$ or $\beta_{i, j}$ __within__ each sample, i.e., over $i$, by a Beta distribution. As I argue in __SECTION__, the distribution of $\beta_{i, j}$ __between__ samples, i.e., over $j$, is of more relevance and interest when performing inference.

* __Show a bunch of different parameterisations of the Beta distribution__

#### Distribution of $\beta$-values within samples {-}

* Stratify by genomic elements

#### Distribution of $\beta$-values between samples {-}

* Stratify by genomic elements
* Means and variations of these distributions and how these relate to differential methylation

### Correlations

Many researchers have observed that DNA methylation is spatially correlated along the genome, (e.g. \cite{Eckhardt:2006gh, Cokus:2008fc, Li:2010fb, Hansen:2011gu, Hebestreit:2013ko, Wang:2011cw, Pedersen:2012vl, Lacey:2013iy, Sofer:2013bk, Liu:dy, Lyko:2010drm, Landan:2012kp, Lister:2009hy}). Just as we can explore the distribution of DNA methylation levels within a sample or between samples, so too can we explore the correlations of DNA methylation levels within a sample or between samples.

#### Correlation of DNA methylation within samples {-}

There are two levels of within-sample correlations of DNA methylation. The lower level is the dependence of DNA methylation at loci on the same DNA fragment, which I call _co-methylation_. I have performed the most comprehensive study to date of this type of within-sample correlation in __CHAPTER__. The higher level is the correlation of the $\beta$-values along the genome. A variety of approaches have been used to estimate these '$\beta$ correlations', although the methods have not always been clearly documented. I clarify these methods and extend these results in __CHAPTER__.

#### Correlation of DNA methylation between samples {-}

There has been less research on between-sample correlations of DNA methylation levels. The most frequently reported between-sample correlation is $cor(\{\beta_{i, j}\}_{i = 1}^{i = N_{loci}}, \{\beta_{i, j'}\}_{i = 1}^{i = N_{loci}})$, which is the correlation of the $\beta$-values for a pair of samples. This has been reported as evidence for the "concordance" or "replicability" of methylation levels for biological replicates (e.g. the same cell line after 4 or 5 cell passagings, \citet{Lister:2009hy}, and three different colon tumour samples \citet{Hansen:2011gu}) and technical  replicates (e.g. different technologies and assays \citet{Zhang:2013uu, Stevens:2013hv, Hansen:2011gu}).

A more detailed measure is the correlation of $\beta$-values for a particular pair of methylation loci between a set of samples. __TODO: Review literature of "co-methylation".__


### Bias

The natural interpretation of $beta_{i, j}$ will be biased if the probability of sequencing a fragment with a methylated site is different from the probability of sequencing a fragment with an unmethylated site. It is very likely that this bias exists but it is difficult to assess without a carefully designed experiment. __TODO: Is there any evidence that this bias exists? Is CGI dropout (as see in Sue Clark's lab) evidence for this?__

One reason to believe that this bias exists is that there is a well-documented GC-bias with Illumina sequencing, whereby fragments with a low or high GC-content are undersampled (__CITE__). In bisulfite sequencing, highly methylated fragments will have more _C_s and fewer _T_s than a lowly methylated version of the same fragment and will therefore have a higher GC-content. Depending on the rest of the sequence in the fragment, this could result in a highly methylated fragment being undersampled (very high GC-content of the entire fragment) or a lowly methylated fragment being undersampled (very low GC-content of the entire fragment).

### Transformations

$\beta$-values are the _de facto_ standard unit for reporting methylation levels due to their natural interpretation. However, they are not necessarily the best unit for statistical inference. This is because a $\beta$-value is an estimate of a proportion and there are a well-known challenges when working with proportion data, such as:

1. Proportions are bound between 0 and 1, inclusive.
2. The estimate of the standard error depends on the estimate of the mean (i.e., $\beta$), through $se(\beta) = \sqrt{\frac{\beta (1 - \beta)}{m + u}}$. Taking the derivative of this, we see that the maximum standard error, $\sqrt(\frac{0.25}{m + u})$, occurs at $\beta = 0.5$ and the minimum standard error, $0$, occurs at $\beta = 0, 1$.
3. We need to more than just the $\beta$-value to have a sense of how precise an estimate it is. Essentially, we need to also know the sequencing coverage of the methylation loci. Consider two CpGs, one with $m = 1, u = 3$ and the other with $m = 100$ $u = 300$. Both CpGs have $\beta = 1/4$ but the second CpG is measured with much greater precision. Assuming the binomial model, the first CpG has $se(\beta) = \sqrt{\frac{1/4 \times 3/4}{4}} = 0.22$ whereas the second CpG has $se(\beta) = \sqrt{\frac{1/4 \times 3/4}{400}} = 0.02$.

To address (2), proportion data are often tranformed via a variance stabilisation transformation. The aim is to make the variance (approximately) independent of the mean. Popular variance stabilisation transformations include:

* The arcsine transformation, $\arcsin{\sqrt{\beta}}$ \citep{ANSCOMBE:1948bw}. A small value, e.g. $0.5$, is added to $m$ and $u$ to avoid $\beta = 0, 1$. 
* The "averaged arcsine" transformation, $\arcsin{\sqrt{\frac{m}{m + u + 1}}} + \arcsin{\sqrt{\frac{m + 1}{m + u + 1}}}$ \citep{Freeman:1950bh}. This transformation does not have a unique inverse \citep{Nunes:2008vj}.

The use of variance stabilising transformations for proportion data seem to have fallen out of favour (e.g. [http://www.esajournals.org/doi/full/10.1890/10-0340.1](http://www.esajournals.org/doi/full/10.1890/10-0340.1)). The noq favoured approach is generalised linear models (__CITE__), in particular the logistic regression model (__CITE__?).


### Regression models

* Logistic regression
	* M-values 
	* The logit transformation, $\text{logit} (\frac{\beta}{1 - \beta})$. A small value, e.g. $0.5$, is added to $m$ and $u$ to avoid $\beta = 0, 1$.
	* The probit transformation, $\Phi^{-1}(\beta)$, where $\Phi$ is the standard Normal cumulative distribution function. A small value, e.g. $0.5$, is added to $m$ and $u$ to avoid $\beta = 0, 1$.
* Beta regression

\cite{Du:2010dc} show improved inferences when using Logit transformed $\beta$-values, which they call M-values. This transformation reduces the heteroscedasticity of $\beta$-values near 0 and 1.


## General TODOs
* Use lower case letters for estimates or hatted versions, e.g. $m$ vs. $\hat{M}$?
* "methylC-seq" or "WGBS" or ???
* "bisulfite-sequencing" or "bisulfite sequencing"
* Notation abuse: e.g. m is defined with respect to m-tuples but also in terms of $m_{i}$. Similarly, $M$ is used to represent methylation patterns but also in terms of $M_{i}$. Perhaps use different typefaces to distinguish them?
* Autocorrelation or correlation?
* Link to all software and corresponding publications
* Decide how to give the total number of a loci and samples in a consistent way. Currently usign $N_{loci}$ for loci and $n$ for samples, which is inconsistent.
* Decide how to denote coverage. Currently using $U + M$ rather than $N$ because $N$ and $n$ are already being used to represent multiple variables.
* Method descriptions are often ambiguous or missing in details. The majority explanations favour words over mathematics and only __WHICH PAPERS__ provide software that implements their analysis methods.
* The Binomial model is really a conditional Binomial model, where the condition is on the sequencing coverage, $M_{i} + U_{i}$. That is, $M_{i} | (M_{i} + U_{i}, B_{i}) = Bin(M_{i} + U_{i}, B_{i})$. __DISCUSS WITH TERRY: What is the correct way to write this Binomial model?
* Discuss distribution of coverage? Perhaps in context of sequencing bias?


