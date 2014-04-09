# A statistical framework for analysing bisulfite sequencing data

## Notes from outline.md
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
		* Lister _et al._ compute both within-read and across-read correlations or dependencies of DNA methylation.
		* Akulenko, R., _et al._ compute gene-level $\beta$ values and compute the correlations between pairs of genes across samples as a function of genomic distance.
	* Dependence/correlation over $j$ (?)
		* Akulenko, R., _et al._ compute gene-level $\beta$ values and compute the correlations between pairs of genes across samples as a function of genomic distance.

## Chapter overview
In this chapter I set out a statistical framework for analysing bisulfite sequencing data. I begin by explaining the various levels of stochasticity in a methylC-seq experiment. Following this, I set out the mathematical notation that I use throughout my thesis and, using this notation, formalise some of the concepts introduced in the previous section. 

I also use this framework to formulate in statistical terminology the key variables and common questions in studies of DNA methylation and compare previous estimators of these key variables.

## One sample
As described in __SECTION__, most analyses to date of WGBS data have focused on the "sample-average" level methylation at individual loci. For this reason, I initially describe the framework in terms of single loci. However, bisulfite sequencing data contain much more information than provided by these univariate, marginal summaries of the data. I introduce _m-tuples_ that can summarise some of this extra information and which I use extensively in __CHAPTERS__.

To begin, I consider the simple experiment of performing methylC-seq on a single sample.
Even with only a single sample, this experiment has a hierarchical structure. I separate this structure into two stages, as illustrated in __CARTOON FIGURE OF METHYLC-SEQ EXPERIMENT__:

1. Pre-sequencing
2. Post-sequencing

Once I have described the framework for a single sample, I extend it to allow for multiple samples. The ideas here are simple, although when extended to their full generality the notation becomes messy. I address complications and limitations of this framework at the end of this sub-section. 

### Single locus analysis

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

Note that in __FIGURE__ I haven't constructed a 2-tuple using the first CpG and the last CpG. In general, I focus on m-tuples where the m methylation loci are neighbours. That is, I focus on m-tuples where the _number of intervening loci_ is zero ($NIL = 0$). There are 3 reasons for this:

1. Quantity: From a sequence containing $M$ methylation loci there are $M - \text{m} + 1$ m-tuples when we restrict ourselves to those m-tuples with $NIL = 0$. In contrast, if we allow $NIL \geq 0$ then there are $\binom{M}{\text{m}} ~ $ m-tuples. (__TODO__: Describe how these terms grow asymptotically).
2. Interpretability: Discussed in __SECTION__
3. Measurability: We cannot observe m-tuples where the methylation loci are far apart due to the read length limitations of the Illumina sequencing technology. This is true even when $NIL = 0$ but is especially the case if we allow $NIL \geq 0$.

When I refer to m-tuples I implicitly mean $NIL = 0$; I will explicitly use the notation $NIL \geq 0$ when I wish to make clear that there may be intervening methylation loci in the m-tuple.

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

By "summing" over the number of reads containing the m-tuple $(i, i + 1, \ldots, i + \text{m} - 1)$ we obtain the number of reads contain each methylation pattern. Here are the definitions for $\text{m} = 2$:

\begin{align*}
	mm_{(i, i + 1)} &= |\{z: z \in \mathcal{R}_{(i, i + 1)}, z = (1, 1)\}| \\
	mu_{(i, i + 1)} &= |\{z: z \in \mathcal{R}_{(i, i + 1)}, z = (1, 0)\}| \\ 
	um_{(i, i + 1)} &= |\{z: z \in \mathcal{R}_{(i, i + 1)}, z = (0, 1)\}| \\ 
	uu_{(i, i + 1)} &= |\{z: z \in \mathcal{R}_{(i, i + 1)}, z = (0, 0)\}| 
\end{align*}

The definitions for $\text{m} > 2$ follow in the obvious manner.

We can extend the $\beta_{i}$ values to m-tuples, although the intuitive interpretation is somewhat lost. Here are the definitions for $\text{m} = 2$:
\begin{align*}
	\beta_{(i, i + 1)}^{mm} &= \frac{mm_{(i, i + 1)}}{mm_{(i, i + 1)} + mu_{(i, i + 1)} + um_{(i, i + 1)} + uu_{(i, i + 1)}} \\
	\beta_{(i, i + 1)}^{mu} &= \frac{MU_{(i, i + 1)}}{mm_{(i, i + 1)} + mu_{(i, i + 1)} + um_{(i, i + 1)} + uu_{(i, i + 1)}} \\
	\beta_{(i, i + 1)}^{um} &= \frac{UM_{(i, i + 1)}}{mm_{(i, i + 1)} + mu_{(i, i + 1)} + um_{(i, i + 1)} + uu_{(i, i + 1)}} \\
	\beta_{(i, i + 1)}^{uu} &= \frac{UU_{(i, i + 1)}}{mm_{(i, i + 1)} + mu_{(i, i + 1)} + um_{(i, i + 1)} + uu_{(i, i + 1)}}
\end{align*}
The definitions for $\text{m} > 3$ follow in the obvious manner.

__TODO__: Change definition of $M$, $m$, etc. when $m = 1$
__TODO__: Mathematically define: reads containing multiple methylation loci; m-tuples; `counts` of methylation patterns.


* Extend $z: z \in \mathcal{G}_{i}$ to $z: z \in \mathcal{R}_{(i, i + 1)}$. And note that we do not know which $h$ each read came from, only that they came from the same DNA fragment.

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
To move from a single sample to $n$ samples simply requires an additional subscript, $j = 1, \ldots, n$. For example, $\mathcal{I}_{j}$ is the set of methylation loci in the $j^{th}$ sample and $\beta_{i, j}$ is the beta value for the $i^{th}$ locus in the $j^{th}$ sample. This defines the three levels in the hierarchy of a typical experiment -- individual molecules ($h$), individual methylation loci ($j$) and individual samples ($j$). A fourth level is how the samples relate in terms of an outcome of interest, such as phenotype. This fourth level might be defined up-front, such as in a designed experiment looking for differences in methylation between pre-defined groups of samples, or the aim of the experiment might be to discover this level.

A common experiment of the first type is the two-group design in which $n_{1}$ samples are from group $1$ and $n_{2}$ samples are from group $2$ ($n_{1} + n_{2} = n$). This can be represented by a design matrix $X = [X_{j}]$, where $X_{j} = 1$ if the sample is from group $1$ and $X_{j} = 0$ if the sample is from group $2$. 

More generally, we can study multi-group experiments via a suitable definition of the design matrix $X$. We can also include covariates by allowing $X_{j}$ to be a row vector, $X_{j} = (x_{1, j}, \ldots, x_{P, j})$, where $x_{p, j}$ encodes the information on the $p^{th}$ covariate for the $j^{th}$ sample.

The second type of experiment is one that aims to cluster individuals based on methylation patterns. This can be thought of as an experiment where the outcome of interest is unknown and we wish to discover it by lumping together samples with "similar" methylation patterns. We might also wish to then correlate these clusters with some external information such as disease-status. 

### Some complications for $n$ samples

In addition to the complications of the previous section, we now have sample-to-sample variation. For example, each individual has their own set of methylation loci, that is, $\mathcal{I}_{j}$ differs for all $j$. Furthermore, sequencing coverage varies from sample-to-sample. This means that even if the samples have exactly the same $\mathcal{I}_{j}$, e.g. the samples are genetically identical, each sample will have a different set of loci with "sufficient" sequencing coverage. Loci without sufficient coverage are effectively missing data.

In practice, we might choose to study $\mathcal{I}^{common} = \cap_{j} \mathcal{I}_{j}$ or some other combination of the $\mathcal{I}_{j}$, such as all methylation loci present in at least some fraction of the $n$ samples. 

A conservative analysis might only analyse those loci where at least some fraction of the $n$ samples have sufficient sequencing coverage. A less conservative analysis might try to impute the missing values based on methylation levels at neighbouring loci (__bsseq__).

#### Complications that aren't merely notational {-}

* $n$ is typically small while $N_{loci}$ is typically large



## Other people's models defined in my framework

# General
* "methylC-seq" or "WGBS" or ???
* "bisulfite-sequencing" or "bisulfite sequencing"
* Notation abuse: e.g. $m$ is defined with respect to m-tuples but also in terms of $m_{i}$. Similarly, $M$ is used to represent methylation patterns but also in terms of $M_{i}$. Perhaps use different typefaces to distinguish them?