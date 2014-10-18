# Paper reviews

These are my reviews of important papers in the field. I do not want to reproduce these verbatim in my thesis but rather summarise these down to key points of similarity and difference between the methods.

## \cite{Eckhardt:2006gh}

\cite{Eckhardt:2006gh} studied DNA methylation patterns on chromosome 6, 20 and 22 in 43 human samples from 12 different tissues. Each of the 43 samples was made up from a pool of up to 3 sex- and age-matched samples.

This study was very labour-intensive work and low-throughput by contemporary standards because it predated the era of cheap, "high-throughput" sequencing (e.g. Illumina sequencing). Instead, the authors used Sanger sequencing of bisulfite-converted PCR amplicons. Briefly, they designed primers for bisulfite-treated DNA from chromosomes 6, 20 and 22. Some of these PCR amplicons were then subcloned into a vector. The non-cloned PCR amplicons and the cloned PCR-amplicons (up to 20 clones per amplicon) were Sanger sequenced using ABI 3730 capillary sequencers. All this results in a fairly sparse sample by modern standards from each chromosome and sample.

Nonetheless, they were perhaps the first to study co-methylation, which they defined as "the relationship between the degree of methylation over distance". The method by which they did this isn't clearly stated in the paper and what follows is my interpretation of what is reported. For a variety of distances between 0 and 20 kb they sampled 25,000 pairs of sequenced fragments (e.g. $z_{h, i}$ and $z_{h', i'}$) where two CpGs were separated by the given IPD (i.e. $pos_{i'} - pos_{i}$ = IPD). For each pair at that distance they recorded whether the methylation measurements were identical, that is, both methylated or both unmethylated (i.e. $\mathbf{1}_{z_{i, h} = z_{i', h'}}$). Then, for each distance, they computed the percentage of pairs that were in agreement, i.e. $\frac{1}{25000}\sum_{25,000 pairs}\mathbf{1}_{z_{i, h} = z_{i', h'}}$.

__TODO: Discuss interpetation with Terry__
__TODO: What are the grey and blue dots in Fig. 3c?__

Based on this analysis they concluded that there was "a significant correlation for co-methylation over short distances ($\leq$ 1,000 bp), [but] it deteriorated rapidly for distances $>$ 2,000 bp".

__TODO: Fix indicator function__

A useful quote regarding mosaicism and epipolymorphism:

> we conclude that the majority (490%) of the observed heterogeneous methylation is caused by mosaicism, although we cannot exclude the additional possibility of heteroge- neous tissue sampling.



## \cite{Cokus:2008fc}
\cite{Cokus:2008fc} published the BS-seq protocol for performing WGBS. Most of the results in this paper used BS-seq of wild-type _Arabidopsis Thaliana_ and DNMT mutants. They also report limited results from low-coverage BS-seq of mouse germ cells but these are not relevant to my discussion.

The authors used the simple $m$ and $u$ read-counting estimators of $M$ and $U$, subject to some filtering and other post-processing of the reads. Most analyses were restricted to loci with at least 5x sequencing coverage.

They visualised the levels of DNA methylation across a variety of genomic elements such as protein coding genes, pseudogenes and repeats. They did not perform an analysis of differential methylation. The authors were also interested in the dependence structure of DNA methylation. For the most part they did this by studying the $\beta$-values but they also performed a "within-read" analysis.

To look at the correlation of "average" methylation, they performed an autocorrelation analysis of the $\beta$-values. It is not clear from the paper whether they restricted these analyses to pairs of methylation loci with $NIL = 0$ or whether they allowed $NIL \geq 0$. They did separate autocorrelation analyses for CpG, CHG and CHH methylation.

To look at the dependence of methylation events from the same DNA fragment, they estimated the "within-read probability of additional methylation of CHH sites within a given distance from a methylated CHH site". They only report the "within-read" results for CHH methylation. This amounts to estimating the conditional probability $Pr(Z_{h, i + s} | Z_{h, i})$ and plotting it as a function of genomic distance between the two loci, $pos_{i + 1} - pos{i}$. It is not clear from the methods whether $s = 1$, i.e. the conditional probability that the __next__ CHH is methylated (which amounts to using pairs of CHH loci with $NIL = 0$), or whether $s \geq 1$, i.e. the conditional probability that __any__ upstream CHH is methylated (which amounts to using pairs of CHH loci with $NIL \geq 0$). Furthermore, these results are limited by the short read lengths (31nt) available at the time of the study.

The authors reported a 10 bp periodicity of CHH methylation, which could be detected in both the $\beta$-value and "within-read" analyses. They noted that this periodicity is consistent with one helical DNA turn and that the mammalian Dnmt3a, which is the homologue to the main enzyme controlling asymmetric (e.g. CHH) methylation in _Arabidopsis_, _DRM2_, had recently been shown to methylate CpG cites $~8-10$ nucleotides apart.

Autocorrelation analyses of the $\beta$-values also identified a periodicity of $~167$nt for CpG, CHG and CHH methylation, although it was strongest for CHG methylation. The authors noted that this "is similar to, but slightly shorter than, estimates of the average spacing of nucleosomes in plant chromatin". From this they postulated that nucleosomes or histone modifications might dictate access to the DNA by DNMTs and that methylated DNA might be more "compact" than unmethylated DNA.

Both the $~10$nt and $~167$nt periodicities can be seen by eye and were confirmed by Fourier analysis of the $\beta$-autocorrelation function.


## \cite{Lister:2008bh}

\cite{Lister:2008bh} published the methylC-seq protocol for performing WGBS. Like \cite{Cokus:2008fc}, \cite{Lister:2008bh} studied wild-type and DNMT-defective mutant _Arabidopsis Thaliana_ strains. Also like \cite{Cokus:2008fc}, much of the analysis is by visual comparison of summarised data.

The authors used the simple $m$ and $u$ read-counting estimators of $M$ and $U$, subject to some filtering and other post-processing of the reads. To identify methylcytosines in each sample, the authors used a simple Binomial model, $Bin(m_{i, j} + u_{i, j}, B_{i, j})$ . They defined $\epsilon_{j}$ as the "error rate" for the $j^{th}$ sample, which is the sequencing error rate + the bisulfite non-conversion rate and is assumed constant across all loci. They estimated $\epsilon_{j}$ from reads aligned to the unmethylated chloroplast genome of each sample and obtained values between $1.3-3.2\%$ per sample.

For each locus and sample the authors tested the hypothesis $H_{0}: B_{i, j} = \epsilon{j}$ vs. $H_{1}: B_{i, j} > \epsilon{j}$. This can be thought of as testing the hypothesis "what is the minimum $m_{i, j}$ I would have to observe to conclude that these weren't all due to "errors"" (__CHECK WITH TERRY, in particular is the alternative one- or two-sided?__). For each sample, all loci were ranked by the false discovery rate-adjusted (FDR-adjusted) P-value from this test (__Not clear from text which FDR procedure was used__). All loci with an FDR < 0.05 were called as methylcytosines.

The authors did not perform any statistical tests of differential methylation nor did they report any results on the dependence structure of DNA methylation.

## \citet{Lister:2009hy}

\cite{Lister:2009hy} was a landmark paper in the study of DNA methylation in humans from WGBS data. The authors studied 4 samples: 2 cell types (IMR90 and H1) and 2 biological replicates per cell type. However, most of the results reported in \cite{Lister:2009hy} were from analyses of data pooled across biological replicates. This completely ignores all biological variability and in general isn't a good idea. In the following description of their statistical analyses, all references to "samples" means "pooled samples".

The authors used the simple $m$ and $u$ read-counting estimators of $M$ and $U$, subject to some filtering of the reads. The authors used the Binomial model from \cite{Lister:2008bh} to identify methylcytosines.  The bisulfite non-conversion rate was estimated from a spike-in control that was included with each sample for sequencing[^spike_in]. It is not clear how they estimated the sequencing error rate. They reported that $\epsilon_{\text{IMR90}} = 0.005$ and $\epsilon_{\text{IMR90}} = 0.0024$. All loci with an FDR < 0.01 were called as methylcytosines.

[^spike_in]: The spike-in control was unmethylated cl857 Sam7 Lambda DNA. This "lambda phage" is a common choice of spike-in control in bisulfite-conversion assays.

With this set of methylcytosines they sought to identify differential methylation and partially methylated domains (PMDs). They used a two-stage analysis to identify differential methylation:

1. Identify differentially methylated cytosines (DMCs)
2. Identify differentially methylated regions (DMRs)

The first step used a Binomial model for each sample, namely $Bin(m_{i, \text{IMR90}} + u_{i, \text{IMR90}}, B_{i, \text{IMR90}})$ and $Bin(m_{i, \text{H1}} + u_{i, \text{H1}}, B_{i, \text{H1}})$. The authors then used a two-tailed Fisher's exact test of the hypothesis $B_{i, \text{IMR90}} = B_{i, \text{H1}}$ vs. $B_{i, \text{IMR90}} \neq B_{i, \text{H1}}$. They performed this test at all loci that were called as methylcytosines in at least one of the cell types and that had at least 3 reads in at least one of the samples. Cytosines were called as differentially methylated if the FDR-adjusted P-value was less than 0.05.

The second step was to identify DMRs, that is, regions containing multiple cytosines that display differential methylation between the IMR90 and H1 samples. In fact, \cite{Lister:2009hy} only sought to "find regions of the genome enriched for sites of higher levels of DNA methylation in IMR90 relative to H1, as identified by Fisher's Exact Test" \cite[Supplementary Material, pp. 26]{Lister:2009hy}. This was only performed for CpG methlylation loci.

The authors used a heuristic approach based on a 1kb sliding window approach and 100bp step size. If the window contained at least 4 differentially methylated cytosines then it was extended in 1kb increments until a 1kb increment was reached that contained less than 4 differentially methylated cytosines. Once the extension had terminated, a region was declared to be a DMR if it contained at least 10 differentially methylated cytosines and was at least 2kb in length. The authors did not report a sensitivity analysis of any of these parameters.

A sliding window approach was also used to identify partially methylated domains (PMDs). This was only performed for CpG methlylation loci. A larger window size of 10kb and step size of 10kb were used. If the window contained 10 mCpGs, each covered by at least 5 reads, and the __average__ $\beta$-value in the region was less than 0.7 then the region was incremented by 10kb. The extension was terminated once the next increment had an average $\beta$-value greater than 0.7 or less than 10 mCpGs and the region was called a PMD.

\cite{Lister:2009hy} also investigated what at first appears to be a "within-read" measure of co-methylation (Fig. 3g, 3h; Supp. Fig. 9). However, on a closer reading I do not think it is a truly "within-read" measure. Specifically, they tabulated the number of mCG (resp. mCHG or mCHH) sites 1-50n bp downstream of a mCG (resp. mCHG or mCHH). I believe that they define methylcytosines based on their binomial test, which means these tabulations are "across-reads" rather than "within-reads". To truly do this as a "within-read" analysis you would tabulate $z_{i, j} = (1, 1)$ for a variety of IPD. The fact that they used 50 bp reads and measured over a 1-50 bp window makes it easy to misinterpret their results as being "within-read".

With that in mind, these results suggest a 8-10 bp periodicity in the co-occurence of methylcytosines. This co-occurence is clearest for CHG and CHH loci in intronic sequences. Other loci and contexts do not show this behaviour or not to the same extend. For example, the graph of mCG co-occurence in exonic sequences is dominated by a 3nt cycle, presumably due to codon structure and selective pressures on coding sequences (Sup. Fig. 9).

To summarise, I take the "co-methylation" results of \cite{Lister:2009hy} with a grain of salt for several reasons:

1. The details of the method are not very clear from the paper. Aside from the "within-read" vs. "across-read" issue, it is not clear whether they use pairs of cytosines with $NIL = 0$ or with $NIL \geq 0$ in their analysis.
2. Estimates of periodicity are based on visual inference from a cubic spline smoothing of the co-occurence patterns. They did not perform the more meaningful Fourier analysis.
3. They only look at the co-occurence of methylcytosines and not the co-occurence of unmethylated cytosines.
4. The number of observations per distance is very small in some contexts. For example, less than 100 observations per IPD are used in the graph of mCHH co-occurence in the "random" context (Sup. Fig. 9).

## \citet{Lister:2011kg}

\cite{Lister:2011kg} is an extension of \cite{Lister:2009hy}. Here the authors studied 15 methylC-seq datas, 4 of which were from previous publications. These samples came a variety of tissues but can be classified as being either from a cell line that is embryonic stem cell (ESC), induced pluripotent stem cell (iPSC), differentiated cell or _in vitro_ differentiated from pluripotent cell (IVD).

\cite{Lister:2011kg} used similar analysis methods as they did in \cite{Lister:2009hy}. Namely, methylcytosines were identified using the Binomial test and DMRs were identified using a sliding window approach. The details of the DMR finder are more complicated due to the larger sample size and multiple comparisons made between the 4 classes of cell type. For two-group comparisons, the average (smoothed) $\beta$-values in each window were tested for a mean difference using a Wilcoxon test. For multi-group comparisons, the Wilcoxon test was replaced by a Kruskall-Wallis one-way analysis of variance. The authors corrected the resulting P-values using the Benjamini-Hochberg method (__CITE__). Putative DMRs were those with an adjusted P-value < 0.01 and were also required to have a mean-difference greater than some threshold.

A sliding window approach was also used to identify partially methylated domains. The authors did not investigate the dependence structure of DNA methylation.

[Mark van de Wiel](http://www.few.vu.nl/~MA.van.de.Wiel/) raised some concerns about the statistical analysis of these experiments [on Nature's website](http://www.nature.com/nature/journal/v471/n7336/full/nature09798.html#comments).

## \citet{Lister:2013et}
\cite{Lister:2013et} used methylC-seq to study 5mc and TAB-seq to study 5hmC in neurons and glial cells from the frontal cortex of human and mouse samples. The authors reported that non-CpG methylation was the dominant form of 5mC in neurons but not in glial cells.

The analysis of 5mC used a different strategy to that of \cite{Lister:2009hy, Lister:2011kg}. Firstly, they did not perform an initial screen for "methylcytosines" using the Binomial test. Secondly, differentially methylated cytosines were identified using a test of the $J \times 2$ contingency table.  Specifically, for the $i^{th}$ CpG a $J \times 2$ table was constructed where row $j$ was the counts of methylated and unmethylated reads for the $j^{th}$ sample, $(m_{i, j}, u_{i, j})$. This table was tested for goodness-of-fit using a root-mean test (http://tygert.com/chi.pdf). DMRs were again constructed using a sliding window approach and were subjected to _post-hoc_ filters.

## \cite{Li:2010fb}
\cite{Li:2010fb} report the methylome of a single sample. They used WGBS to study 5mC from peripheral blood mononuclear cells (PBMCs) from an Asian man whose genome had also been used to create the Han Chinese reference genome.

The authors used the simple $m$ and $u$ read-counting estimators of $M$ and $U$, subject to some filtering of the reads. The authors used the Binomial model from \cite{Lister:2008bh} to identify methylcytosines. As they only had the one sample, many of the analyses were descriptive. For example, they looked at the distribution of $\beta$-values 20 different genomic elements such as CGIs, UTRs and repetive sequences.

The authors looked at the autocorrelation of $\beta$-values across a range of genomic elements. They reported a $~170$nt periodicity in the autocorrelation plot of CpG methylation, similar to that found by \cite{Cokus:2008fc} in _Arabidopsis_. They did not find evidence of any smaller periodicities such as the 8-10 bp periodicity reported by \cite{Cokus:2008fc} and \cite{Lister:2009hy}.

The $~170$nt periodicity can be seen in the autocorrelation plots and was confirmed by Fourier analysis. They found that the autocorrelation of $\beta$-values is stronger when considering pairs of $\beta$-values from the same DNA strand than those from opposite strands. Furthermore, they found this autocorrelation was stronger in some genomic elements than others.

The authors compared the PBMC methylome to the IMR90 methylome from \cite{Lister:2009hy} to identify tissue-specific DMRs (tDMRs). DMR testing was done by forming regions containing 5 CpGs and comparing methylation levels between PBMC and IMR90 using Fisher's exact test, presumably by aggregating read counts across all CpGs in the window. Regions with a P-value $< 10^{-20}$, at least a two-fold difference in methylation between PBMC and IMR90 were declared  tDMRs. Neighbouring tDMRs were joined if they were 'consistent'.

Because the genome of the individual had already been sequenced, this allowed them to study allele-specific methylation (ASM). They extracted all reads overlapping a heterozygous SNP and then looked at methylation levels from reads containing each allele.

## \cite{Hansen:2011gu} and \cite{Hansen:2012gr}

 \cite{Hansen:2012gr} describes a method for identifying DMRs between two groups of samples. This method was first used in \cite{Hansen:2011gu} to identify DMRs in colon cancer tumours. As part of their study, \cite{Hansen:2011gu} performed WGBS on 3 colon cancer samples and their matched normal mucosa, as well as of 2 colon adenomas[^adenoma]. The statistical method they developed, which they call `BSmooth`, are implemented in the R/Bioconductor package `bsseq`. I will describe `BSmooth` and then explain how it was used, in a modified form, in \cite{Hansen:2011gu}.

[^adenoma]: An adenoma is a benign tumour formed from glandular structures in epithelial tissue (__SOURCE: google.com__)

`BSmooth` proposes the following general framework for performing a two-group differential methylation analysis:

1. Smooth raw $\beta$-values.
2. Compute at each methylation locus a test statistic that quantifies the difference between the two groups (in the case of `Bsmooth`, the mean difference in methylation level between the two groups).
3. Find contiguous runs of extreme test statistics; call these putative DMRs.
4. Filter the list of putative DMRs to produce the final list of candidate DMRs. These filters might be biologically-motivated or be designed to remove known false-positive calls from the DMR-calling algorithm.

 \cite{Hansen:2012gr} show that because `BSmooth` explicitly accounts for biological variability it outperforms an analysis using Fisher's exact test of data that has been pooled by group. `Bsmooth` can also make use of the paired design of the experiment described in \cite{Hansen:2011gu}.

### Step 1 {-}

The analysis begins with a table for each sample of the simple $m$ and $u$ read-counting estimators of $M$ and $U$, subject to some filtering of the reads, for each CpG. For each sample, the resulting "raw" $\beta$-values are smoothed in order to remove noise due to low sequencing coverage.

* Is smoothing of the raw $\beta$-values still useful when you have high-coverage sequencing data (__DISCUSS WITH TERRY__)?

A binomial local likelihood smoother was chosen because they model $M_{i, j}$ as $Binom(M_{i, j} + U_{i, j}, B_{i})$ and it is "local" because "methylation levels are strongly correlated across the genome" \citep{Hansen:2012gr}. The authors cite \citet{Eckhardt:2006gh} as evidence that DNA methylation levels are similar at proximal CpGs).

Under the binomial model for $M_{i, j}$, $\beta_{i, j} = \frac{m_{i, j}}{m_{i, j} + u_{i, j}}$ is an unbiased estiamtor of $B_{i, j}$ with standard error $se(\beta_{i, j}) = \sqrt(\frac{\beta_{i, j}(1 - \beta_{i, j})}{M_{i, j} + U_{i, j}})$ (__CHECK WITH TERRY: the standard error should be defined in terms of estimates not parameters, i.e. $\beta$ instead of $B$, correct?__). `BSmooth` assumes that for each sample that the underlying methylation level, $f_{j}(i)$, is a smoothly varying function of the position in the genome, $i$.

A window size for smoothing, $w_{i}$, is defined for each $i$. The default window size is defined as one that contains at least 70 CpGs and is at least 2000kb wide, that is, $w_{i} = (pos_{i}, pos_{i'})$, where:
\begin{equation*}
pos_{i'} = \left\{
  \begin{array}{l l}
    pos_{i} + 2000 & \quad \text{if } NIL(pos_{i}, \ldots, pos_{i} + 2000) \geq 70 \\
    pos_{i} + min(\omega: NIL(pos_{i}, \ldots, pos_{i} + \omega) = 70 - 1 = 69) & \quad \text{otherwise}
  \end{array} \right.
\end{equation*}
Note that the end of the window may not be a CpG if $w_{i} = 2000$ but will be if $w_{i} > 2000$.

The raw $\beta$-values are weighted according to the binomial likelihood and a tricube kernel. The binomial likelihood weights points inversely to their standard error, $se(\beta_{i, j})$, and the tricube kernal gives greater weight to those $\beta_{i, j}$ near the centre of the window.

### Step 2 {-}

\citet{Hansen:2012gr} describe `BSmooth` in the context of a simple two-group linear model. The model could be extended to deal with multiple groups but this has not been done.

Let $X_{j} = 1$ if the $j^{th}$ sample is a case and $X_{j} = 0$ if a control. There are $n_{0}$ controls and $n_{1}$ cases ($N = n_{0} + n_{1}$). `BSmooth` assumes that the samples within each group are biological replicates. Within each window, $\log(\frac{f_{j}(i)}{1 - f_{j}(i)})$ is approximated by a second degree polynomial. `BSmooth` transforms the resulting model parameters(__?__) to fit the model $f_{j}(i) = \alpha_{i} + \beta_{i} X_{i} + \epsilon_{i, j}$. The model allows for locus-dependent variation, $\sigma_{i}^{2}$, that, like $f_{j}(i)$, is assumed to be a a smoothly varying function of the position in the genome, $i$.

__DISCUSS WITH TERRY: Are they fitting the model to $\log(\frac{f_{j}(i)}{1 - f_{j}(i)})$, i.e. a logistic model, and then transforming or fitting directly to $f_{j}(i)$, i.e. a linear model?__

In this model, $\alpha_{i}$ represents the base level of methylation at the $i^{th}$ locus and $\beta_{i}$ is the true difference between the two groups. The $\epsilon_{i, j}$ represent biological variability within the locus-dependent variance, $\sigma_{j}^{2}$.

`BSmooth` fits this model to the __smoothed $\beta$-values__. The parameters are estimated as empirical averages:
\begin{align*}
	\hat{\alpha}_{i} &= \frac{1}{N} \sum_{j = 1}^{N} f_{j}(i) \\
	\hat{\beta}_{i}  &= \frac{1}{n_{1}} \sum_{j: X_j = 1} f_{j}(i) - \frac{1}{n_{0}} \sum_{j: X_j = 0} f_{j}(i).
\end{align*}

To estimate $\sigma_{i}$, firstly the within-group standard deviations are computed. Secondly, in an attempt to "improve precision", these standard deviations are truncated at the $75^{th}$ percentile and then smoothed using a running mean (window size of 101 observations). I denote the resulting estimate of $\sigma_{i}$ by $\hat{\sigma}_{i}$.

__TODO: it's not clear from the text whether smoothing+flooring are done on the within-group sds or the grand sd__.

Finally, a t-statistic is computed at each locus as $t_{i} = \frac{\beta_{i}}{\hat{\sigma}_{i} \sqrt(1 / n_{0} + 1 / n_{1})}$. This tests the hypothesis that the average methylation level at this locus is the equal in both groups, or, formally, $H_{0}: \beta_{i}^{0} = \beta_{i}^{0}$ vs. $H_{1}: \beta_{i}^{0} \neq \beta_{i}^{0}$, where $\beta_{i}^{k}$ is the average methylation level in group $k = 0, 1$. __TODO: Don't use $\beta$ in model definition as this gets confused with $\beta$-values__

The authors recommend that t-statistics are only computed at loci with "some coverage in most or all samples".

### Step 3 {-}

Regions with $|t_{i}| > q_{t}^{0.95}$, where $q_{t}^{0.95}$ is the 95th percentile of the empirical distribution of $\{t_i\}_{i = 1}^{M}$, and where all differences are in the same direction, are called putative DMRs.

### Step 4 {-}

 \citet{Hansen:2012gr} recommend _post-hoc_ filtering of these putative DMRs. These filters include an algorithm to filter out very small DRMs and those that did not display a large mean-difference ($> 0.1$) across the DMR; an algorithm to merge neighbouring DMRs if they were within a pre-defined distance of another DMR and displayed a consistent methylation-difference pattern; and an algorithm to clasify the final list of DMRs according to biological significance.

### `BSmooth` use in \citet{Hansen:2011gu} {-}

All WGBS samples were sequenced at very low coverage ($~5 \times$), which initially motivated the smoothing-approach taken by `BSmooth`. The authors only report looking for DMRs in the tumour-normal comparison and not the normal-adenoma or adenoma-tumour comparisons. However, they compared the adenoma methylation profiles to those of the tumour and normals at several DMRs.

Unlike in 'normal' samples, DNA methylation is highly variable in tumour samples. This means that in a two-group experiment, the within-group variation in the tumour samples will often dwarf the between-group variation. To avoid this, \citet{Hansen:2011gu} use only the normal samples to estimate $\sigma_{i}$. They argue that estimating the standard deviations from the normal samples is equivalent to assuming that the tumour samples are not biological replicates (__DISCUSS ARGUMENT WITH TERRY. Also, do they not do a Welch t-test because this does not account for pairing? Is there a "paired Welch t-test"?__).

 \citet{Hansen:2011gu} used two different window sizes for two different analyses -- a large window size of to detect low frequency changes and a small window size of to detect high frequency changes. The large window size required 500 CpGs per window and a minimum window of 2 kb and the small window size required 70 CpGs per window and a minimum window size of 1 kb[^bsmooth_typo1]. They also used a lower t-statistic cutoff for the low-frequency analysis.

[^bsmooth_typo1]: The supplementary material of \cite{Hansen:2011gu} says that the large window size was 2 kb but I think this is a typo and they meant 20 kb. They add to the confusion in \citet{Hansen:2012gr} where they suggest that the small window size is $\geq2$ kb with at $\geq70$ CpGs and that the large window is $\geq40$ kb and $\geq500$ CpGs.

The analysis also includes a correction factor to allow for high-frequency DMRs within low-frequency DMRs, such as a small hypermethylated region within a larger hypomethylated block. To do this, they modified the linear model to be $f_{j}(i) = \alpha_{i} + \beta_{i}^{s} X_{i} + \beta_{i}^{l} X_{i} + \epsilon_{i, j}$, where $\beta_{i}^{s}$ represents small DMRs and $\beta_{i}^{l}$ represents large DMRs. $\beta_{i}^{l}$ is much more slowly varying than $\beta_{i}^{s}$.

Similarly, the t-statistics are separated into $t_{i} = t_{i}^{s} + t_{i}^{l}$. Any site with differential methylation should have a large value of $|t_{i}|$. However, in order to identify small DMRs there needs to be a correction for $t_{i}^{l}$. \citet{Hansen:2011gu} estimate $t_{i}^{l}$ by linearly interpolating $t_{i}$ across all genomic positions, not just methylation loci. This function is evaluation on a 2kb grid and smoothed using a Huber likelihood model in 25 kb windows[^bsmooth_typo2]. The smoothed function is then evaluated at each methylation loci and these values are taken as estimates of $t_{i}^{l}$. They then use the 'corrected' values of $t_{i}^{s} = t_{i} - t_{i}^{l}$ to identify small DMRs, using the method described in steps 2-4.

[bsmooth_typo2]: \citet{Hansen:2012gr} reports that this is a 50 kb window

### `BSmooth` use in \citet{Hansen:2013eo}
__TODO__

## \citet{Hebestreit:2013ko}

`BiSeq`, published in \citet{Hebestreit:2013ko}, is another R/Bioconductor package that performs DMR detection and uses smoothing of the raw $\beta$-values. It It is designed for use with CpG methylation from RRBS data and takes a 5-step approach to DMR detection:

1. Define CpG clusters
2. Smooth methylation data within CpG clusters
3. Model and test group effect within CpG clusters
4. Apply hierarchical testing procedure:
	* Test CpG clusters for differential methylation and control weighted FDR on clusters
	* Trim rejected CpG clusters and control FDR on single CpG sites
5. Define DMR boundaries

### Step 1 {-}

A CpG cluster is defined as a group of CpGs with the following properies (default parameter values in brackets):

1. All CpGs have 'sufficient' sequencing coverage in 'most' samples (75% of samples have 'sufficient' coverage)
2. All CpGs are within a maximum $IPD$ of the next CpG ($100$ bp)
3. There are at least a minimum number of CpGs in the cluster ($20$ CpGs)

These properties are used to define the cluster boundaries but all $\beta$-values within the cluster boundaries are used in subsequent analyses.

### Step 2 {-}

Within each cluster, raw $\beta$-values are smoothed using a binomial local likelihood with a triangular kernel. This is in contrast to `BSmooth` where smoothing is done across the entire genome and not separately within clusters. Also in constrast to `BSmooth`, `BiSeq` uses a much smaller smoothing window of 80 bp, which is fixed across the genome.

### Step 3 {-}

\citet{Hebestreit:2013ko} use the beta distribution to model the smoothed $\beta$-values, that is, smoothed $\beta_{i} = Beta(\mu_{i}, \phi_{i})$. Under this model,
\begin{align*}
E(\beta_{i}) &= \mu_{i} \\
Var(\beta_{i}) &= \frac{\mu_{i} (1 - \mu_{i})}{1 + \phi_{i}}.
\end{align*}

They claim this is the appropriate distribution since $\beta \in [0, 1]$ for both the smoothed and raw $\beta$-values[^biseq_typo1]. This assumption is debatable, as discussed in __SECTION__.

[^biseq_typo1]: Not $\beta \in (0, 1)$ as claimed in the paper.

\citet{Hebestreit:2013ko} develop a regression model of the smoothed $\beta$-values. `BiSeq` regresses the smoothed $\beta$-values against $P$ covariates via a probit link function, $g(\mu) = \Phi^{-1}(\mu)$, where $\Phi^{-1}(.)$ is the standard Normal cumulative distribution function. Thus, $g(\mu_{i}) = X_{i} \beta_{i}$, where $X_{i,j} = (X_{i, 1}, \ldots, X_{x, P})^{T}$ is a $P$-vector of covariates and $\beta_{i} = (\beta_{i, 1}, \ldots, \beta_{i, P})^{T}$ is a $P-vector$ of regression parameters.

__TODO: Don't use $\beta$ in model definition as this gets confused with $\beta$-values__

Maximum likelihood estimators of $\beta$ and $\phi$ are obtained via numerical methods since no closed form solutions are available[^biseq_typo2]. While $\mu_{i}$ is presumably allowed to vary for each $i$, it is not clear whether $\phi_{i}$ is allowed to vary across loci. However, \citet{Hebestreit:2013ko} do suggest that $\phi_{i}$ might be estimated separately within each group when doing a tumour-normal comparison to allow for the increased within-group methylation variability of tumour samples \citep{Hansen:2011gu}.

[^biseq_typo2]: __CHECK WITH TERRY: is $\beta$ or $\mu$ really the regression parameter? I suppose it doesn't really matter since one is a function of the other.

`BiSeq` performs a Wald test (__CITE__) at each CpG as a preliminary test of differential methylation.

### Step 4 {-}

`BiSeq` uses a hierarchical FDR procedure developed by __CITE BENJAMINI and HELLER__. This procedure first tests for significance at the cluster-level and then at CpGs within significant clusters. This reduces the number of hypotheses since only loci within significant clusters are tested, therefore increasing the statistical power.

As\citeauthor{Hebestreit:2013ko} say, the first step is to "$\ldots$ detect CpG clusters containing at least one differentially methylated location and to control a size-weighted FDR on clusters". For each CpG the P-value from the Wald test is transformed to a z-score, $z_{i}$. Each cluster is then assigned an average z-score, $\bar{Z}_{c}$, where $c = 1, \ldots, C$ is the cluster index. The standard deviation of $\bar{Z}_{c}$, $\hat{sigma}_{\bar{Z}_{c}}$, is estimated using an approach that accounts for the correlation of P-values between two loci $i$ and $i'$. This correlation is estimated from a semivariogram of the $z_{i}$. It is recommended that the variogram is estimated from re-sampled data to ensure that the null hypothesis holds and thus that variance of the z-scores follow the standard Normal distribution.

For each cluster, a P-value is obtained as $\Phi^{-1}(\frac{\bar{Z}_{c}}{\hat{sigma}_{\bar{Z}_{c}}})$. Those clusters with a P-value satisfying the Benjamini-Heller FDR criteria are retained for further investigation.

The second step is to trim the clusters of CpGs that are not themselves differentially methylated. A P-value is computed for each CpG within the significant clusters; this P-value is conditional on the cluster being rejected. These conditional P-values are then subjected to a secondary FDR criteria, proposed by Benjamini and Heller, and those sites with P-values greater than the cutoff are trimmed from the clusters.

__DISCUSS WITH TERRY: I don't fully understand how this FDR procedure works. E.g. it would surely miss any significant CpG within an insignificant cluster, rigth?__


### Step 5 {-}

The result of steps 1-4 are significantly differentially methylated CpGs. While `BiSeq` uses clusters to find these DMCs, these clusters are not themselves the DMRs. Rather, DMRs are defined as adjacent DMCs that occur within the same cluster and that display consistent methylation differences, that is, the DMCs are not switching from positive to negative differences.

__NOTE TO SELF: I like the FDR procedure but I'm not convince beta regression is necessary. Perhaps we can get away with the cheaper linear/logistic regression but make use of the cluster + FDR approach?__

## \citet{Akalin:2012cm}
\citet{Akalin:2012cm} developed `methylKit`, an `R` package for processing bisulfite-sequencing data (__LINK TO SOFTWARE__). The methods describedin the paper are applicable to WGBS, although all examples in the paper use RRBS data and the accompanying website for the software suggests that WGBS data might be too large for the software to handle [`v0.9.2` https://code.google.com/p/methylkit/](https://code.google.com/p/methylkit/). `methylKit` contains many useful utility functions for processing and visualising bisulfite-sequencing data.

\citet{Akalin:2012cm} does not describe a method for _identifying_ DMRs, rather it describes two methods for _testing_ pre-defined regions

\citet{Akalin:2012cm} describe two methods implemented in `methylKit` for testing for differential methylation at either individual cytosines or at __pre-defined__ regions (using some aggregated/average measure of methylation in the region). It does not describe a method for _identifying_ DMRs, such as post-hoc combining DMCs into regions.

The tests of differential methylation are logistic regression of the raw $\beta$-values or Fisher's exact test of a table of $m$ and $u$ counts. The default in `methylKit` is to use logistic regression if there are multiple samples per group and Fisher's exact test otherwise. \citet{Akalin:2012cm}  

Once DMC/DMR testing has been performed, P-values are corrected for multiple hypothesis testing. `methylKit` implements both the sliding linear model (SLIM), which adjusts for correlations amongst nearby P-values \citep{Wang:2011cw}, and the Benjamini-Hochberg (__CITE__) method for P-value correction.

## \citet{Feng:2014iq}
\citet{Feng:2014iq} describe an empirical Bayes hierarchical model to identify DMCs from bisulfite sequencing data. The method is implemented in the R/Bioconductor package, `DSS` (__CITE__).

Let $k = 1, \ldots, K$ denote which group the $j^{th}$ sample belongs to. \citep{Feng:2014iq} describe a case-control experiment (i.e. $K = 2$), and only consider CpGs. `DSS` uses the following beta-binomial model[^dss_notation]:

[^dss_notation]: Note that in order to make the notation in this section consistent with the general framework described in this chapter, there are subtle differences to the index parameters used in \citet{Feng:2014iq}.

\begin{align*}
M_{i, j, k} | B_{i, k}, M_{i, j, k} + U_{i, j, k} &= Binomial(M_{i, j, k} + U_{i, j, k}, B_{i, k}) \\
B_{i, k} &= Beta(\mu_{i, k}, \theta_{i, k}) \\
\theta_{i, k} &= log-Normal(m_{0, k}, r_{0, k}^{2})
\end{align*}

In words, for each CpG the true level of methylation, the prior parameters and the hyperparameters are assumed to be identical within each group.

Parameter estimates are obtained by the following algorithm:

1. For each CpG, estimate the dispersion parameter, $\theta_{i, k}$, using a method  of moments estimator.
2. Estimate the hyperparameters, $m_{0, k}$ and $r_{0, k}^{2}$, as the mean and variance, respectively, of the empirical distribution of the logarithm of the $\theta_{i, k}$.
3. Estimate the group-wise mean level of methylation, $\mu_{i, k}$, by $\hat{\mu}_{i, k} = \sum_{j: j \in group_{k}} \frac{m_{i, j, k}}{m_{i, j, k} + u_{i, j, k}}$.
4. For each CpG, compute relevant quantities, such as the mean and standard deviation, of the conditional posterior distribution of the dispersion parameters, $Pr(\theta_{i, k} | m_{i, j, k}, m_{i, j, k} + u_{i, j, k}, \mu_{i, k})$. In practice, these quantities are computed using the Newton-Raphson method after plugging in estimates of $m_{0, k}, r_{0, k}^{2}$ and $\mu_{i, k}$.

The effect of the this procedure is to shrink the CpG-, group-wise dispersion estimates, $\theta_{i, k}$, towards the group-wise prior mean, $m_{0, k}$. In this form, each CpG site is allowed to have its own group-wise dispersion although `DSS` also allows the user to select a common dispersion across groups. \citet{Feng:2014iq} do not describe any further moderating of the $\theta_{i, k}$, such as the abundance-dependent trended dispersion estimates available in the differential gene expression analysis software, `edgeR` (__CITE__).

Now that all parameters have been estimated, each CpG is tested for differential methylation. Namely, for a two-group experiment, at each CpG `DSS` uses a Wald test of the hypothesis $H_{0}: \mu_{i, 1} = \mu_{i, 2}$ vs. $H_{0}: \mu_{i, 1} \neq \mu_{i, 2}$. Due to the hierarchical structure of the model, it is not straightforward to derive the null distribution of the resulting test statistics. Instead, \citet{Feng:2014iq} use simulation studies to argue that the test statistics can be safely approximated by the standard Normal distribution.

While \citet{Feng:2014iq} focuses an empirical Bayes model for detecting differential methylation at individual CpGs, the authors propose a simple thresholding algorithm to identify DMRs. This algorithm calls DMRs as regions greater than a specified minimum size (100bp), containing a specified minimal number of CpGs (3) and containing multiple DMCs that satisfy a P-value threshold (at least 80% of DMCs with P-values < a user-specific minimal P-value). As the authors note, there is no accounting for the spatial correlation amongst P-values and say that improved DMR calling is left for future work.

## \cite{Sun:2014fk}

The `MOABS` software, published in \citet{Sun:2014fk}, uses a Beta-Binomial empirical Bayes hierarchical model of methylation at individual CpGs. The model is therefore similar in spirit to that proposed in \citet{Feng:2014iq}, although the model is not as well described[^moabs_errors1]. For example, it is not made explicit how the hyperparameters of the Beta distribution are estimated.

[^moabs_errors1]: The EB description is poorly written and I think is wrong in several places, which makes it confusing. (__DISCUSS WITH TERRY__)

For a two-group experiment, `MOABS` does not perform hypothesis testing to identify DMCs. Rather, for each methylation loci, `MOABS` estimates a $95\%$ posterior probability interval, $PI(a, b)$ (__TODO__ Is this a credible interval in Bayesian parlance?) of the difference in mean methylation levels between group 1 and 2[^moabs_errors2].

[^moabs_errors2]: The paper calls this a confidence interval rather than a posterior probability interval. This seems to mix up the hypothesis-testing and Bayesian frameworks.

They then define a _credible methylation difference_ ($CDIF$) as:

\begin{equation*}
PI(a, b) = \left\{
  \begin{array}{l l}
    a & \quad \text{if } a \geq 0 \\
    0 & \quad \text{if } a < 0 < b \\
    b & \quad \text{if } b \leq 0 \\
  \end{array} \right.
\end{equation*}

The CDIF effectively converts all posterior intervals to a single number. \cite{Sun:2014fk} write, "In practice, $CDIF$ represents the conservative estimation of the true methylation difference, i.e. for $97.5\%$ of chance the absolute value of true methylation difference is greater than or equal to that of $CDIF$." They interpret $CDIF = 0$ to mean that there is no significant difference in methylation betweent the two groups[^moabs_errors3].

[^moabs_errors3]: Again, this seems to mix up the hypothesis-testing and Bayesian frameworks.

DMCs are those cytosines with a $CDIF > \text{cutoff}$, where $\text{cutoff}$ is chosen so as to control the FDR. The FDR-based $\text{cutoff}$ is estimated by permutating sample labels and re-running the analysis (__CHECK WITH TERRY: Necessary to do this permutation multiple times?__).

The procedure could also be used to identify methylation changes at pre-defined regions, rather than individual methylation loci, by aggregating/averaging methylation levels across each region.

For a two-group experiment, `MOABS` uses a hidden Markov model (HMM) to group cytosines into DMRs. Note that all cytosines, not just DMCs, are used by the DMR-detection HMM. To do this, `MOABS` uses a first-order HMM where the observed values are the group-wise difference in mean methylation at each locus, $B_{i}^{2} - B_{i}^{1}$. There are 3 hidden states -- hypomethylation ($B_{i}^{2} - B_{i}^{1} < -\text{cutoff}$), no difference (|$B_{i}^{2} - B_{i}^{1}| < \text{cutoff}$) and hypermethylation ($B_{i}^{2} - B_{i}^{1} > \text{cutoff}$). The $\text{cutoff}$ is the same at that used for identifying DMCs (__DISCUSS WITH TERRY: I think this is true. Is this a good idea?__).

`MOABS` is written in `C++` (__TODO: Is it a command line program or is there an R interface?__) and also provides several other useful functions including:

* A first-order, 2-state HMM to detect hypomethylated and hypermethylated regions from a single sample
* A script to downsample BAM files in regions containing predicted copy number variants (CNVs). \cite{Sun:2014fk} claim that "this will result in CNV bias free `BAM` files for downstream analysis".

## \citet{Stockwell:2014fq}

`DMAP`, published in \citet{Stockwell:2014fq}, does not allow detection of DMCs, rather it focuses solely on identifying DMRs. `DMAP` provides a function to test for differential methylation in sliding windows, _a la_ \cite{Lister:2009hy}, which the authors recommend for analysing WGBS data. For analysing RRBS data, `DMAP` proposes the use of the MspI fragments as the natural unit for calling DMRs. MspI fragments are 40-220 bp long and are identified by scanning the reference genome for the recognition motif, _C'CGG_. `DMAP` includes an algorithm for resolving methylation calls from reads that overlap adjacent MspI fragments.

Within each window, whether it be a fixed-width sliding window or based on MspI fragments, a Fisher's exact test is performed between two samples[^dmap]. If there are more than two samples then the minimum P-value from all pair-wise comparisons is reported as evidence of differential methylation. `DMAP` also provides a $\chi^2$-test for inter-individual variability of DNA methylation in multiple samples.

[^dmap]: It is not clear whether `DMAP` only uses the Fisher's exact test for two __samples__ or also for two __groups__, e.g., by first aggregating all read counts within each group and then doing a Fisher's exact test.

`DMAP` can also perform an F-test to compare DNA methylation between two groups, for example, cases and controls. This test is implemented outside of a regression framework and so does not allow for additional covariates to be included in the model.

Regardless of which test is used, DMRs are called as those regions with a P-values below a given threshold and satisfying other criteria such as a minimum number of CpGs in the window.

## \citet{Lacey:2013iy}

\cite{Lacey:2013iy} details a modeling and inference framework for RRBS data, as well as a method for simulating RRBS data. Here I only describe the modelin framework and reserve the description and discussion of the simulation framework for __CHAPTER__.

The modeling framework is based on a two-group experiment. At each methylation loci the authors fit a logistic regression model of $\beta_i$ against a indicator variable denoting from which group the sample comes. The raw P-values from the $\chi^2$ test of the residuals of this model are then transformed by a procedure dubbed `QUASASS` (quantile adjustment for score and sample size). The authors claim that `QUASASS` adjusts for the differences in read coverage and sample size that would otherwise reduce the sensitivity in low-coverage regions and reduce the specificity in high-coverage regions. __DISCUSS WITH TERRY: What is `QUASASS` actually doing?__. The authors recommend that `QUASASS`-adjusted P-values, rather than the raw P-values, are used in all downstream analyses.

To identify DMRs, \citeauthor{Lacey:2013iy} seek to model the distribution of P-values at neighbouring methylation loci. They use a Uniform Product distribution to model the P-values under the null hypothesis of no differential methylation at each loci. This model assumes that P-values at neighbouring loci are independent, which seems unlikely given the spatial correlation of $\beta$ values.

DMRs are regions where the first methylation loci has a P-value $\leq 0.05$ and are then built up using an iterative procedure. This iterative algorithm continues to add neighbouring methylation loci (up to a user specific threshold with default value of $50$) and computes the P-value of the region under the Uniform Product model. The maximum number of additional loci that retains a region-wise P-value $\leq 0.05$ is defined as the end of the region provided the following conditions are also met:

* That the end site by itself has a P-value $< 0.05$
* All differences in methylation are 'in the same direction'
* All methylation loci in the region are in "sufficiently close proximity"

The final list of putative DMRs is then pruned to remove regions containing only a single loci and overlapping DMRs are merged.

\citet{Lacey:2013iy} state that an advantage of this iterative procedure is that, "DMRs can include internal sites that are not significantly differentially methylated" as well as sites that are missing data.

An R package implementing the algorithm is available from [http://rrbs-sim.r-forge.r-project.org/](http://rrbs-sim.r-forge.r-project.org/).

### Simulation model {-}

__This section to move to chapter on simulating WGBS data.__

It is important to note that \citet{Lacey:2013iy} seek to simulate RRBS data and not WGBS data. This has several implications. Firstly, the RRBS assay enriches for small, CpG-dense regions of the genome such as CGIs. As such the distribution of IPDs is skewed towards zero and is less bimodal than that from WGBS data (__CHECK__), although it still has a long right-tail. Secondly, these CpG-dense regions are enriched for methylation loci with methylation levels near zero. Finally, the sequencing coverage of loci is right-skewed, even bimodal, in RRBS data.

All parameters used in the simulation model were estimated from a single normal myotube cell line (MTCTL2; __CITE Tsumagari et al., 2013__). \citeauthor{Lacey:2013iy} state that this sample is representative of the methylation profiles that they are typically interested in as part of their research emphasis. It is therefore unclear how generalisable these parameters estimates are.

Rather than simulate data at methylation loci from a reference genome, the authors model the distribution of IPDs for CpGs. Specifically, they model IPDs in human RRBS data. For this they use a two-component Normal mixture distribution, where one component models CpG-dense regions of the genome and the second component models CpG-sparse regions of the genome. They further impose a hidden Markov model on the distribution of IPDs so that a small IPD is more likely to be followed by another small IPD, in an attempt to model the clusters of CpG sites observed in RRBS data. It is unclear to me what advantage, if any, modelling the distribution of IPDs has over simply taking the empirical distribution of IPDs from a reference genome.

Unlike `methsim`, \citeauthor{Lacey:2013iy} simulate $\beta$-values and sequencing coverage rather than simulating the reads themselves. They use a mixture of two gamma distributions to model the distribution of log-transformed sequencing coverage. This model is applied independently to bins of methylation loci (bin width = $50$ bp) and thus all loci within each bin have the same sequencing coverage for a given sample.  As this model is independently applied to each bin, it does not quite account for the fact that sequencing coverage is highly correlated at loci separated by short distances ($< 50-200$ bp) due to each read potentially containing multiple methylation loci/covering multiple bins (__TODO: show correlation of sequencing coverage elsewhere in thesis and reference here__). They did, however, model the correlation of sequencing depth at each site across samples (median correlation $= 0.8$) via a Normal copula. The simulated values are then exponentiated (as the model is of logarithmic coverage) and rounded to the nearest integer (as the model produces non-integer variables).

The simulation procedure aims to capture the spatial correlation of $\beta$-values. This correlation was estimated by fitting a Guassian model to the empirical variogram for sites randomly selected from chromosome 11 of the MTCTL2 data. \citet{Lacey:2013iy} report that the variogram shows "a strong correlation for sites in close proximity, decaying to near independence at distances beyong 3000 bp"[^lacey]. The spatially correlated $\beta$-values were then simulated via a two-step process. Firstly, a sample of independent $\beta$-values were sampled from a Beta distribution, with parameters estimated from chr11 of the MTCTL2 data assuming independence of loci. Secondly, correlations were induced via the model proposed by __CITE Zaykin et al. and describe the model__, which uses the correlations estimated from the fitted variogram to induce correlations in the Beta random variates. The parameters of the underlying Beta distribution are then re-estimated to better approximate the distribution of $\beta$-values from the MTCTL2 chromosome 11 data. This algorithm was run several times to arrive at the final parameter values.

[^lacey]: My results show a non-zero correlation beyond 3 kb, which is hardly independence.

They also argue, however, that the correlations of differences in $\beta$-values between two groups will be less than correlations of the $\beta$-values themselves, and provide some evidence for this. To do this they fit exponential variogram models to the raw and `QUASASS`-adjusted P-values from non-differentially methylated loci and report $r \leq 0.4$ for sites with $IPD > 50$ bp and $r < 0.2$ for sites with $IPD > 1000$ bp. It is not clear whether only pairs of sites with $NIL = 0$ where used in this analysis.

\citet{Lacey:2013iy} also propose an algorithm for inducing DMRs in a two-group experiment, where the first $\frac{n}{2}$ samples are cases and the second $\frac{n}{2}$ samples are controls. The DMR-construction algorithm is as follows:

1. Specify the proportion of the methylation loci that are in DMRs, the expected difference in methylation levels for loci in these DMRs ($\delta_{B}$) and the length of these DMRs.
2. Regions satisfying the minimum length requirement are identified from the data and a subset are sampled so that the constraint on the proportion of methylation loci in DMRs is satisfied.
3. Compute the median $B_{i, j}$ for control samples in each region.
	* If the median $B_{i, j}$ in the controls is less than $\delta_{B}$ then the region is assigned as hypermethylated relative to the controls.
	* If the median $B_{i, j}$ in the controls is greater than $1 - \delta_{B}$ then the region is assigned as hypomethylated relative to the controls.
	* Otherwise the region is randomly assigned as hypomethylated or hypermethylated in controls.
4. If the region is assigned as hypomethylated, then the methylation levels of the cases in that region are replaced with $min(B_{i, j} - \delta_{B}, 0.01)$. If the region is hypermethylated, then the methylation levels of the cases in that region are replaced with $max(B_{i, j} + \delta_{B}, 0.99)$.

Under this model all DMRs are of the same length, although they may contain a different number of methylation loci, and all DMRs have the same expected difference in methylation levels between cases and controls.

To summarise, \citeauthor{Lacey:2013iy} propose the following algorithm to simulate RRBS data, including DMRs, for $n_{loci}$ loci and $n$ samples:

1. Simulate $n_{loci} - 1$ IPDs. All samples are assumed to have the same set of resulting methylation loci, i.e. $\mathcal{I}_{1} = \ldots = \mathcal{I}_{n}$.
2. Bin the $n_loci$ loci in 50 bp bins. For each bin, simulate sequencing coverage across the $n$ samples using the Normal copula model of sequencing coverage. Let $n_{i, j} = m_{i, j} + u_{i, j}$ by the simulated sequencing coverage at the $i^{th}$ locus for the $j^{th}$ sample.
3. Simulate spatially correlated methylation levels, $B_{i}$. __The same $B$-values are used for all samples, i.e. $B_{i, j} = B_{i}$.__
4. Construct DMRs by modifying the methylation levels $B_{i, j}$ of cases within these regions.
5. Generate observed $\beta$-values, $\beta_{i, j} = \frac{m_{i, j}}{n_{i, j}}$ using the binomial model, $m_{i, j} = Bin(n_{i, j}, B_{i, j})$, where outside of DMRs $B_{i, j} = B_{i}$ and inside DMRs $B_{i} = B_{i, case}$ or $B_{i} = B_{i, control}$, as appropriate.

## \citet{Akulenko:2013fl}

\citet{Akulenko:2013fl} investigated the co-methylation of genes between samples in a study of $344$ samples from The Cancer Genome Atlas (TCGA). Most of these samples ($317/344$) were breast cancer samples while the remainder ($27/344$) were from matched normal tissue. This study used data from the Illumina Infinium HumanMethylation27k BeadChip (__CITE__), which contains approximately $27,578$ probes that measure single-CpG DNA methylation at CpGs in $14,475$ (__CITE: Original paper not Akulenko__). Most (__HOW MANY__) of these probes are in gene promoters.

After pre-processing of the raw data, each gene was assigned the average $\beta$-value of all $\beta$-values within that gene for each sample. This reduced the number of $\beta$-values from 25,578 to 13,313. After some further filtering to remove probes that were known to be prone to technical artefacts, they computed the Pearson correlation ($r$) of all pairs of gene-level $\beta$-values ($88,611,328$ pairs) using all 344 samples.

This analysis identified $13,643$ pairs of $\beta$-values with $|r| > 0.9$ and $377,547$ pairs of $\beta$-values with $|r| > 0.75$. The overwhelming majority ($377360/377547, 99.95\%$) of pairs with $|r| > 0.75$ were driven by outlier observations, whereby a small number of samples had $\beta$-values greatly different to the rest of the samples.

\citeauthor{Akulenko:2013fl} also investigated how these pair-wise correlations vary as a function of the genomic distance between the two probes in each pair. They reported that "co-methylation level only weakly anti-correlated with genomic distance ($r = -0.29$)", however, this analysis only used a tiny subset of the dataset. Namely, they only performed this analysis for the 74 of the 187 "significant" pairs of $\beta$-values where both genes in the pair were on the same chromosome.

There are several limitations to this study:

1. The low resolution of the 27k array, which leads to summarising a gene's methylation level from measurements of one or two CpGs.
2. They only used a tiny subset of the available data to explore the relationship between co-methylation and genomic distance.
3. Once the "outliers" were removed, almost all the 'co-methylation' disappears. I think these outlier samples need closer investigation rather than just removal.


## \citet{Sofer:2013bk}

The [`Aclust`](http://www.hsph.harvard.edu/tamar-sofer/packages/) software, published in \citet{Sofer:2013bk}, uses generalised estimating equations to detect differential methylation. `Aclust` first clusters sites and then tests for exposure effects on clusters (e.g. case/control status). This in contrast to most other DMR methods, which test for exposure effects at individual methylation loci and then cluster loci. While `Aclust` was designed for analysing methylation array data, in principle it should also work for bisulfite-sequencing data.

`Aclust` uses the following model, where $E_{j}$ denotes the exposure of the $j^{th}$ sample and $X_{j}$, a $p \times n$ matrix, denotes $p$ additional covariates for each of the $n$ samples, such as the batch or age of the sample[^aclust_notation]:

[^aclust_notation]: Note that this differs from the original notation, particularly in the use of index variables.

\begin{equation*}
\beta_{i, j} = \alpha_{i} + E_{j} \alpha_{E_{i}} + X_{j}^{T} \alpha_{X_{i, j}} + \epsilon_{i, j}.
\end{equation*}

This model allows each of the $i = 1, \ldots, n_{loci}$ methylation loci to have a unique baseline methylation value ($\alpha_{i}$), an exposure effect ($\alpha_{E_{i}}$) and covariate effects ($\alpha_{X_{i}}$). Note that the covariates themselves can vary across loci. The errors, $\epsilon_{i, j}$, are assumed to follow a zero-mean distribution. The test of interest is generally of $H_{0}: \alpha_{E_{i}} = 0$ vs. $H_{1}: \alpha_{E_{i}} \neq 0$.

Rather than test every loci, however, `Aclust` first clusters the loci and instead fits a related model to the clusters themselves. Suppose that loci $i = 1, 2, 3$ are determined to be in a cluster, $c$. Then, the model for cluster $c$ is:

\begin{equation*}
\beta_{i, j} = \alpha_{i} + E_{j} \alpha_{E_{c}} + X_{j}^{T} \alpha_{X_{i, j}} + \epsilon_{i, j}.
\end{equation*}

This model allows for each of the $i = 1, 2, 3$ loci in the $c^{th}$ cluster to have a unique baseline methylation value ($\alpha_{i}$) and covariate effects ($\alpha_{X_{i}}$). However, the exposure effect, $E_{c}$, is assumed to be constant for all loci in the cluster. The errors, $\epsilon_{i, j}$, are assumed to follow a zero-mean distribution with a covariance matrix[^aclust_typo]. The covariance matrix could include covariances between errors $\epsilon_{i, j}$ and $\epsilon_{i, j'}$. The test of interest is generally of $H_{0}: \alpha_{E_{c}} = 0$ vs. $H_{1}: \alpha_{E_{c}} \neq 0$.

[^aclust_type]: The $\epsilon$ are defined in two different ways in the paper (p2).

`Aclust` performs a form of agglomorative nested clustering (__CITE Izenman, 2008__) of adjacent methylation loci. Initially, each cluster, $c_{l}$, is comprised of a single methylation loci. Then, clusters are iteratively merged if the _distance metric_ of the two clusters is less than a minimum value, $\bar{\mathcal{D}}$, (default 0.25) and the $IPD$ of the last site of the first cluster and the first site of the second cluster is less than a minimum distance (default 1 kb). The algorithm terminates when no more clusters can be merged.

The distance metric for two methylation loci, $(i, i')$, is $dist(i, i') = 1 - cor(\{(\beta_{i, j}, \beta_{i', j})\}_{j = 1}^{n})$, with the Spearman correlation as the default. The distance metric between two clsuters may use _single_, _average_ or _complete_ distance of the all $dist(i, i')$, where $i$ is from the first cluster and $i'$ is from the second cluster. The _single_ distance requires that only at least one pair of methylation loci have $dist(i, i') < \bar{\mathcal{D}}$, the _average_ distance requires that the mean distance between all pairs of loci is $< \bar{\mathcal{D}}$, and the _complete_ distance requires that the distance between all pairs of loci is $< \bar{\mathcal{D}}$. For a fixed $\bar{\mathcal{D}}$, the _complete_ distance produces smaller clusters than the _average_ distance, than in turn produces smaller clusters than the _single_ distance. The authors recommend _average_ or _complete_ linkage for identifying DMRs, based on their simulation study.

The authors also recommend an initial merging step prior to clustering, which merges all loci within a given $IPD$ (default 99 bp). This step is recommended in order to avoid `Aclust` breaking larger clusters of generally highly correlated loci into multiple smaller clusters due to an intervening methylation loci that is not as correlated with the rest of the loci in the larger cluster. Clusters may be filtered out if they do not contain a minimum number of methylation loci.

Once the final set of clusters is formed, the effect of the exposure variable on each cluster is tested within the GEE framework. `Aclust` uses a robust sandwich variance estimator and requires the user to supply a (possibly mis-specified) working covariance matrix. A Wald test is used to test for the effect of the exposure. Finally, the resulting P-values are corrected for multiple testing via the Benjamini-Hochberg procedure (__CITE__).

### Simulation methodology {-}

\citeauthor{Sofer:2013bk} adapt the methodology of __CITE: Gaile et al.__ to simulate correlated $\beta$-values in a two-group experiment. They simulate Illumina 450k beadchip data with $n_{1} = 40$ cases and $n_{2} = 40$ controls. The simulated data are based on a block re-sampling of data from 539 breast invasive adenocarcinoma samples assayed on the Illumina 450k beandchip obtained from The Cancer Genome Atlas (TCGA). By sampling blocks of CpGs from the same sample rather than sampling individual CpGs, the correlations of the $\beta$-values are preserved.

Firstly, they segment the genome into blocks ($n_{blocks} = 2861$) based on the distance between assayed CpGs, where a new block is formed if the next CpG is more than 10 kb from the previous CpG. Singleton CpGs were added to the nearest block. This segmentation means that methylation levels at CpGs in the same block are correlated, although this is imposed somewhat indirectly by segmenting on $IPD$ rather than the correlations themselves.

Based on the 539 TCGA samples, they then selected a small number of CpGs (10), called targets, that displayed substantial variability across the 539 samples, had "substantial correlation with neighbouring sites" and were not in the same block. They then simulated sample-specific blocks of $\beta$-values by sampling these from the 539 TCGA samples. If a block contained one of the targets, then sampling was done in such a way that one group had that block preferentially sampled from those samples' blocks with a high $\beta$-value at that target and the other group had that block preferentially sampled from those samples' blocks with a low $\beta$-value at that target. This induces differential methylation between the two-groups at that target and, more generally, across several CpGs in the block due to the correlation amongst $\beta$-values. Blocks did not contain targets were sampled from the 539 TCGA samples without any sampling weights.

Due to the correlation structure of the $\beta$-values, it is likely, though not guaranteed, that other CpGs in the blocks containing targets also display differential methylation.

## \citet{Su:2012hl}

\citet{Su:2012hl} describes the `CpG_MPs` software (http://bioinfo.hrbmu.edu.cn/CpG\_MPs) for identifying differential methylation patterns from bisulfite sequencing data. `CpG_MPs` scans the genome for contiguous runs of $\beta$-values $\leq 0.3$ ($\geq 0.7$) and calls these as unmethylated (methylated) "hotspots" (regions)). It then extends these regions by adding neighbouring CpGs, provided that there is no more than 1 Cpg with a $\beta$-values $> 0.3$ ($< 0.7$) in the extension. Finally, `CpG_MPs` computes summary statistics, such as the mean and standard deviation, of the $\beta$-values in these hotspots (__I don't know when these values are even used by `CpG_Mps`__).

`CpG_MPs` also assigns methylated hotspots a value of $1$ and unmethylated hotspots a value of $-1$. Hotspots that overlap across samples are called an overlapping region (OR) and are assigned a score, $-1 \leq u \leq 1$, based on an aggregation of the hotspot scores. Overlapping regions are labelled as "conservatively unmethylated regions" if $u = -1$, DMRs if $ -1 < u < 1$ and "conservatively methylated regions" if $u = 1$. __Then there's a whole bunch of other weird arbitrary cutoffs, definitions and use of Shannon Entropy to further refine these regions__.

## \citet{Liu:dy}

\citet{Liu:dy} identified clusters of CpGs that had correlated $\beta$-values and where the level of methylation in the cluster was influenced by genetic variation (SNPs), which they called these _GeMes_. They used from three studies ($n = 247, n = 91$ and $n= 305$, respectively), where DNA methylation was measured using the Illumina 450k beadchip and each sample was also genotyped using a SNP array. I only describe what they discovered about $\beta$-value correlations and do not discuss the algorithm for finding GeMes or the GeMes themselves.

Spatial correlations of $\beta$-values along the genome across multiple samples are similar to linkage disequilibrium (LD) between SNPs. \citeauthor{Liu:dy} explicitly compared the strength of these correlations in the same individuals. They reported that the $\beta$-value correlations, which started at around $0.5$, were reduced to less than $0.25$ when $IPD > 500$ bp and were more-or-less zero when $IPD > 2000$ bp. In contrast, SNP LD correlations in the same individuals, which also started at around $0.5$, were reduced by half when $IPD > 3000$ bp. Spatial correlations of $\beta$-values used only the top $25\%$ most variably methylated CpGs and were smoothed using a cubic spline. They also drew heat maps, like LD block maps, of $\beta$-value correlations to display their results.

## \citet{Xie:2011cy} and \citet{He:2013cj}

\citet{Xie:2011cy} define _methylation entropy_ (ME) to quantify the heterogeneity of methylation patterns at m-tuples within a single sample. `DMEAS`, published in \citet{He:2013cj} and available from [http://sourceforge.net/projects/dmeas/files/](http://sourceforge.net/projects/dmeas/files/), is software that implements an algorithm to compute ME.

m-tuples can have identical $ME$ but very different average levels of methylation ($\beta$-values) at each of the $m$ loci. For example, an m-tuple that has only fully methylated reads and an m-tuple that has only fully unmethylated reads both have an $ME = 0$ because there is no heterogeneity at either m-tuple. Similarly, m-tuples can have different $\beta$-values but identical $ME$.

ME is a modified version of Shannon entropy and is defined as:

\begin{equation*}
ME = \frac{e}{m} \sum_{p = 1}^{2^{m}} -\frac{n_{p}}{N} \log_{10}(frac{n_{p}}{N}),
\end{equation*}

where $m$ is the size of the m-tuple, $N$ is the number of reads covering the m-tuple and $n_{p}$ is the number of reads with the $p^{th}$ methylation pattern[^dmeas_problem]. $e$ is never properly defined in the paper but is described as the "entropy for code bit". In a supplementary file available on the `DMEAS` website, it is stated that $e = \frac{\ln(10)}{\ln(2)}$ (__CITE http://ufpr.dl.sourceforge.net/project/dmeas/DMEAS%20user%20guide.pdf__).

[^dmeas_problem]: The description of the `DMEAS` software in \citet{He:2013cj} allows for methylation loci with an unknown methylation state. Under this model there are now $3^{m}$ possible methylation patterns per m-tuple and the summatation index should be updated accordingly.

I think that this means it is the entropy value of the $p^{th}$ methylation pattern, where an unmethylated CpG is 0 and methylated CpG is $1$[^dmeas_problem2]. This amounts to $e$ being equal to the number of methylated CpGs in the m-tuple. For example, $e = 4$ if the pattern contains 4 methylated CpGs and $e = 2$ if the pattern contains 2 methylated CpGs. This is then normalised to account for the size of the m-tuple, $0 \leq \frac{e}{m} \leq 1$.

[^dmeas_problem2]: The `DMEAS` software further defines a CpG of unknown methylation state as 2. This seems to me like it would screw things up because $ME > 1$ is possible if an m-tuple contains a CpG with unknown methylation state.

ME is minimal ($ME = 0$) when all reads mapping to an m-tuple have the same methylation pattern and is maximal ($ME = 1$) when all $2^{m}$ methylation patterns are observed at equal frequency.

The statistical significance of an observed $ME$ is assessed using a permutation test. The null hypothesis is that the observed methylation patterns are purely stochastic. Specifically, for each m-tuple they simulate $10,000$ datasets with the same $N$ and average $\beta$-value as the observed data and compute the $ME$. The observed $ME$ is compared to this permutation distribution to compute a P-value. They also devise a test of "allele-specific methylation", although strictly speaking this is a test of whether there are two or more methylation subpopulations in the sample because it does not make use of genetic data.

As described  in \citet{He:2013cj}, `DMEAS` has very limited functionality because it can only examine 4-tuples of CpGs with $\geq 16$x sequencing coverage. Furthermore, `DMEAS` is only available as compiled software for Windows or as a Perl script that is provided in a `PDF` document.

## \citet{Li:2013gn}

[`eDMR`](https://github.com/ShengLi/edmr), published in \citet{Li:2013gn}, is software to identify DMRs from two-group experiments using bisulfite sequencing data. The publication focuses on CpG methylation from enrichment assays, such as RRBS, and it is unclear whether the software scales to WGBS data. The example in the dataset is an RRBS experiment of ten accute myeloid leukemia (AML) samples and five normal bone marrow (NBM) samples.

`DMR` is provided as an `R` package and is intended to be compatible with the `methylKit` R package. Indeed, it first uses `methylKit` to identify DMCs using Fisher's exact test at CpGs with at least $10\times$ coverage in at least 3 NBM samples. By using Fisher's exact test they are ignoring all within-group variability of methylation levels.

`eDMR` "aims to optimize the threshold for determining DNA methylation regions and to perform statistical significance testing". The first step is to define clusters of CpGs. They do this using this somewhat convoluted algorithm:

1. Compute the empirical distribution of $IPD$s for all CpGs with $\geq 10\times$ coverage[^edmr_problem1]
2. Fit a mixture of Normal distributions to this empirical distribution.
3. Determine the best separation point between the two components of the mixture distribution. This step incorporates a cost function, the role of which I don't really understand. This separation point, $D$, is used to define the boundaries of regions. $D = 183.5$ in the AML vs. NBM example.
4. Scan the genome for CpGs with $IPD > D$[^edmr_problem2]. The first CpGs of such a pair is the end of the previous region and the second CpG of such a pair is the start of the next region. The average level of methylation within each region in each group is the average of the $\beta$-values within that region in each group.
5. eDMRs are then defined as those regions with satisfying the following default criteria: contains at least 1 DMC, contains at least 3 CpGs, and the absolute mean methylation difference is greater than $20\%$. The exact cut-offs can be specified by the user.

[^edmr_problem1]: It is not clear for how many samples these sites must have $\geq 10\times$ coverage.
[^edmr_problem2]: It's not clear if all CpGs are used or just those with sequencing coverage $\geq 10 \times$.

Each eDMR is then assigned a P-value by combining the P-values from the Fisher's test used to identify DMCs in that eDMR. P-values are combined using the Stouffer-Liptak method (__CITE__). These final P-values are then adjusted using, I believe, the Benjamini-Hochberg FDR procedure.

## \citet{Lyko:2010dr}

__This should go in the co-methylation section.__

\citet{Lyko:2010dr} performed WGBS of honey bee (_Apis mellifera_) brains. It was known that young bees fed large amounts of royal jelly developed into queens (fertile females) whereas those fed smaller amounts develop into drones (males) or workers (infertile females). Furthermore, \citet{Kucharski:2008gu} had shown that a similar effect could be achieved by silencing expression of the DNA methyltransferase, DNMT3. Therefore,\citet{Lyko:2010dr} investigated the hypothesis that honey bees fed large amounts of royal jelly (queens) have different brain methylomes to those bees fed smaller amounts of royal jelly (workers).

As a part of their study, \citeauthor{Lyko:2010dr} investigated autocorrelation of $\beta$-values at CpGs in honey bees. It is not clear whether they only considered pairs with $NIL=0$ or whether they used all pairs ($NIL \geq 0$).

They found that

> As in the human and Arabidopsis genomes, methylation in Apis shows evidence of periodicity, although due to a much lower density of modified CpGs in this species the periodicity of 10 nucleotides (one helical DNA turn) is not obvious. However, a 3-base periodic pattern is clearly detectable, reflecting a preferential methylation of CpGs occupying the first and second position of the arginine codons.

However, the absolute value of these correlations are __tiny__, with the maximum being less than $0.015$. I suspect there is an error in these plots (Supplementary Figure 7, __DISCUSS WITH TERRY__).

They did not stratify their analysis by genomic elements, such as CGIs.

## \citet{Peng:2012dh}

\citet{Peng:2012dh} descibe an algorithm for identifying allele specific methylation (ASM). Their method assigns each read to one of $N_{h}$ "epigenomes" based on the pattern of methylation along the read. Each "epigenome" corresponds to what I call a haplotype, and the model assumes that the number of haplotypes in the sample is equal to the ploidy of the sample, e.g., $N_{h} = 2$ for a human sample.

Basically, \citeauthor{Peng:2012dh} assume a $H$-component mixture model and assign reads to each component, i.e. haplotype, using an EM algorithm. The output is a $N_{loci} \times N_{H}$ matrix, $A$, where $A[i, h]$ is the EM-estimate of methylation at the $i^{th}$ loci in the $h^{th}$ haplotype.

Unlike other methods for detecting ASM, this algorithm does not rely on methylation loci being close to heterozygous SNPs and therefore, at least, theoretically, can detect ASM in any region of the genome.

In practice, this method will only work well in regions with a high density of methylation loci and where the $N_{h}$ epigenomes are quite distinct from one another. The example described in the paper focuses on PMDs, regions where we know there are multiple methylation patterns, in _Arabidopsis_, an organism with a "reasonably high" density of methylation loci.
Furthermore, the example dataset is rather artificial because it is a synthetic dataset made by combining reads from two _Arabidopsis_ methylomes.

__There is no software available that implements the proposed algorithm.__

## \citet{Qu:2013ji}

`MLML`, published in \citet{Qu:2013ji}, is software to jointly estimate 5mC and 5hmC levels from samples sequenced with any two of BS-seq, oxBS-seq and TAB-seq. Depending on the exact combination of technologies, the estimation of either 5mC or 5hmC levels is based on a "subtraction" of read counts between the two experiments. A naive "subtraction" can result in estimates that are $< 0$ or $> 1$.

`MLML` implements an EM algorithm to estimate both 5mC and 5hmC levels under a Binomial-mixture likelihood. The resulting estimates are guaranteed to be "non-negative, and never sum over one". __It's not immediately clear why `MLML` does not guarantee that the estimates sum to exactly one.__

## \citet{TricheJr:2013tj}

\citet{TricheJr:2013tj} investigated different regression models for identifying DMCs and DMRs from Illumina methylation beadchip data. They report that a Beta-regression model[^beta_regression] "yielded more validated hits" than models that used Wilcoxon-Mann-Whitney tests, Student's t-tests or Welch's t-tests.

[^beta_regression]: A "Beta regression model" means explicitly modelling the $\beta$-values as Beta random variables. This is distinct from a "regression of the $\beta$-values, where other distributional assumptions or transformations of the $\beta$-values might be used in a standard linear model framework.

Beta regression (__CITE Ferrari and Cribari-Neto__) does not fit into the generalise linear model framework (__CHECK WITH TERRY__) because the mean and variance terms are coupled to both depend on covariates.  A software implementation of the proposed Beta regression model is __AVAILABLE WHERE__.

\citeauthor{TricheJr:2013tj} caution that when the number of samples per group is small ($< 25$) or when the design is unbalanced that control of the Type-I error rate may be lost.

## \citet{Xu:2013eg}

\citet{Xu:2013eg} propose a test to identify DMCs from bisulfite-sequencing data in two-group experiments. Their model is based on a method for the analysis of clustered binary data, published in __CITE Rao and Scott, 1992__.

For the $i^{th}$ individual, the methylation level at the $j^{th}$ locus is modelled as a binomial random variable, $M_{i, j} = Bin(M_{i, j} + M_{i, j}, B_{i, j}$. Within each group, the observed $\beta$-values are adjusted by a group-specific "design effect". The design effect for the $i^{th}$ locus in the $l^{th}$ group, $d_{i, l}$, is the ratio of the estimated $\beta$-value variability while accounting for within-group variability to without accounting for within-group variability (__I THINK__). The design effect is used to correct both the average withing-group count of methylated reads and the average within-group sequencing coverage.

The test statistic for the $i^{th}$ locus is then a $\chi^2$ statistic on 1 degree of freedom comparing the each group's design-adjusted average $\beta$-value to the design-adjusted grand average $\beta$-value. __This model seems convoluted and I wonder whether it fits into a simpler GLM analaysis of deviance framework__.

__There is no software available that implements the proposed algorithm.__

### Simulation study {-}

\citet{Xu:2013eg} includes a simple simulation study. The "true" group-wise methylation levels at each methylation loci, $B_{i, l}$, are simulated from a beta, a (truncated) Normal or a (truncated) Normal-mixture distribution. Thus all samples within each group have the same "true" methylation level and the methylation level at each locus is independent. Sequencing coverage for each individuals was sampled from a (rounded) Normal distribution.

## \citet{Zhang:2012id}

\citet{Zhang:2012id} propose a Dirichlet process beta mixture model (DPBMM) to cluster methylation profiles. The Dirichlet process effectively means there is an 'infinite' number of components to the mixture distribution. The resulting model is analytically intractable and so they propose a Gibbs sampler to compute the relevant posterior distributions. DPBMM does not require that the number of clusters be pre-specified as this is estimated as part of the Gibbs sampler.

MATLAB code implementing the algorithm is available from [https://sites.google.com/ site/bdpmmmethy/home](https://sites.google.com/ site/bdpmmmethy/home).

## \citet{Zhang:2011dp}

`QDMR`, published in \citet{Zhang:2011dp}, is another method for identifying DMRs based on a modified form of Shannon Entropy. `QDMR` uses $\beta$-values and is so applicable to both microarray- and sequencing-based bisulfite-treatment assays. The publication uses MeDIP-chip and RRBS data.

The methylation entropy of a locus is computed using all samples, without regard to an 'group' variable. The methylation entropy of a region, rather than a single methylation locus, can be computed by using the sample-wise 'averages' of $\beta$-value across the region.

The `QDMR` software is available from [http://bioinfo.hrbmu.edu.cn/qdmr](http://bioinfo.hrbmu.edu.cn/qdmr).

## \citet{Landan:2012kp}

\citet{Landan:2012kp} investigated the stochastic variability of methylation patterns at CpG 4-tuples. Borrowing concepts from population genetics, they call a group of adjacent methylation loci an _allele_ and each possible methylation patterns at the allele an _epiallele_. In my framework, an allele is an m-tuple and an epiallele is a (section) of a haplotype.

\citeauthor{Landan:2012kp} defined the _epipolymorphism_ of an allele as "the probability that two epialleles randomly sampled from the [allele] differ from each other". Mathematically, the epipolymorphism of a 4-tuple is defined as $p_{poly} = 1 - \sum_{l = 1}^{l = 2^{4}} p_{l}$, where $p_{l}$ is the relative frequency of the $l^{th}$ epiallele in the sample. They also defined the average methylation level of a 4-tuple, $\zeta_{(i, i + 1, i + 2, i + 3)}$, as follows:

1. Compute the methylation level of each read overlapping that 4-tuple = $0, 0.25, 0.5, 0.75, 1$.
2. Compute the average methylation content of all reads overlapping that 4-tuple, $\zeta_{(i, i + 1, i + 2, i + 3)}$. $\zeta_{(i, i + 1, i + 2, i + 3)}$ is equivalent to the average $\beta$-value in that m-tuple, $\frac{1}{4} \sum_{i' = i}^{i' = i + 3} \beta_{i'}$, __provided that all reads overlapping any of the methylation loci in that m-tuple overlap all methylation loci in that m-tuple.__

They then plotted the level of epipolymorphism against the average level of methylation for all 4-tuples (with "suffient" coverage) and compared the resulting curves to two theoretical model of methylation variability. These two models were _bimodal epipolymorphism_ and _maximum epipolymorphism_. The bimodal epipolymorphism model assumes that each read, $z$, has $\zeta_{z} = 0$ or $\zeta{z} = 1$, that is, each read is either fully unmethylated or fully methylated. Mathematically, 'bimodal epipolymorphism' can be written as $2m (1 - m)$, where $m$ is the average methylation level of the 4-tuple. The maximum epipolymorphism model was modelled as $1 - (1 - 2(m(1 - m)))^4$ (__TODO: Figure out what the maximum model corresponds to__).

In the ENCODE data they found that differentiated, somatic tissues had a higher frequency of epipolymorphism than did embryonic cells. From this they concluded that germline and pluropotent cells are able to establish or maintain a more coherent methylation state and that somatic cells accumulate substantial stochastic variation during somatic development.

\citeauthor{Landan:2012kp} also found that DMRs identified between the H1 embryonic stem cell line and somatic cell lines had a higher degree of epipolymorphism than non-DMRs[^landan_dmrs]. They found a similar pattern in the bisulfite tumour vs. normal padlock probe bisulfite sequencing data from \cite{Hansen:2011gu}.

[^landan_dmrs]: Sample-specific DMRs were simply defined as regions with $\geq 20\%$ difference in average methylation between the (duplicate-averaged) H1 methylome and the sample of interest. No details are provided on how large the regions had to be to qualify as DMRs.

\citeauthor{Landan:2012kp} then studied the evolution of DNA methylation patterns in an _in vitro_ model system tracking immortalised fibroblasts over 300 generation Two lines, A and B, were tracked and sampled at several time points for profiling with MeDIP-seq. Unfortunately, at least from my perspective, MeDIP-seq does not provide single-base-resolution of DNA methylation and so I do not discuss these result any further.

They did, however, also profile 45 cancer-related CGIs in the the model system at multiple time points with "ultra-deep" ($> 10,000 \times$ coverage) bisulfite-sequencing. By sampling the over time they were able to investigate the evolution of stochastic epipolymorphism. Using this data they also analysed the correlation of DNA methylation at pairs of CpGs in these 45 regions. They visualised these using plots analogous to linkdage-disequilibrium plots. They did not look at how these correlations vary as a function of distance.

From the analysis of this small set of regions they concluded that:

> [the] correlation between the methylation states of pairs of CpGs was generally very low. This lack of correlation suggests that methylation dynamics are typically independent for different CpGs, making the methylation state of one CpG (whether high or low) uninformative on the methylation state of nearby CpGs.

__My co-methylation results and extensive $\beta$-value correlations suggest otherwise.__

The correlation of methylation at a pair of CpGs was computed from the $2 \times 2$ contigency table, shown __BELOW__.

__TODO: Represent as contingency table.__ Let $n_{11}$ be the number of reads that are methylated at both CpGs, $n_{00}$ be the number of reads that are unmethylated at both CpGs and $n_{01}$ and $n_{10}$ be the number of reads that have different methylation states at each CpG ($N = n_{00} + n{01} + n{10} + n_{11}$). The average methylation level of each CpG is $m_{1} = \frac{n_{11} + n_{10}}{N}$ and $m_{2} = \frac{n_{11} + n_{01}}{N}$, respectively. The variance of methylation at each CpG is $v_{1} = m_{1} \times (1 - m_{1})$ and $v_{2} = m_{2} \times (1 - m_{2})$, respectively. The correlation of methylation at the pair of CpGs is then $\frac{n_{11} / N - m_{1} \times m_{2}}{\sqrt(v_{1} + v_{2})}$. __TODO: Is this the same as correlation of $beta$-values using only those reads that overlap both CpGs?__

## \citet{Zhang:2013uu} and \citet{Stevens:2013hv}

`M&M`, published in \citet{Zhang:2013uu}, is a statistical framework for jointly analysing MeDIP-seq and MRE-seq to detect DMRs in a two-group experiment. MeDIP-seq is an enrichment-assay and MRE-seq is a restriction enzyme assay. `M&M` is a window-based algorithm (default = 500 bp non-overlapping windows). `M&M` was compared to `MEDIP`, which only uses MeDIP-seq data, and WGBS of the same samples. Each method was run on biological duplicates to estimate the corresponding false positive rate.

`methylCRF`, published in \citet{Stevens:2013hv}, is a software implementation that extends the `M&M` framework. Specifically, it uses a conditional random field model of counts-per-window whereas `M&M` assumes independence across windows. It gives single-base estimates for DNA methylation, despite the input being data from enrichment-assays. \citeauthor{Stevens:2013hv} provide evidence that these predicts single-base estimates are very close to observed values from WGBS of the same sample.

\citet{Stevens:2013hv} plot $\beta_{i}$ vs. $\beta_{i+1}$ for $IPD \leq 750$ bp to argue that the distribution of methylation is highly non-random (Fig 2b).

## \citet{Hon:2013jra}

### Summary of method for identifying tsDMRs

The authors aim to identify tissue-specific differentially methylated regions (tsDMRs). They study 17 different tissues from mouse, although they exclude the placental methylome from most analyses, which reduces the total number of tissues to 16.

To identify tsDMRs they use a Hidden Markov Model (HMM) on the $\chi^2$ test statistics of different methylation levels across the 16 tissues. Specifically, the $\chi^2$ statistic was calculated as follows:

> For a given genomic interval, let $m_{t}$ be the number of methylated cytosines sequenced and $d_{t}$ be the depth of sequencing for tissue $t$. To capture the deviation of the methylomes from that expected if all tissues were uniformly methylated, we used a $\chi^2$ test statistic. We denoted the uniform methylation level in this interval as $f = \sum m_{t} / \sum d_{t}$. Then, the expected number of methylated cytosines sequenced for each tissue was $e_{t} = f d_{t}$. Thus, the $\chi^2$ test statistic was calculated by $\chi^2 = \sum(m_{t} - e_{t})^{2} / e_{t}$, with the degrees of freedom equal to one less than the total number of tissues.

It is not clear that they correctly adjusted the degrees of freedom when the number of available tissues (i.e. those with sufficient sequencing coverage) was less than 16 for a given region.

> To identify tsDMRs, we calculated $\chi^2$ values for each CpG unit, __defined as three consecutive CpGs__

The genomic intervals used in their analyses are groups of 3 CpGs, which they at times refer to as _CpG units_. I assume that these CpG units do not overlap but this isn't made clear. This $\chi^2$ test statistic differs from that used in my analysis of the _Avy_ methylome data because I test at every CpG whereas they test at every group of 3 CpGs.

Then, a Hidden Markov Model (HMM) was used to identify putative tsDMRs:

> ..., [we] employed an HMM as implemented by [`pmtk3`](https://github.com/probml/pmtk3) to segment the genome, using methods as previously described[^58] with several alterations. Briefly, we trained a four-state HMM, with each state consisting of a mixture of two Gaussians, with the Baun-Welch [sic] algorithm. States were estimated using the forward-backward algorithm, with the highest-valued $\chi^2$ state denoting a prefiltered set of tsDMRs[^wot].

[^58]: Dixon, J.R. et al. Topological domains in mammalian genomes identified by analysis of chromatin interactions. Nature 485, 376–380 (2012).

[^wot]: "with the highest-valued $\chi^2$ state denoting a prefiltered set of tsDMRs"; what does this even mean?

__There is not enough detail to reproduce the HMM used in their analysis__. They don't describe what the four states in the HMM actually are. Nor do they describe the "several alterations" made to the method in Dixon et al. Moreover, based on a quick scan, Dixon et al. only contains a generic description of a 3-state HMM, which isn't very helpful in trying to understand the methods in the present paper. In summary, the HMM isn't described in sufficient detail to reproduce it.

The putative tsDMRs were then filtered according to variation across tissues:

> To select for tsDMRs with a large magnitude of variance in tissue methylation, we determined the intersection point of tissue-specific and non–tissue specific standard deviation distributions (Fig. 2e, right) and removed sites with smaller standard deviation than this intersection point (here $11.88\%$ mCG). The final set of tsDMRs was represented by three genomic loci: the left and right CpG boundaries were called by the HMM, along with the central CpG site having the highest $\chi^2$ value, which represented an estimate of the methylated base showing the greatest tissue specificity.

No particular comments on this step.


They report an FDR for each region:

> To estimate the false discovery rate, we repeated this analysis on ten random permutations of the data set, each of which consisted of random assignment of methylated base calls within each tissue while maintaining baselevel sequencing depth.

They claim that the intial HMM segmentation identifies DMRs at the unbelievable (in the literal sense of the word) FDR value of $3.96 \times 10^{-5}$. My instinct is that the number of permutations is insufficient to accurately estimate the FDR and I take the reported FDRs with a large handful of salt.

Given a putative tsDMR, the tissues driving that variability eredentified as those with the highest specificity relative to the Shannon entropy of the region. This basically amounts to checking which of the tissues had the "most different" methylation levels within putative tsDMRs.

### Other questions

#### Are these tissues from the same mouse? Do they have multiple mice per tissue? How much of the variation is mouse-specific (sequence-driven) rather than tissue-driven? {-}

The Methods says that, "Mouse tissues were isolated from a female C57Bl/6 mouse at least 4 weeks old (Charles River) at 14.5 d of pregnancy". I _think_ this means a single mouse, although the wording is somewhat ambiguous. If multiple mice were used it does not say how many mice were used nor which tissues were collected from which mouse. Later on they say that, "Previously published mouse methylomes (ES and NPC) were downloaded and remapped". So, at least two of the tissues were from different mice.

While they did use methylomes from Cast/129 and 129/Cast hybrids to replicate some of the putative "cortex vestigal enhancers" identified in the C57Bl/6 methylomes, if the majority of the tissues are $n=1$ from a single mouse then there is no sense of the variability of these tsDMRs across mice.

### How could this be adapted to the _Avy_ methylome data?

The basic idea would be to run a HMM on the $\chi^2$ statistics. There are two key things to consider:

1. What goes into the HMM, i.e. what statistic?
2. What sort of HMM?

#### What goes into the HMM? {-}

You could do this on CpG-specific $\chi^2$ statistics (what we have now) or on $\chi^2$ statistics computed from groups of CpGs (e.g. groups of 3 as done in the present paper). There's no reason to believe that groups of 3 is the "correct" number to use if grouping, and I think the amount of grouping would depend on the region's characteristics (e.g. CpG density, spatial correlation of methylation of the region).

They suggest grouping because they have low coverage (~$8\times$) sequencing data, whereas the _Avy_ methylome data has a higher sequencing coverage, so in that sense grouping may be unnecessary for the _Avy_ methylome data.

#### What sort of HMM? {-}

As discussed above, the HMM used in this paper is not described in sufficient detail to replicate it. In particular, I don't know what the four states are. But if I was to try to use a HMM for this problem in the _Avy_ methylome data I would start with two states: "no difference in methylation across the 5 samples" vs. "difference in methylation across the 5 samples".

There still remain several other important choices to parameterise the HMM such as the transition probabilities and emission distributions, which again aren't described in sufficient detail in the present paper.

In summary, implementing a HMM would require a fair bit of exploratory work to identify a good parameterisation.

## \citet{Knijnenburg:2014ke}

\citeauthor{Knijnenburg:2014ke} propose a segmentation algorithm, `MSR` (multiscale signal representation), to represent a vast array of genomic signals that are measured on different scales. Such signals include DNA methylation.

An example of segmentation DNA methylation data, comparing a colon tumor to its matched normal mucosa, uses the log-ratio of tumor to normal DNA methylation levels. While this is a "differential methylation" signal, as the authors' note, `MSR` does not explicitly test for differential methylation as its purpose is genome segmentation:

> The `MSR` method is complementary to other segmentation techniques, such as `Segway` and `ChromHMM`: whereas the `MSR` segments one genomic signal at multiple scales, these HMM-based models use multiple genomic signals to provide one segmentation, thereby dividing the genome into a number of functional states called chromatin domains.

This example, like all others in the paper, is effectively an $n=1$ analysis. __It is not clear how this method would be generalised to include biological variation.__

While the authors propose a method to test for significant signal, it is not clear whether or how this method accounts for the fact that the same dataset is used to discover and test the features. This score, called _SFC_ (significant fold change), seems similar in spirit to the CDIF proposed in \cite{Sun:2014fk}.

The code is written in Matlab and available from [https://github.com/tknijnen/msr/](https://github.com/tknijnen/msr/). It includes a runtime environment that allows users who don't have Matlab installed to run the `MSR` software.

## \citet{Hovestadt:2014fm}

\citeauthor{Hovestadt:2014fm} performed WGBS on 34 human and five murine medulloblastoma tumors, and eight human and three murine normals. They are interested in integrating DNA methylation with other measures, such as RNA expression and histone modifications. Much of the analysis of the WGBS data is "novel", and is based on binning and thresholding. To identify sub-group-specific differential methylation they do the following:

> For all patterns described earlier, the average methylation level per gene and sample was calculated for the entire WGBS cohort. Analysis of variance was used to determine genes differentially methylated between all four medulloblastoma subgroups and combined controls (adjusted P value , 0.001). The Benjamini–Hochberg procedure was used to adjust for multiple testing. Subgroup-specific differentially methylated genes were determined by applying a post-hoc test on genes previously determined as being differentially methylated. (R package: `multcomp`, functions: `mcp` and `glht`). Individual medulloblastoma subgroups were required to be significant against all other subgroups combined, and against control samples separately (P value , 0.001; in essence eight comparisons were made). This testing procedure was deemed appropriate because of large sample sizes and continuous methylation levels resulting of averaging of multiple CpGs.

## \citet{Ushijima:2003gqa}

This paper uses a clever strategy of starting with a single cell, and then growing clonal populations to track the fidelity of DNA methylation. As it is an older paper, they used Sanger bisulfite-sequencing, which means they only studied a small number of regions (albeit at base pair resolution).

> We determined the methylation status of each CpG site on each DNA molecule obtained from clonal populations of normal human mammary epithelial cells. Methylation pattern error rates (MPERs) were calculated based upon the deviation from the methylation patterns that should be obtained if the cells had 100% fidelity in replicating the methylation pattern. Unmethylated CGIs in the promoter regions of five genes showed MPERs of 0.018–0.032 errors/site/21.6 generations, and the fidelity of methylation pattern was calculated as 99.85%–99.92%/site/generation. In contrast, unmethylated CGIs outside the promoter regions showed MPERs more than twice as high (P < 0.01). Methylated regions, including a CGI in the MAGE-A3 promoter and DMR of the H19 gene, showed much lower MPERs than unmethylated CGIs. These showed that errors in methylation pattern were mainly due to de novo methylations in unmethylated regions.

For some of the regions, they identified two major clonal sub-populations, which they interpreted as being the initial methylation pattern of each allele in the original single cell.

They also show evidence that in their experiment the the fidelity errors occur at least an order of magnitude more frequently than failures of the bisulfite conversion process:

> These [bisulfite conversion efficiency] values showed that the MPERs in CGIs in the promoter regions are 10-fold more than the unconversion rates.

In theory, could estimate epiallele frequencies from Figure 3. The analysis methods are very simple proportions of 'CpGs that changed methylation state'.

## \citet{Capra:2014uh}

This is the first paper I've seen that uses phylogenetic approaches to studying DNA methylation dynamics. Specifically, they seek to infer DNA methylation changes from RRBS data of the haematopoetic lineage. The main result of interest is that they try to incorporate the "vertical" correlation of DNA methylation dynamics between precursor and dependent instead of the "horizontal" correlation of neighbouring CpGs within the same sample. The vertical approach beats the (simple) horizontal approach in inferring missing data.

This approach discretisation of $\beta$-values into lowly-, intermediate- and highly-methylated categories.

Putting DNA methylation into a phylogenetic framework gives (in principle) access to a whole raft of existing analysis methodologies.

## \citet{Lee:2014cv}

Figure 3 describes several scenarios that would lead to DNA methylation heterogeneity.

## \citet{Gokhman:2014cp}

This paper compares DNA methylomes from a Denisovan sample and a Neandertal to methylomes from present day humans. They exploit the natural deamination of cytosines (mC -> thymine and C -> uracil) to infer the methylation state of cytosines from the C->T conversion ratios. __IMPORTANT:__ These C->T ratios are __not__ the same as $\beta$-values, but they can be scaled by a number, $\frac{1}{p_{m}}$, so that they have a similar interpretation.


They smooth these C->T ratios using a sliding window, where the window size is chosen as that which maximises the correlation with methylomes from present day humans (43-89 CpGs, depending on Denisovan or Neandartal and autosome or X-chromosome).

The authors advocate smoothing of DNA methylation:

> Such noise reduction is advisable in methylation samples, as smoothing over neighboring CpGs increases accuracy for low-coverage positions and complies with the high correlation in methylation between adjacent positions and the mostly regional function of DNA methylation

and cite the `BSmooth` paper in support of this position.

They use a rather weird way of removing potential PCR duplicates:

> To remove PCR duplicates that are characterized by abnormally high coverage, we fitted the (strand-specific) coverage distribution to a mixture of (two) Gaussians. Positions with a coverage that was higher than three standard deviations from the mean of the main distribution were removed from the analyses. In the Neandertal, the threshold was set to 46 reads and in the Denisovan the threshold was set to 30 reads.

I don't really understand why they used this and not the standard Picard `MarkDuplicates`-type approach; perhaps it's because the sequencing protocol removes uracils, which would lead to variation in read lengths?

### Identification of DMRs

There's a convoluted procedure to estimate what amounts to the average $\beta$-value of each window, and its associated standard deviation, which are used to form a z-stat. P-values are FDR-adjusted and only windows with FDR $< 0.05$ are retained. Furthermore, only windows with a greater than $60\%$ change in methylation are retained for further analysis.

## \citet{Rijlaarsdam:2014hp}

`DMRforPairs` tries to find DMRs in 1 vs. 1 comparisons - an ambitious aim. I understand the desire for such comparisons but I wouldn't believe the results of most $n = 1$ studies.

## \citep{Chen:2014jb}

This paper proposes a non-parametric test (_Cuzik_ test) of differences in DNA methylation __between multiple groups__, that also accounts for the well-known association of age and DNA methylation. The model requires that the age variables is discretised into bins. There is no discussion of how to choose bins (they use 5-year bins in the real data analysis). Furthermore, my intuition is that leaving age as a continuous variable is preferable.

The Cuzik test assumes there is a trend of differential methylation with age, either increasing ($W_{1}$) or decreasing ($W_{2}$). As the trend is unknown, the authors propose using $W = max(W_{1}, W_{2})$ and then make use of an asymptotic bound to obtain approximate P-values. The description of the Cuzik test sounds rather different [on Wikipedia](http://en.wikipedia.org/wiki/Cuzick%E2%80%93Edwards_test) - I'm not sure what to make of that.

The authors make the reasonable observation that if there is a trend (even if it's not truly linear) then there test, which is designed to detect such a trend, has more power than do tests that aren't designed to detect a trend.

### Simulation study

$\beta$-values are simulated from marginal distributions (Uniform, truncated Normal or Beta), which do not include the spatial correlation of $\beta$-values. They specifically simulate scenarios where there is a trend-with-age of the $\beta$-values, thus ensuring there test looks good, and do not consider a setting where there is no trend.

## \citep{Park:2014ho}

This paper describes the `methylSig` software. `methylSig` uses a Beta-Binomial model to test for DMCs across an arbitrary number of groups (although the proposed test is a likelihood ratio test of group $k$ vs. group $k'$). The dispersion parameter in the Beta-Binomial model, $\theta_{i}$ (assumed constant across groups), is treated as a nuisance parameter that is first estimated and then plugged into the likelihoods to estimate (using maximum likelihood) the parameters of interest, $\mu_{i, k}$ and $\mu_{i, k'}$. For small sample sizes, they use the $t$-distribution, rather than the 'natural' $\chi^{2}_{1}$-distribution of a likelihood ratio statistic. The dispersion parameters are not shrunk, e.g. via Empirical Bayes.

They also propose incorporating local information, via kernel smoothing, to estimate the parameters $\mu_{i}$ and $\theta_{i}$.

The examples all use ERRBS (enhanced RRBS) data and so it isn't clear whether `methylSig` scales to WGBS data.

To increase power, the authors suggest "tiling" observations. Tiling amounts to pooling data, i.e. summing all measurements within a tile (25 bp for ERRBS data, larger for WGBS)) and then using the tiles as the observations rather than the CpGs.

The plots of the results are very busy and I don't like them much.

### Identifying differentially methylated transcription factors

Rather than looking for general DMRs, the authors focus on identifying differentially methylated transcription factors. There are two approaches:

1. Test which TFs contain a significant number of DMCs.
2. Test which TFs are hypo- or hyper-methylated (combining all CpGs in the TF).

### Simulation study

The simulation study is based on ERRBS data from a study of acute myeloid leukemia. It uses the coverage distribution, the locations of CpGs and the estimated dispersions from the real data. However, the $\beta$-values are simulated from Beta-Binomial distributions with group-specific means to induce DMCs. This model does not incorporate correlation amongst the $\beta$-values, although the simulation models allows for DMCs to cluster into DMRs, which is similar to inducing correlations amongst the $\beta$-values, albeit restricted to CpGs in DMRs.

The authors advocate for using `methylSig` with dispersions estimated locally via kernel smoothing. They note that `BSmooth` has low power to detect "independent DMCs", but this is not surprising because `BSmooth` is designed to detect DMRs.

## \citet{Smallwood:2014kn}

This is the first paper that I know of with single-cell bisulfite sequencing data (there is also a scRRBS protocol). The technique, _scBS-seq_, is a variation on PBAT-seq. The single-cell DNA underwent five cycles of random priming, extension and purifying, which means there were up to five 'copies' of each allele.

It is __not a genome-wide technique__, since it only gives "accurate" results for $48.4\%$ of CpGs.

They profile ovulated metaphase II oocytes (MII) because:

> MIIs are an excellent model for technical assessment as they: (i) can be individually handpicked to ensure that only one cell is processed; (ii) represent a highly homogeneous population, which allows discrimination between technical and biological variability; and (iii) present a distinct DNA methylome comprising large-scale hypermethylated and hypomethylated domains.

They contrast the MIIs with embryonic stem cells (ESCs):

> ESCs grown in serum conditions exist in a state of dynamic equilibrium characterized by transcriptional heterogeneity and stochastic switching of transcriptional states9–12, and emerging evidence from immunofluorescence and locus-specific studies suggests that 5mC heterogeneity exists in ESCs13.

A fascinating fact from the supplementary material, which I don't understand:

> The overall higher mapping efficiency of oocytes versus ESCs can be explained by the amount of DNA in each cells (4n for MII oocytes and 2n for ESCs), resulting in a relatively lower contribution of spurious sequences in MIIs.

They sequenced multiple single cells (51 across from four different 'populations') at very shallow depth ($\approx 20 \times 10^{6}$ 100bp PE reads/sample). Furthermore, mapping efficiences were very low ($\approx 3.9 \times 10^{6}$ reads/sample or $20\%$), which is largely due to low-complexity sequences (specifically, poly-Ts that are inherant to the protocol). They also sequenced "bull cell mass" samples, i.e. pools of the same cells. They could largely recapitulate the $\beta$-values of the pools using 12 single-cells.

## \citet{Robinson:2014tx}

Mark has written a review of methods to detect DMRs. It is fairly high-level (as it should be) and could basically act as my description of calling DMCs and DMRs! For example:

> It is therefore no surprise that BB assumptions are made in several recently proposed packages, such as BiSeq [30], MOABS [29], DSS [28] and RADMeth [33]. Similarly, empirical Bayes (EB) methods fit naturally for modeling and inference across many types of genomic data and DNA methylation assays are no different. MOABS and DSS both implement hierarchical models and use the whole dataset to estimate the hyperparameters of the beta distribution; RADMeth and BiSeq use standard maximum likelihood without any moderation. While BiSeq and RADMeth do not moderate parameter estimates, they provide facilities for complex designs through design matrices, which MOABS and DSS do not currently offer.

Mark notes that neither Fisher's exact test, nor logistic regression, control for within-group variability:

> While this strategy (Fisher's exact test) may be sufficient in comparing cell lines, we stress that the use of FET should be avoided in the general case; most systems have inherent biological variation and FET does not account for it.

> Likewise, using the binomial distribution, such as a logistic regression framework (e.g., methylKit; 32), also does not facilitate estimation of biological variability, unless an overdispersion term is used.

This quote fits with my intuition:

> Ultimately, our intuition suggests that moderated t/F-statistics on the normalized log-ratios of intensities seems most rigorous. -- Highlighted 25/07/2014

Mark emphasies the fact that there is a __big__ difference between methods that work on pre-defined regions and methods that seek to identify such regions:

> ... one must distinguish between methods that operate on predefined regions, with those that define regions of DM. The latter is considerably more difficult because ensuring control of the false discovery rate (FDR) at the region-level is non-trivial; in particular, controlling false discoveries at the site-level does not give a direct way to controlling false discoveries at the region-level when the region itself is also to be defined.

## \citet{Liu:2012ge}

I believe `Bis-SNP` is the best methylation calling software due to it's comprehensiveness. Unfortunately, it can only call methylation at 1-tuples and I need to be able to call methylation at m-tuples. `Bis-SNP` is built on `GATK (v1)` and so includes base-quality recalibration and indel realignment. Importantly, `Bis-SNP` works with any `BAM` file of bisulfite-sequencing data, which means that the user can choose his/her favourite alignment software.

As the authors state:

> Bis-SNP is a practical tool that can both (1) improve DNA methylation calling accuracy by detecting SNPs at cytosines and adjacent positions, and (2) identify heterozygous SNPs that can be used to investigate mono-allelic DNA methylation and polymorphisms in cis-regulatory sequences.

Figure 1 explains how it performs SNP calling from bisulfite-sequencing data.

Calling "C-strand" SNPs depends on the underlying methylation state ($\beta$, which defaults to the standard $\beta$-values) and bisulfite-conversion errors (underconversion = $\alpha$, which defaults to $0.25$) and overconversion = $\gamma$, which defaults to $0$).

Rather than filter out positions from an M-bias plot, `Bis-SNP` "walks" from the 5' end of each read and excludes all positions prior to the first $T$ that is mapped to a reference $C$. This is to remove the so-called "5' bias" of the Illumina bisulfite-sequencing but will do nothing for 3' bias.

## \citet{Dolzhenko:2014bo}

> Existing methods based on Fisher’s Exact Test and HMMs are appropriate for comparing a pair of samples at a time (coming either directly from the experiment or obtained by pooling other samples); however, they lack the ability to account for variability of methylation levels between replicates.

> Unlike BSmooth, [BiSeq] can be used for experiments that go beyond comparing two groups of samples, but it requires a set of candidate regions that may exhibit differential methylation.

> The beta- binomial, which has first been used for modeling WGBS proportions by Molaro and others [17], is a natural choice for describing methylation levels of an individual site across replicates as it can account for both sampling and epigenetic variability.

> The beta-binomial regression is fit separately for each target site...[parameters] are estimated using the method of maximum likelihood

> [To find DMRs] the p-values are transformed using weighted Z test (also known as Stouffer-Liptak test), employing an approach proposed by Kechris and others [23].

Correlations of z-scores are computed using the method implemented in `comb-P`.

> We recommend calculating correlation and subsequently combining the pvalues of sites located within 200 bp. In our experience correlation typically becomes much weaker beyond this point, so the users do not generally need to alter this parameter. However, it may be appropriate to increase the value of this parameter when analyzing very noisy data.

Sounds fair.

> As explained by Rice [33], Fisher’s method is best suited for testing the existence of at least one significant test among the ones being combined, while the Z test is more appropriate in situations requiring the consensus among all of the combined tests, suggesting that the Z test is more appropriate for our purposes.

### Some discussion of Rice

> In a wide variety of biological applications there is a need to combine the results from independent tests for which the raw data cannot be pooled.

True, but the tests in BS-seq are __not independent__. So, what are the consequences of this?

> Fisher's testing procedure represents a test against broad alternatives. It specifically tests whether at least one component test is significant, and can yield a significant combined test statistic when the component tests, on balance, strongly support $H_0$. This is an undesirable characteristic when asking whether a group of tests collectively supports the same $H_0$.

### Simulation

6 vs. 6 two-group experiment.

1. Sample coverage from real data
2. Sample number of methylated reads from a binomial distribution
3. Simulate methylation levels from distributions in [Epigenome-wide association studies for common human diseases](http://www.nature.com/nrg/journal/v12/n8/full/nrg3000.html). These distributions are hypothetical.

__No correlations of methylation levels in simulation.__

## \citet{Xie:2014ez}

A Bayesian test to replace the "Lister method" of identifying methylcytosines. I still don't think identifying methylcytosines is particularly useful.

### Simulation

They use the [`BSsim`](http://122.228.158.106/BSSim) software for simulation.

## \citet{Huh:2014ki}

`Bis-class` is a Bayesian method for clasifying sites as methylated or unmethylated. The authors claim it is most useful for lowly methylated genomes with low coverage, which is where simpler methods struggle.

The classifier incorporates information about the global and local methylation levels in the sample via a prior distribution. The "localness" is incorporated via a kernel with the window size set to be the distance where the spatial correlation between methylation levels drops below 0.2. __TODO: Do they describe how they compute these correlations__

## \citet{Chodavarapu:2010iqa}

> We identified 10-base periodicities in the DNA methylation status of nucleosome-bound DNA (Arabidopsis) and found that nucleosomal DNA was more highly methylated than flanking DNA. These results indicate that nucleosome positioning influences DNA methylation patterning throughout the genome and that DNA methyltransferases preferentially target nucleosome-bound DNA. We also observed similar trends in human nucleosomal DNA, indicating that the relationships between nucleosomes and DNA methyltransferases are conserved.

> Similar to findings in animal and fungal high-throughput nucleosome sequencing studies5–8, we found 10-base periodicities in WW (W 5 A or T) dinucleotides, and SS (S 5 G or C) dinucleotides that were 5 bases out of phase with the WW dinucleotides (Fig. 1d and Supplementary Figs 3–5).

> Notably, all three types of methylation showed a 10-base periodicity on nucleosomal DNA, which was in phase with the WW dinucleotides and out of phase with the SS dinucleotides (Figs 1d, 2a–c and Supplementary Fig. 6). These methylation preferences were not correlated with preferences for CG, CHG or CHH sequences at these locations (Supplementary Fig. 7). Because DNA methyltransferases access the major groove, this methylation would be on DNA that is on the outside of the nucleosome (minor groove facing the histones) and thus more accessible to the DNA methyltransferases. This in turn indicates that DNA might be in part methylated as the DNA is still bound to nucleosomes, leading to the observed 10-base pair periodicity.

> We previously reported a 10-nucleotide periodicity in CHH methylation data when performing autocorrelation analysis on the methylation pattern of the whole genome3. Our previous interpretation was that the structure of the DRM2 enzyme might be responsible for this pattern, because the orthologous Dnmt3 enzymes in mammals are known to act as heteromeric complexes in which two methyltransferase active sites have a spacing equivalent to roughly 10 nucleotides of DNA13. However, the current data in which we see 10-base pair periodicities for all types of methylation indicate a more general explanation: that nucleosomes are to some extent dictating access to the DNA and therefore setting the register of methylation for all DNA methyltransferases. Nucleosomal preferences could also partially explain the sequence preferences that we observed previously for CHG and CHH methylation3. Highly methylated cytosines tended to be followed immediately by A/T but not C, consistent with our finding that DNA methylation is out of phase with CC dinucleotides.

> This suggests the possibility that DNA methylation, which frequently exists in the transcribed regions of active genes26,28, could have a conserved role in exon definition or splicing regulation.

## \cite{Law:2010jr}

> the symmetric CG and CHG contexts (in which H = A, T or C) and the asymmetric CHH context

__This is not correct, CHG is not symmetric__. For example, a CCG on the + strand is a CGG on the - strand, i.e. a CG motif not a CHG motif.

```sh
+   CAG		CCG		CTG
-   GCT		GGC		GAC
```

> a model in which the two DNMT3A active sites are separated by approximately one helical DNA turn, which suggests that each tetrameric complex could simultaneously methylate two cytosine residues at a defined spacing of 8–10 bp

> In vivo, the spacing of CG dinucleotides at many DMRs is also consistent with an ~8–10 bp periodicity51, as is the finding that CG dinucleotides at an 8 bp spacing are overrepresented across the human genome68,69 and, to a lesser extent, across the mouse genome68. As DNMT3A seems to be a non-processive DNA methyltransferase70, the formation of an oligomer could help to explain the observed periodic pattern of DNA methylation.

## \cite{Jurkowska:2011dt}

> Active genes usually show hypomethylation around the transcriptional start site (TSS) and high levels of methylation in the gene body, which increase with gene expression.[30, 31, 46] The methylation of the gene body is believed to block aberrant transcription initiation inside the gene and therefore to help in avoiding the production of truncated mRNAs and proteins.

> Dnmt3L hence has a dual, gender-specific role: it is essentialfor the establishment of genomic imprints in oocytes and forthe silencing of the dispersed repeated sequences in malegerm cells.

> In the crystal structure, the active sites of the two central Dnmt3a subunits are separated by one DNA helical turn, which corresponds to around 10 bps, suggesting that two CG sites on __opposite strands__ separated by 10 bps could be methylated by Dnmt3a in one binding event. Additionally, Dnmt3a polymerises on DNA in vitro, forming Dnmt3a-DNA filaments, which positions the adjacent methyltransferase molecules that interact with the same DNA strand at a distance of 8–10 bps,[169] suggesting that these sites could become comethylated by Dnmt3a (Figure 3 B). Indeed, in vitro methylation experiments demonstrated that there is a correlation of methylation between two sites localised ~ 10 bps apart both in the same strands and in opposite strands.

> Additionally, the presence of 5–10 % of hemimethylated sites in the genome[210, 211] cannot be explained by the simple copy model alone, in which every hemimethylated CG site is faithfully and quickly remethylated by Dnmt1 after each round of replication. These considerations led to the development of stochastic models for DNA methylation,[212,213] based on differential equations taking into account the probabilities for methylation of a hemimethylated site (fhm) and methylation of an unmethylated site (fum).


> In summary, it appears that the “maintenance DNA methylation” refers to the preservation of average levels of DNA methylation at certain regions, but not to an accurate copying of site-specific DNA methylation patterns.

## \cite{Lacey:2009by}

> Excessively long context ranges might render the model ineffective, while choosing multiple breakpoints along a sequence would assume a piecewise model that is not likely to be biologically realistic.

I should address the "piecewise model ... unrealistic" issue in my discussion of `methsim`.

# Review of simulation methods

### \cite{Feng:2014iq}

1. Simulate group-specific $B_{i, j}$ from a Beta distribution with hyperparameters given by:
  - Group-specific means sampled from real data.
  - Group-specific dispersion parameters sampled from real data.
2. Sample sample-specific $u_{i, j} + m_{i, j}$ from same RRBS experiment.
3. Simulate sample-specfic $\beta_{i, j}$ from Beta-Binomial distribution that incorporates (1) and (2).

### \cite{Lacey:2013iy}

Simulate RRBS and not WGBS. All parameters estimated from a single sample.

0. Simulate CpGs from a point-process (HMM).
1. Simulate group-specific $B_{i, j}$.
  - (a) Simulate independent Beta-distributed RVs
  - (b) Induce correlations via variogram and re-estimate Beta-parameters.
  - Iterate (a) and (b) until stable.
2. Simulate sequencing-coverage as Gamma-mixture across bins. Includes across sample correlations of sequencing depth.
3. Simulate sample-specfic $\beta_{i, j}$ from Beta-Binomial distribution that incorporates (1) and (2).

### \cite{Sofer:2013bk}

Simulate microarray data by sampling from 539 TCGA 450k arrays.

A. Sample blocks of probes to preserve spatial correlation.
B. Induce differential methylation by weighted sampling of blocks containing "targets" and hopefully across other CpGs in the block since they are generally correlated.

### \cite{Xu:2013eg}

1. Simulate group-specific $B_{i, j_{k}}$
  - Beta
  - Truncated Normal
  - Truncated Normal-mixture
2. Simulate sample-specific $u_{i, j} + m_{i, j}$ from a round Normal distribution.
3. Simulate sample-specific $\beta_{i, j}$ from Binomial distribution that incorporates (1) and (2).

### \cite{Chen:2014jb}

No simulation of true values or sequencing coverage, but directly simulate $\beta$-values.

1. Simulate group-specific $B_{i, j_{k}}$
  - Uniform
  - Truncated Normal
  - Beta

Differential methylation is from a trend-with-category model.

### \cite{Chen:2014jb}

Simulate ERRBS and not WGBS.

1. Simulate group-specific $B_{i, j_{k}}$ from Beta distribution
2. Sample sample-specific $u_{i, j} + m_{i, j}$ from ERRBS experiment.
3. Simulate sample-specific $\beta_{i, j}$ from Binomial distribution that incorporates (1) and (2).

### \cite{Dolzhenko:2014bo}

1. Simulate group-specific $B_{i, j_{k}}$ from Beta distribution
2. Sample sample-specific $u_{i, j} + m_{i, j}$ from real data.
3. Simulate sample-specific $\beta_{i, j}$ from Binomial distribution that incorporates (1) and (2).

### \cite{Chen:2013eh}

No simulation of true values or sequencing coverage, but directly simulate $\beta$-values.

1. Simulate group-specific $B_{i, j_{k}}$
  - Beta
  - Normal distribution
  - Normal-mixture
