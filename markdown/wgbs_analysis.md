# Analysis of bisulfite-sequencing data

## Notes from outline.md (plus other notes added along the way)

* Chapter Overview
* QC, mapping, post-processing
* Methylation calling
        * Reference-based calling
* Biases
        * M-bias
        * Cell-composition, age, global hypo- or hyper-methylation, SNPs, alignment artefacts
* Univariate tests of $\beta_i$

## Chapter overview

While there are several bisulfite-sequencing assays (see __CHAPTER__) and many different software tools, there are four fundamental steps when analysing bisulfite-sequencing data:

1. Data quality control checks
2. Read mapping and post-processing of mapped reads
3. Methylation calling
4. Inference

__FIGURE__ provides a more detailed view of these four steps (__Figure should be a flow chart with the various sub-steps included, e.g. Inference = DMC, DMR, ASM, meth-SNP detection, etc.__)

Steps 1, 2 and 4 will be familiar to anyone who analyses high-thoughput sequencing data, but each requires a "twist" to work with bisulfite-sequencing data. Step 3 is obviously unique to assays of DNA methylation but there are similarities to variant calling from DNA-seq. Throughout this chapter I concentrate on the methylC-seq protocol, which is a directional protocol and the standard whole-genome bisulfite-sequencing assay. All data used in my thesis were generated using this protocol (__CHECK__). Other assays, particularly targetted assays such as RRBS and Methyl-Seq, require some additional tweaks.

In this chapter I describe steps 1-2 in broad terms. An excellent review of these fundamental steps, including comparisons of different software, is given by \citet{Krueger:2012ks}. My thesis has focused on steps 3 and 4 of this process. In this chapter I review previous work on steps 3-4 and defer my contributions on these topics to subsequent chapters.

## Data quality control checks
The first step in any analysis of high-throughput sequencing data is to perform a quality control check of the data. Much of this is done visually by comparing summary graphs of the current sample(s) to previous 'good' samples. As such, much of data quality control checking relies on the judgement of the analyst.

The `FastQC` software (__CITE__) is a very useful tool for performing this first step. It produces summary graphs of many key measures such as base quality scores, read length distribution and sequence contamination. `FastQC` is a general purpose tool for performing quality control checks of high-throughput sequencing data. This means that some of its output is not so useful, or at least must be cautiously interpreted, for bisulfite-sequencing data. For example, `FastQC` will report a __warning/error__ if the GC-percentage, the percentage of read bases that are guanine or cytosine, is __less than what value__. Due to the bisulfite-treatment, bisulfite-sequencing data will naturally have a very low percentage of cytosine bases (__around what%__) and therefore a low GC-percentage. Therefore, a low GC-percentage is no cause for concern.

More important is the identification and removal of contaminating sequences. FastQC will screen a subset of the reads against a list of known, common contaminants. Examples of such contaminants include the sequencing adapters, PCR primers and __another example__.

Illumina adapter sequences are ligated to each DNA molecule in the library in order to perform Illumina sequencing. The sequencer can "read into" the adapter sequence, particularly when using paired-end sequencing of short DNA fragments such as those created in WGBS libraries. This means that some reads are a chimera containing sequence of interest (from the sample) and junk (from the adapters). This junk needs to be removed for two reasons:

1. Reads containing adapter contamination will generally not map to the reference genome, meaning these reads go to waste.
2. If they do map, then this will result incorrect inferences - the garbage in, garbage out maxim.

Using a tool such as `Trim Galore!` (__CITE__) or `trimmomatic` (__CITE__), the reads can be _trimmed_ to remove these contaminants. Reads might also be trimmed to remove low quality positions, which are common at the 3' end of reads, although this isn't as essential as trimming to remove contaminants.

## Read mapping and post-processing of mapped reads
Read mapping is complicated by the bisulfite-treatment of the DNA. Following bisulfite-treatment, the DNA fragments are now mostly composed of three bases rather than four, which means there are many more sequence mismatches between a read and its true mapping location. Simply using standard read mapping software and allowing for more mismatches would result in many reads mapping to multiple locations in the reference genome. Instead, a field of read mapping software dedicated to bisulfite-sequencing data has developed.

These bisulfite-sequencing read mappers take one of two approaches:

1. "Methylation-aware" mismatch penalties.
2. _In silico_ bisulfite-conversion of reads and reference genomes.

While "methylation-aware" mappers provide the highest efficiency, these suffer from a bias whereby methylated reads are preferentially mapped over unmethylated reads (__CITE FELIX__). This biases downstream inference and means that these mappers are generally less popular. I will instead focus on the _in silico_ bisulfite-conversion mappers, called as such because of a key step taken by these in the mapping process.

_In silico_ bisulfite-conversion mappers convert all cytosines to thymines (resp. guanines to adenines) of the forward (resp. reverse) strand from the reference genome. They then take each read and create two _in silico_ bisulfite-converted versions of it[^2_reads]: the CT-read replaces all residual thymines with cytosines and the GA-read replaces all residual guanines with adenines. The CT-read is mapped against the CT-genome and the GA-read is mapped against the GA-genome using a standard mapping tools such as Bowtie2 (__CITE__) or BWA (__CITE__). This process is illustrated in __FIGURE__.

[^2_reads]: Two versions are made because we don't know _a priori_ from which of the two strands the read originated.

Depending on the exact settings used, the mapper reports the "best" location of each read with respect to the two reference genomes. It reports the original sequence of the read in the output file so that the methylation status of each position can be inferred by comparing it to the corresponding reference sequence.

_In silico_ bisulfite-conversion mappers avoid the bias inherent in the "methylation-aware" mappers because all reads, regardless of methylation status, "look the same" to the mapper. However, they do suffer from a slight loss in mapping efficiency (__CITE FELIX__).

__TABLE__ lists some popular bisulfite-sequencing read mappers with the underlying mapping software listed in brackets:

* Bismark (Bowtie1 or Bowtie2, __CITE__)
* BWA-meth (BWA, __CITE__)
* BSMAP (SOAP, __CITE__)
* __OTHERS__

Each of these aligners can report the output in the standard `SAM` format (__CITE__), although each mapper does so in a slightly different way, which it difficult to develop downstream analysis tools.

### PCR duplicates

PCR amplification of the input DNA is a common step in creating a library for high-throughput sequencing. PCR amplification is often required to ensure that there is a sufficient amount of DNA for the sequencer to properly work. Unfortunately, it can introduce biases into the library that result in some molecules being over- or under-represented compared to their "true" frequency. This means that when we sequence the library that we might sequence multiple fragments that are all copies of the same original piece of DNA, which gives a biased sampling of our sample's genome. These multiply-sequenced fragments are called _PCR duplicates_.

In bisulfite-sequencing data, PCR duplicates containing methylation locus $i$ result in $\beta_{i}$ being a biased estimate of $B_{i}$. This is because the sequenced reads do not accurately represent the true methylation levels of the sample.

Generally, there is no way to tell based on sequencing data if a read is truly a PCR duplicate. However, it is relatively easy to identify __possible__ PCR duplicates (which are almost always inaccurately referred to as "PCR duplicates"[^pcr_dup]). Software to identify PCR duplicates include the `MarkDuplicates` function that is a part of the `Picard` library (__CITE__), the `rmdup` function that is a part of the `SAMtools` library (__CITE__) and `SAMBLASTER` \citep{Faust:2014hf}. These software use mapped reads from a `SAM/BAM` file as input. There are subtle differences in how each method calls PCR duplicates. As an example, this is how `Picard`'s `MarkDuplicates` function works:

[^pcr_dup]: The distinction between possible PCR duplicates and true PCR duplicates is rarely made, probably because the phrase is so clunky. Possible PCR duplicates are almost always referred to as PCR duplicates with the implicit assumption that the reader is aware that these very likely include falst positive calls. Consistent with the literature, I will use the term PCR duplicates when I refer to reads identified as PCR duplicates by some software. I will use _true_ PCR duplicates when I need to distinguish the two concepts.

> Essentially what it does (for pairs; single-end data is also handled) is to find the 5' coordinates and mapping orientations of each read pair. When doing this it takes into account all clipping that has taking place as well as any gaps or jumps in the alignment. You can thus think of it as determining "if all the bases from the read were aligned, where would the 5' most base have been aligned". It then matches all read pairs that have identical 5' coordinates and orientations and marks as duplicates all but the "best" pair. "Best" is defined as the read pair having the highest sum of base qualities as bases with Q >= 15 (Source: [http://sourceforge.net/apps/mediawiki/picard/index.php?title=Main_Page#Q:_How_does_MarkDuplicates_work.3F](http://sourceforge.net/apps/mediawiki/picard/index.php?title=Main_Page#Q:_How_does_MarkDuplicates_work.3F)).

Using these software will result in false positive calls because reads may satisfy this (or similar) criteria yet come from separate DNA fragments in the original sample. The false positive rate is a particular problem when a small "genome", which might be a targetted region of a genome such as exons or CGIs, is sequenced at high coverage. This can be thought of as an example of the pigeonhole princple, which states that if we have $m$ containers (positions where a read can start) and $n > m$ items (reads), then at least one container must contain more than one item (at least one position must have more than one read starting there).

We could make this mathematically more precise, but it doesn't give us a simple answer to the question, "should we remove possible PCR duplicates from bisulfite-sequencing data?". The unsatisfactory answer is, "it depends". A rough guide is that provided the average/median sequencing coverage of the "genome" is less than the fragment length[^pcr_fragment_length] then we expect few false positive calls.

[^pcr_fragment_length]: For single-end data, the "fragment length" in these calculations is the read length.

In practice, this means that for whole-genome sequencing data we can be fairly confident that possible PCR duplicates are _true_ PCR duplicates. However, for targetted sequencing, such as RRBS or amplicon sequencing, we are much less confident and may remove some of our signal if we remove possible PCR duplicates. Instead, for RRBS we might exclude regions with an "abnormally" high sequencing coverage \citep{Krueger:2012ks}. For amplicon sequencing we often can't afford to exclude possible PCR duplicates if, for example, the aim is to identify rare epialleles by very deep sequencing of a small region.

### M-bias

Ideally, the probability that a base is called methylated should be independent of the sequencing cycle. \citet{Hansen:2011gu} found that this is not the case and that in fact there is considerable bias towards the start (5') and end (3') of reads. They called this bias the _M-bias_.

M-bias can be identified by plotting the read-position methylation level ($rpml$), which is the proportion of reads that are methylated at each read-position, as a function of read-position (see __FIGURE__). These $rpml$ are computed separately for each methylation type and, for paired-end sequencing data, separately for `read_1` and `read_2`. If there is no M-bias then this plot should be a horizontal line. A "bend" or "spike" in this plot is evidence of M-bias. Furthermore, these lines, which indicate the average level of methylation in the sample for that methylation type, should be at the same height for `read_1` and `read_2`, although we see this is often not quite the case (see __FIGURE__).

If a sample was processed in multiple batches, then M-bias estimation (and methylation calling) should perhaps be performed separately for each batch and then combined, or in some manner that is "batch aware". For example, two libraries with DNA derived from the same cell line, but with separate library preparations and sequencing runs, will likely suffer from batch effects due to differences in the library preparations or differences with the sequencing runs. Unfortunately, the person analysing the data doesn't always know all the sample processing steps that may have introduced such batch effects and so these can be hard to deal with.

The strongest source of M-bias in Illumina WGBS data is at the 5' end of `read_2`, which sequences the 3' end of the DNA fragment. Because the DNA fragment is often shorter than the sum of the read lengths, the 3' end of the fragment often contains adapter sequence and other "junk" sequence. The adapter sequence may contain cytosine bases, which will be misinterpreted as evidence of methylation \citep{Krueger:2012ks}. Similarly, "fill-in cytosines" are used in the construction of RRBS libraries to repair the ends of DNA fragments after cleavage by MspI; these would also be misinterpreted as evidence of methylation \citep{Krueger:2012ks}. Another source of M-bias is incomplete or uneven bisulfite-conversion (__CITE__).

#### Estimating M-bias {-}

Estimating the M-bias and incorporating the results into the methylation calling can be done using two different strategies:

1. Compute the M-bias from the filtered data, then call methylation events incorporating both the filters and the M-bias. This strategy requires two passes over the `SAM/BAM` file - one to compute the M-bias and one to do the methylation calling. This is the approach used by `bismark_methylation_extractor`, `comethylation` and `Bis-SNP`.
2. Call methylation events from the filtered data but retain the read-position of each methylation event. Compute the M-bias from this first file and then filter out methylation events that suffer from M-bias. This strategy requires only a single pass over the `SAM/BAM` file but requires additional information to be stored alongside the methylation calls which is then followed by a pass over the first file containing the methylation calls. This is the approach taken by `BSmooth`.

I find the first strategy conceptually simpler, and easier to program, and so use it in my methylation calling software, `comethylation`. In theory, the first pass over the `SAM/BAM` could be "skipped" by computing the M-bias simultaneously with sorting the `SAM/BAM` and marking PCR duplicates. For example, using (idealised) Unix pipes:

```
# 'align' is read mapping software
# 'mark_duplicates' is PCR duplicate marking software. It creates 'out.bam'.
# 'mbias' computes mbias. It produces 'mbias.out' and (passes through) 'out.bam'.
align ref.idx r1.fq r2.fq | mark_duplicates | mbias --filters
# call_meth is methylation calling software. It produces 'call_meth.out' based on 'mbias.out' and 'out.bam'
call_meth --filters mbias.out out.bam
```

A somewhat subtle point is that M-bias should only be estimated from reads that are actually going to be used for methylation calling. Suppose M-bias is highly correlated with some other quality metric, such as base quality, so that positions with M-bias also have low base quality[^bq_mbias]. If you already intend to ignore all positions with a base quality less than 20 in your methylation calling then it makes sense to also ignore these positions when computing M-bias; otherwise you will overestimate the effect of M-bias and unnecessarily exclude read positions in your methylation calling. Unfortunately, the software I use for estimating M-bias, `bismark_methylation_extractor`, does not allow the user to exclude certain reads or read-positions from its calculations.  

[^bq_mbias]: This is very often the case since the 3' end of reads are typically of lower quality and also frequently suffer from M-bias.

```{r}
# TODO: Figure of M-bias computed from raw vs. filtered data
```

```{r}
Figure showing several M-bias curves with different patterns, e.g. 5' and 3' junk, within-read spikes/troughs, etc.
```

#### Pre-trimming reads destroys the one-to-one relationship between read-position and sequencing cycle {-}

Trimming reads prior to alignment (_pre-trimming_), such as using `Trim Galore!` to remove adapter sequence from reads, destroys the one-to-one relationship between the sequencing cycle and the read position. This causes a minor problem when computing M-bias because we no longer know whether the read position is identical to the sequencing cycle. Soft- or hard-clipping reads of their adapter sequence __during__ the alignment avoids this issue, because the clipping information (should be) preserved in the `CIGAR` string[^cigar].

[^cigar]: Many downstream tools, including the current version of `comethylation`, do not properly handle the information in the `CIGAR` string, particularly for soft-clipped reads. The current version of `comethylation` will skip a read (with a warning) that contains a `I` (insertion), `D` (deletion), `S` (soft-clip) or `H` (hard-clip). __FIX THIS__. This is a shortfall of the downstream tools and not of aligner-based clipping _per se_, but is nonetheless an issue in practice. __TODO: Update what changes are made to `comethylation`. The `bwa-meth` and `LAST` alignment software both perform well without pre-trimming of reads because they can soft-clip reads on the fly whereas other bisulfite-sequencing aligners, such as `Bismark`, cannot soft-clip reads and so require that reads are pre-trimmed (__CITE `bwa-meth` pre-print or paper__).__

The M-bias plot is based on the read-position from the aligned data and not the sequencing cycle (which isn't directly available in the `SAM/BAM` file). If the reads have been pre-trimmed, then each read-position in the M-bias plot will therefore contain data from multiple sequencing cycles, which can amplify or mask M-bias signals.

For example, suppose we performed $100$ bp single-end sequencing and pre-trimmed the first 20 bp of $90\%$ of the reads. Then, read-position $80$ will comprised of $10\%$ sequencing cycle $80$ and $90\%$ sequencing cycle $100$. It is very likely the sequencing cycle $100$ suffers more from M-bias than does cycle $80$, and so this will appear in the M-bias plot as M-bias at read position $80$. When using `bismark_methylation_extractor` this may necessitate that read-positions $80-100$ are ignored in downstream analyses, which would unnecessarily throw away $20\%$ of the data!

This problem cannot be avoided if reads are pre-trimmed because the trimming information is not preserved. \cite{Hansen:2012gr} suggest a separate M-bias plot for each read-length, which will help mitigate the effect of confounding between read-position and sequencing-cycle[^read_trimming]. However, if trimming is performed during the alignment then all the necessary information is retained and the x-axis of M-bias plots would be "sequencing cycle" rather than "read-position".

[^read_trimming]: Generally, this would also require that methylation calling is performed separately on each read length because most methylation callers are unable to deal with different M-bias profiles for different read lengths.

#### Identifying M-bias {-}

In practice, the M-bias curves for each sample are visually inspected to look for evidence of M-bias, i.e. read-positions that are "too far away" from the vast majority of read-positions. A minor problem with this approach is that the user may (unconsciously) use a different cutoff for each sample. For example, it is more difficult to determine an appropriate cutoff when there is a gradual slope in the M-bias plot than when there is a clear jump (__TODO: An example where there last few read-positions are more than `cutoff` from the median but due to a gradual slope vs. an example where a single position is more than `cutoff` from the median due to a spike__).

This motivated me to include a function, `trim.mbias` (__TODO: Check names used in final version. In particular, `trim.mbias` needs to be re-named to reflect that it identifies, rather than trims, M-bias.__), to identify M-bias in a more systematic way. `trim.mbias`, along with the related functions `read.mbias`[^read.mbias] and `plot.mbias`[^plot.mbias], are part of my `cometh` package.

[^read.mbias]: `read.mbias` reads in the `M-bias` files produced by `Bismark` and returns these as a `data.frame`.
[^plot.mbias]: `plot.mbias` creates a plot of M-bias stratified by the methylation type and, if the data are paired-end, whether the read is `read_1` or `read_2`.

`trim.mbias` identifies read-positions that display M-bias. This is done separately for each methylation type because the level of methylation varies widely between CG and non-CG methylation types. For paired-end sequencing experiments, it is also done separately for each of `read_1` and `read_2` because M-bias is very different for these two mates (see __FIGURE__) and also because `read_2` often has a slightly lower average level of methylation than does `read_1` (see __FIGURE and provide explanation__).

To do this, `trim.mbias` uses a simple normalisation of the read-position methylation levels ($rpml$). Specifically, the median level of methylation across all read-positions is subtracted from the read-position methylation level to create normalised read-position methylation levels (for example, $nrpml_{CpG}^{read_1} = rpml_{CpG}^{read_1} - median(rpml_{CpG}^{read_1})$, for CpG methylation in `read_1`).

`trim.mbias` identifies read-positions where the $rpml$ differs by more than a given value (`cutoff`) from $nrpml$. While it computes these statistics separately for each methylation type and read type, a common `cutoff` is used for all methylation types and read types. I tend to use a value of `cutoff = 0.02` (__TODO: confirm__), meaning that if the $median(rpml_{CpG}) = 0.75$, then any read-position with $nrpml_{CpG} < 0.73$ or $nrpml_{CpG} > 0.77$ is flagged as showing evidence of M-bias. It is currently left to the user to decide whether a read-position should be excluded if it displays evidence of M-bias in a single methylation type or only if it displays evidence across all methylation types. I tend to exclude read-positions that display evidence of M-bias in the methylation type I am working with, which is typically CpG methylation.

The results of `trim.bias` are returned as a `data.frame`:

```{r}
# TODO: Insert code snippet showing a call to `trim.mbias`
```

The results can also be visualised by specifying the `plot = raw` or `plot = normalised` arguments. `plot = raw` gives the same M-bias plot as `plot.mbias` except that outliers are now flagged. `plot = normalised` plots the $nrpml$ at each read position along with the a horizontal line indicating the chosen `cutoff`.

```{r}
# TODO: Insert code snippet showing plots returned by `trim.bias` with both `plot = raw` and `plot = normalised`
```

#### What to do about M-bias? {-}

Now that we've found read-positions with M-bias, what can we do about it? Typically, read-positions showing evidence of M-bias are excluded when calling methylation events. In fact, a slightly cruder procedure is often used whereby all read-positions preceding or following those with M-bias are excluded. Both `bismark_methylation_extractor` and `Bis-SNP`, two popular methylation callers, use this method. Generally, this is okay because M-bias tends to occur as runs of read-positions at the 5' and 3' ends of reads.

However, occasionally there are spikes in the M-bias plot, which indicate specific read-positions that we would like to exclude. __FIGURE__ shows an example of such a spike that occurred at __WHICH__ read-position, likely due to a boost in the laser for sequencing cycles $101-150$ (__CHECK: If can't find a source then rephrase__). Using `bismark_methylation_extractor` we would be forced to either retain this position or ignore read-positions $101-150$, effectively ignoring one third of the sequencing data, much of it unaffected by M-bias. My `comethylation` software is more flexible by allowing the user to exclude individual read-positions or ranges of read-positions via the `--ignore-read1-positions` and `--ignore-read2-positions` flags.

```{r, echo = FALSE}
# TODO: M-bias plot with an internal spike
```

Ideally, we would be able to remove the effects of M-bias by accounting for it in methylation calling rather than by simply excluding those read positions entirely. In several datasets I have analysed, more than __HOW MUCH__ of the entire sequencing data is lost to M-bias.

One idea is to weight methylation calls by the level of M-bias. A downside to this approach is that this would turn an otherwise binary methylation call into a continuous value between $0$ and $1$, with an attendant loss in interpretability. However, we could still compute $\beta$-values, and use these in downstream inferences, without much loss in interpretation. An additional downside to this approach is an increased computational complexity and cost, and this is why I have not further pursued this idea.

### Other biases

There are several other sources of bias when making methylation calls. Some are what I call technical biases - those driven by the sequencing protocol and processing of the data - while others I call biological biases - those driven by measuring something other than what you thought you were measuring.

Apart from M-bias, technical biases include sequencing and alignment errors, which can both give rise to spurious or otherwise incorrect methylation calls in obvious ways.

Biological biases include cellular heterogeneity \citep{Jaffe:2014gy} and sequence variation at or nearby methylation loci, such as CpGs that are also SNPs or with SNPs nearby. I discuss the latter in the section __METHYLATION CALLING__.

#### Cellular heterogeneity {-}

Strictly speaking, cellular heterogeneity exists once you study DNA from more than a single cell, as shown in a recent papers using single-cell bisulfite-sequencing \citep{Guo:2013iha,Smallwood:2014kn}, however in this context it means when you study DNA from a sample that is made up of different types of cells.

For example, many studies of DNA methylation use whole blood as the  sample tissue due to the ease with which it can be obtained. However, whole blood contains a mixture of cell types, each of which has a distinct methylation profile. \citet{Houseman:2012km} exploited this fact to estimate the cell composition of whole blood samples, in particular the differences in white blood cell (leukocyte) proportions between different sub-populations (e.g. cases and controls).

In contrast to \citet{Houseman:2012km} where this cellular heterogeneity was useful, in studies exploring the relationship between DNA methylation changes and phenotypes, this cellular heterogeneity can seriously bias the results and must be properly accounted for in any analysis \citep{Jaffe:2014gy,Houseman:2014ds}. \citet{Jaffe:2014gy} provide evidence that several reported relationships between age and DNA methylation are likely due to changes in the cell composition of whole blood with age and not due to DNA methylation changes _per se_.

Methods to estimate the cellular heterogeneity bias and adjust for it are available (e.g. \citet{Jaffe:2014gy,Houseman:2014ds,Zou:2014gc}), although they have mostly been applied to DNA methylation arrays and not bisulfite-sequencing data. This is not to say that these problems don't exist for sequencing data, merely that these have not been as well-explored.

## Methylation calling

Methylation calling is the process of calling each sequenced methylation locus as being either methylated or unmethylated[^no_call], as well as determining the _context_ or _type_ of each methylation event (i.e., CpG, CHG or CHH) based on the sequencing data and a reference DNA sequence. In principle, this is a simple process, however, this belies some complications, which I discuss in this section.

[^no_call]: A third possibility is making the call that the "methylation locus" is not in fact a methylation locus. For example, if the sequenced base at a methylation locus is an _A_ or a _G_ then this may be evidence that position is not in fact a methylation locus.

Most bisulfite-sequencing alignment software includes methylation calling software (e.g. Bismark includes `bismark_methylation_extractor`), which is run after the alignment and post-processing of the `SAM/BAM` file. An alternative is `Bis-SNP` \citep{Liu:2012ge}, which performs methylation calling, and also variant genotyping, from bisulfite-sequencing data aligned with the user's choice of alignment software. These aforementioned methods perform methylation calling of 1-tuples. In contrast, `comethylation` can perform methylation calling at m-tuples. `MethPat` (__CITE__) is also capable of performing methylation calling at m-tuples.

### Methylated or unmethylated?

At each cytosine in a reference sequence ($\mathcal{I}$), the base from each read aligned to that locus is compared - if the sequenced base at a methylation locus is a _C_ then call it as methylated, if the sequenced base at a methylation locus is a _T_ then call it as unmethylated, 


Determining whether a position is methylated or unmethylated requires that the sequenced base is compared to a reference sequence

All of these software perform reference-based methylation calling. That is, they require the specification of a DNA reference sequence that the aligned bisulfite-sequencing data are compared against to infer the methylation state of each sequenced locus. This is typically the reference genome used in the alignment step, which implicitly defines $\mathcal{I} :=\mathcal{I}_{ref}$, where $\mathcal{I}_{ref}$ is the set of methylation loci in the reference genome. Due to DNA variation between the sample and the reference genome, this assumption is obviously not true but it is often good enough when paired with a post-hoc filtering of the methylation calls (described below).


A more refined set of methylation calls can be made by incorporating sample-specific genetic variation into the methylation calling procedure. For example, `Bis-SNP` is calls genetic variation simulataneously with methylation calls at positions defined $\mathcal{I}_{Bis-SNP}$. Alternatively




As mentioned in __SECTION__, these results can be _post-hoc_ filtered to remove problematic sites to produce $\mathcal{I} \subset \mathcal{I}_{ref}$. There reference genome


This is typically the reference genome used in the alignment step, although this may be refined by including known sample-specific genetic variation


### Determining the context or methylation type


The reference sequence may also used to infer the _context_ or _type_ of each methylation event, that is, CpG, CHG or CHH. This is used by Bismark, `comethylation` and `MethPat`. Alternatively, the context could be inferred from the sequencing data rather than the reference genome by looking at the upstream bases


__UP TO HERE__

__TODO: Re-write.__

1. __Describe methylation calling__
2. __Compare reference-based methylation calling to Bis-SNP and sample-specific methylation calling.__
3. __Describe complications due to sequence variation at or nearby methylation loci.__
4. __Describe post-hoc filtering approaches.__

If using a reference-based methylation caller, such as `Bismark`, then we are implicitly defining $\mathcal{I} :=\mathcal{I}_{ref}$, where $\mathcal{I}_{ref}$ is the set of methylation loci in the reference genome. Due to DNA variation between the sample and the reference genome, this assumption is not true. As mentioned in __SECTION__, these results can be _post-hoc_ filtered to remove problematic sites to produce $\mathcal{I} \subset \mathcal{I}_{ref}$.

Alternatively, we might use a more sophisticated methylation caller, such as `Bis-SNP` to define $\mathcal{I}$.

However, it is often "good enough" for genome-level analyses, provided that the sample and the reference are not too far genetically diverged, but may lead to false conclusions when performing finer-scale analyses.


When relying on a reference genome, there are two different classes of problematic loci:

1. Reference-specific loci: Methylation loci that exist in the reference genome but not in the sample's genome
2. Sample-specific loci: Methylation loci that exist in the sample's genome but not in the reference genome

Reference-specific methylation loci can lead to false methylation calls whereas sample-specific methylation loci will be missed by a reference-based methylation caller or might be misintepreted as genetic variation rather than methylation.


if, for example, a CpG in the reference is a TpG in the sample.


Many methylation callers ignore these problematic loci and perform _reference-based_ methylation calling using $\mathcal{I}_{ref}$ (Bismark, __others?__). These calls may then be _post-hoc_ filtered to remove known problematic sites so that $\mathcal{I} \subset \mathcal{I}_{ref}$.

For example, it is commont to remove all loci that are also sites of known genetic variation, such as _single nucleotide polymorphisms_ (SNPs).


Overall, defining $\mathcal{I} :=\mathcal{I}_{ref}$ isn't a big problem and  However, for certain loci it is clearly an issue. [@Liu:2012ge] developed the `Bis-SNP` software, which performs more sophisticated methylation calling, that is, refining the definition of $\mathcal{I}$.

__TODO: Read Bis-SNP paper and summarise__

__NOT TRUE (SEE BISSNP)__: It is more-or-less impossible to identify sample-specific loci if DNA sequence data of the sample are not available. The reason is that it is generally not possible to distinguish sample-specific methylation from sample-specific genetic variation from bisulfite-sequencing data alone. __FIGURE__ shows why this is so difficult. These loci lead to false negative methylation calls, we miss these completely, and could lead to false positive genetic variation calls because we misinterpret methylation as genetic variation.

The effect of the reference-specific methylation loci can be more readily moderated. Reference-specific methylation loci lead to false methylation calls if, for example, a CpG in the reference is a TpG in the sample. Then all reads from the sample that map to this locus will have the T base, which can be misinterpreted as being an unmethylated cytosine. There is less of a problem if, for example, the sample has an ApG at this position; an A is not evidence for or against methylation at this loci and so a methylation caller should not be influenced by it.

A standard technique to remove most reference-specific loci is to ignore all positions in the reference genome that are known to be common sites of genetic variation, such as _single nucleotide polymorphisms_ (SNPs) (__CITE__). This is a conservative technique that will remove the vast majority of reference-specific loci from further analysis, regardless of whether the sample has a genetic variant at that position. This obviously requires a database of known variation for the organism being studied, which is the case for frequently studied organisms such as humans and mice.

A less conservative technique is to genotype the sample at these SNPs. __FIGURE__ shows how this can be done using only the bisulfite-sequencing data, provided the data are "directional". Briefly, suppose these is a CpG in the reference genome that is an ApG our sample.

To remove those reference-specific loci that are not found in databases we might identify loci in the sample that display a large number of non-C/T bases (resp. non-G/A bases) at a C (resp. G) on the forward (resp. reverse) reference strand.

* What if the downstream base is mutated
* Array-specific problems; see "On the analysis of the illumina 450k array data: probes ambiguously mapped to the human genome" and "Guthrie card methylomics identifies temporally stable epialleles that are present at birth in humans."

## Inference

### Identifying methylcytosines

* i.e. Lister's approach

### DMCs

### DMRs

### ASM

### Methylation entropy

### Epipolymorphism

### Epiallele detection

### Others

## TODOs
* See https://github.com/brentp/clustermodel/README.md for a nice, brief summary of the main DMR calling approaches.
