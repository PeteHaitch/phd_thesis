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

__FIGURE__ provides a more detailed view of these four steps (__Figure should be a flow chart with the various sub-steps inluded, e.g. Inference = DMC, DMR, ASM, meth-SNP detection, etc.__)

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

Each of these aligners can report the output in the standard `SAM` format (__CITE__), although each mapper does so in a slightly different way, which makes difficult the development of downstream analysis tools.

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

* M-bias plots assume the first base of the mapped read is the first base of the sequenced read. __This isn't true if the read had 5' trimming__. This means, e.g., `--ignore_r1 5` might be ignoring cycles 1-5, 2-6, 3-7, etc., depending on whether the read had 5' trimming.

### Other biases

* Cell-composition 
* Age
* Global hypo- or hyper-methylation
* SNPs (Discussed in __Methylation calling__)
* Alignment artefacts

## Methylation calling
__TODO: Re-write.__

1. __Describe methylation calling__
2. __Compare reference-based methylation calling to Bis-SNP and sample-specific methylation calling.__ 
3. __Describe post-hoc filtering approaches.__


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
