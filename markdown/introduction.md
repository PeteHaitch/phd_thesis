# Introduction 
## Chapter Overview
In this chapter I summarise the basic biology necessary for understanding my thesis, such as DNA, genomics and epigenetics. I introduce DNA methylation, describe some of its biological roles and discuss the types of questions asked in experiments studying DNA methylation. I describe assays for studying DNA methylation, with a particular focus on the whole-genome bisulfite-sequencing assay known as _methylC-seq_. I also outline the typical bioinformatics strategies for analysing bisulfite-sequencing data. Finally, I give an overview of the layout of my thesis.

## Basic biology
I describe the basic biology of DNA and the role of DNA in genetics and epigenetics.
 
### DNA
_DNA_, a much-needed abbreviation of deoxyribonucleic acid, is the molecule that encodes the genetic instructions for all known living organisms. __INSERT DNA DEFINITION FROM TEXTBOOK__ 

#### The DNA double helix {-}

The double helix is the most common, and most famous, structure of DNA (__REFERENCE__). In the double helix, the two strands of DNA run in opposite directions to each other and are therefore _anti-parallel_. One strand is called the 5' strand (pronounced "5 prime strand") - also known as the sense strand, Crick strand or top strand - and the other strand called the 3' strand (pronounced "3 prime strand") - also known as the antisense strand, Watson strand or bottom strand. Along each strand of the double helix are the four DNA nucleobases (_bases_): adenine (_A_), cytosine (_C_), guanine (_G_) and thymine (_T_). These bases form _complementary_ base pairings, _A_ with _T_ and _C_ with _G_, along the DNA double helix [__FIGURE OF DNA BASE PAIRING__](http://en.wikipedia.org/wiki/File:DNA_Structure%2BKey%2BLabelled.pn_NoBB.png).

#### DNA replication {-}

DNA _self-replicates_ to produce copies of itself. DNA replication occurs during cell division so that each daughter cell contains the same DNA as its parent. During DNA replication, the two strands are separated and each strand's _complementary_ DNA sequence is copied by an enzyme called _DNA polymerase_. It is because the two strands of DNA are complementary that ensures the daughter cell contains the same DNA sequence as the parent cell[^dna_replication].

[^dna_replication]: This of course ignores errors in the replication process. Such errors are very rare events but because DNA replication happens __HOW FREQUENTLY__ these events do occur. There are error-correcting processes that reduce the chance that such an error is retained in the daughter sequence, however, these are not perfect. Hence errors in DNA replication are one source of what are known as _mutations_ in the DNA.

DNA would not be quite so interesting if all it did was self-replicate. But it has an extra trick up its sleeve; it also acts as the template for the creation of new molecules that ultimately lead life as we know it. 

#### Genes and DNA transcription {-}

A gene is a sequence of DNA that is _transcribed_ to produce a functional product called ribonucleic acid or _RNA_ (Source: [http://sandwalk.blogspot.com.au/2007/01/what-is-gene.html](http://sandwalk.blogspot.com.au/2007/01/what-is-gene.html)). 

It is important to note here than not all DNA is transcribed into RNA. That is not to say that untranscribed DNA is not important. For instance, there are untranscribed _regulatory sequences_ of DNA that determine whether a nearby gene is transcribed. There is also _junk DNA_, _garbage DNA_, __etc.__ (__SOURCE__).

Conversely, not all transcribed DNA is a gene. DNA transcription is "permissive" and there are many DNA sequences that are transcribed by accident or in error.

#### RNA and RNA translation {-}

RNA, like DNA, is a nucleic acid. However, unlike DNA, RNA is usually single-stranded and has uracil, rather than thymine, as one of its four nucleobases.

RNA is _translated_ into proteins, which are the the molecules that "act out" the function of a cell. In eukaryotes it is  _messenger RNA_ (mRNA) that is processed and translated into proteins. Just as not all DNA is transcribed, not all mRNA is translated into protein. Furthermore, there are three other classes of RNA inside cells in addition to mRNA: transfer RNA (tRNA), ribosomal RNA (rRNA) and broad category called "small RNAs" (Source [http://sandwalk.blogspot.com.au/2013/05/scientific-authority-and-role-of-small.html](http://sandwalk.blogspot.com.au/2013/05/scientific-authority-and-role-of-small.html)).

The _genetic code_ describes how triplets of mRNA, or DNA, bases encode for _amino acids_. Amino acids are in turn the building blocks of _proteins_. Each triplet of bases describes one of the $3^4 = 64$ possible  _codons_. There is redundancy in the code. For instance, `GCU`, `GCC`, `GCA` and `GCG` all code alanine. __FIGURE/TABLE__ shows the complete genetic code.

#### The Central Dogma of molecular biology {-}
__KILL?__
Francis Crick, who along with James Watson and Rosalind Franklin discovered the double helix structure of DNA, proposed the Central Dogma of molecular biology (Source: [http://sandwalk.blogspot.com.au/2007/01/central-dogma-of-molecular-biology.html](http://sandwalk.blogspot.com.au/2007/01/central-dogma-of-molecular-biology.html)):
> ... once (sequential) information has passed into protein it cannot get out again (F.H.C. Crick, 1958)

Crick re-stated this in 1970:
> The central dogma of molecular biology deals with the detailed residue-by-residue transfer of sequential information. It states that such information cannot be transferred from protein to either protein or nucleic acid. (F.H.C. Crick, 1970)

Crick further added (1970):
> It is not the same, as is commonly assumed, as the sequence hypothesis, which was clearly distinguished from it in the same article (Crick, 1958). In particular, the sequence hypothesis was a positive statement, saying that the (overall) transfer nucleic acid → protein did exist, whereas the central dogma was a negative statement saying that transfers from protein did not exist.

Lewin’s GENES VIII defines the Central Dogma of Molecular Biology as,
> The central dogma states that information in nucleic acid can be perpetuated or transferred but the transfer of information into protein is irreversible. (B. Lewin, 2004)

### Genetics
Within a cell, DNA is packaged into units called _chromosomes_. A chromosome is a structure made up of DNA, RNA(__?__) and proteins. The combination of DNA and proteins is called _chromatin_.  The fundamental repeating unit in eukaryotic chromatin is the _nucleosome_. A nucleosome is comprised of DNA wrapped around a core of _histones_, which are proteins that act as kind of "spool" that the DNA "thread" wrap around. The beautiful compaction of DNA into chromosomes is illustrated in __CLASSIC FIGURE OF DNA -> NUCLEOSOMES -> CHROMATIN -> CHROMOSOMES__.

Chromatin is classified into two classes: _heterochromatin_ or _closed chromatin_ is tightly packed whereas as _euchromatin_ or _open chromatin_ is less tightly packed. DNA polymerase has easier access to euchromatin and therefore expressed genes are usually found in euchromatic regions of the genome. 

Humans are _diploid_ organisms, meaning that we have two copies of each chromosome in a typical cell[^somatic_cell]. We inherit one chromosome of each pair from our mother and one from our father. A typical human cell has 23 pairs of chromosomes: 22 _autosomes_ and 1 pair of _sex chromosomes_. This 23 pairs of chromosomes comprise what is known as the "human genome"[^mtDNA]

[^somatic_cell]: A sperm or egg cell is haploid and has a single (recombined) copy of each chromosome. 

[^mtDNA]: A typical human cell also contains the the mitochondrial DNA (mtDNA), which is maternally inherited. The mitochondrial chromosome is very small but there are 100-10,000 copies of it in each cell (__SOURCE__). The mtDNA is also considered to be part of the human genome. 

The length of a genome is the number of DNA base pairs (bp) in a haploid copy of that genome, that is, its length considering only a single copy of each chromosome. The human genome is over three billion base pairs long (__SOURCE__). 

Every person, except genetically identical twins, has their own unique genome. However, the genomes of any two randomly selected people are identical at __WHAT PERCENTAGE__ of sites.  This enabled the definition and construction of a human _reference genome_, which was completed in __WHEN__ (__SOURCE__). The vast majority of the human genome, as well as the vast majority of the genomes of other organisms, is non-coding. Much of the non-coding DNA is junk DNA. Indeed, over __WHAT PERCENTAGE__ of the human genome is comprised of repetitive DNA. 

* Each cell has the same genome but different gene expression and different protein expression

Genetic variation can lead to _phenotypic_ variation. A phenotype is the composite of an organism's observable characteristics or traits, such as its morphology, development, biochemical or physiological properties, phenology, behavior, and products of behavior (such as a bird's nest) (Source: [http://en.wikipedia.org/wiki/Phenotype](http://en.wikipedia.org/wiki/Phenotype)). In medical genetics, the phenotype being studied is usually a patient's disease-status or some proxy such as blood pressure. 

Genetic variation is not the only source of phenotypic variation. Environmental factors such as diet and pollution will also affect an individuals's phenotype. More recently, the role of _epigenetic_ variation has been a hot, and contentious, topic.

### Epigenetics
"Epigenetics" is a real Humpty Dumpty phrase: each author seems to believe that, "When I use a word, ... it means just what I choose it to mean - neither more nor less" [@carroll1897through]. I am sure to continue in this grand tradition.

Conrad Waddington coined the phrase in 1942 as a portmanteau of the words "epigenesis" and "genetics" (Source: [http://en.wikipedia.org/wiki/Epigenetics](http://en.wikipedia.org/wiki/Epigenetics)). He meant epigenetics as the study of how _genotypes_ give rise to phenotypes during development [@waddington1957strategy].

A popular contemporary definition of epigenetics is attributed to the epigeneticist (__CHECK SOURCE AND QUOTE__) Art Riggs: "[epigenetics is] the study of mitotically and/or meiotically heritable changes in gene function that cannot be explained by changes in DNA sequence[^mitotic_meiotic]" (p1 @russo1996epigenetic). This is in line with the _epi_ prefix being derived from the Greek word for "upon, near to or in addition" (__INSERT OXFORD DICTIONARY REFERENCE__), hence the idea that epigenetics encodes information "on top of" the DNA sequence. However, it is quite different to Waddington's original definition.

[^mitotic_meiotic]: Mitotically heritable means heritable during cell division and meitotically heritable means heritable during sexual reproduction.

More recently, the definition of epigenetics has taken on a "more biochemical flavour" ([@Daxinger:2010bu]). It is used to desribe marks whose heritably is yet to be proven, such as histone modifications (__SOURCE__). Sir Adrian Bird, an esteemed British geneticist, attempts to unite these definitions in the following ([@Bird:2007em]):

> The following could be a unifying definition of epigenetic events: the structural adaptation of chromosomal regions so as to register, signal or perpetuate altered activity states. This definition is inclusive of chromosomal marks [e.g. histone marks], because transient modifications associated with both DNA repair or cell-cycle phases and stable changes maintained across multiple cell generations qualify.

Regardless of which definition you subscribe to, an _epigenetic mark_ is the modification that causes this "epigenetic change". Actually, "causes" may be too strong a phrase as much of current epigenetics research is in identifying associations rather than causations and the question of whether the epigenetic mark is a cause or consequence of the ascribed function is oft-debated.

The _epigenome_ of a cell is the set of epigenetics marks present on the cell's genome. In contrast to the genome, which is identical between cells within an individual, the epigenome is highly variable between cells within an individual. Indeed, we can identify variation in a single epigenetic mark within cells of the same cell type from the same individual (__EXAMPLE__).

Methylation, acetylation, ubiquination and phosphorylation of histones, which alter the chromatin structure, are examples of the "more biochemical flavour" of recent epigenetics. But the classic epigenetic mark is DNA methylation, which is described in the next section.

## DNA methylation
DNA methylation is a chemical modification of DNA that can impart information on top of the DNA sequence. It is heritable during mitotic cell division, which means that it is faithfully copied across to the daughter cell during cell division[^dnam_copying]. It therefore fits into Art Riggs' definition of an epigenetic modification [@russo1996epigenetic]. So pervasive is DNA methylation that [@Lister:2009boa] dubbed it the "fifth base" of the DNA code. 

[^dnam_copying]: In practice this copying is not as faithful as, say, the copying of DNA from the parent to the daughter strand. Furthermore, the faithfulness of this copying will be different in different conditions, such as in a healthy, well-differentiated liver cell compared to a cancerous liver cell. Nevertheless, the copying of DNA methylation is faithful enough for most biologists to consider it as a mitotically heritable mark, most of the time (__TODO__: Is this just true of 5mC in mammals or more generally?)

What is referred to as DNA methylation is typically methylation of the cytosine nucleotide. Cytosine methylation is by far the most common form of DNA methylation in the animal and plant kingdoms[^other_DNAm], however, I will continue to use the terms _DNA methylation_ and _cytosine methylation_ interchangeably as is standard in the literature. 

[^other_DNAm]: Two additional forms of DNA methylation, N6-Methyladenine (m6A) and N4-methylcytosine (m4C), are primarily found in prokaryotes (__SOURCE: http://en.wikipedia.org/wiki/DNA_methyltransferase#cite_note-21__) (or is it bacterial (__SOURCE: http://www.ncbi.nlm.nih.gov/pubmed/16479578__))

A German chemist, W.G. Ruppel, first identified a methylated nucleic acid in 1898. Ruppel was studying _tuberculinic acid_, the poison of _Mycobacterium tuberculosis_[^tuberculosis] , and discovered that it contained a methylated nucleotide (__SOURCE__). In 1925, Johnson and Coghill isolated this methylated nucleotide as a product of hydrolysis of tuberculinic acid [@Johnson:1925js]. However, Johnson and Coghill's results were disputed for over twenty years by other researchers who were unable to replicate the original findings [@VISCHER:1949ty]. 

In 1945, Hotchkiss ultimately proved Johnson and Coghill correct when he isolated 5-methylcytosine from nucleic acid prepared from cow thymus [@HOTCHKISS:1948va]. Using paper chromatography, Hodgkiss demonstrated that methylated cytosine existed and was distinct from conventional cytosine and uracil.

[^tuberculosis]: _Mycobacterium tuberculosis_ was then known as _Tubercle bacillus_.

In the next section I focus on DNA methylation in mammals, in particular humans, but include brief descriptions of DNA methylation in other organisms. 

### DNA methylation in mammals
DNA methylation is essential for normal development in mammals, such as human and mouse (__REFERENCE__). My thesis focuses on DNA methylation in humans (_Homo sapiens_) but much of the material in this section is true for other mammals.

The typical site of DNA methylation is at the C5 carbon position of a cytosine nucleotide. This modified cytosine is called _5-Methylcytosine_ or _5mC_. [__FIGURE__](http://en.wikipedia.org/wiki/File:5-Methylcytosine.svg) shows the structure of 5mC.

At the level of single-stranded DNA, methylation is a binary event: a cytosine is either methylated or unmethylated. For double-stranded DNA, methylation is usually symmetric across the two strands (_strand symmetric_) although _hemimethylation_ does occur. Within a diploid cell a particular cytosine may be unmethylated or methylated on both homologous chromosomes or methylated on one chromosome and unmethylated on the other. While the former that is more common, examples of the latter case, such as _allele-specific methylation_ and _genomic imprinting_, are of considerable biological interest (__SOURCE/EXAMPLES__).

In mammalian genomes, most cytosines are unmethylated except for those at the _CpG dinucleotide_, which is usually strand-symmetric. I explain the importance of CpG dinucleotides in the next section. 

### CpG dinucleotides

A _CpG dinucleotide_, or more simply a _CpG_, is a cytosine followed by a guanine in the linear DNA sequence (__FIGURE__). The "p" stands for the phosphate backbone of DNA and some authors omit it in favour of simply calling it CG methylation. Approximately 60-80% of CpGs are methylated in mammals, meaning that the cytosine in the dinculeotide is a 5-Methylcytosine (__REFERENCE__). I will use the abbreviation _mCpG_ to refer to a methylated CpG, _uCpG_ to refer to an unmethylated CpG and simply _CpG_ when the methylation state is irrelevant.

CpGs are underrepresented in the human genome. The GC-content of the human genome, which is defined as the percentage of nucleotides that are either guanines or cytosines, is approximately __WHAT %__. Under a random sampling model of nucleotides we would expect approximately __WHAT %__ of the human genome to be comprised of CpGs. In fact, only __WHAT %__ of the dinucleotides are CpGs. __FIGURE__ shows the frequency of each dinucleotide in the human reference genome (hg19)[^cpg_frequency].

[^cpg_frequency]: See __APPENDIX__ for the `R` code used to compute these values for the human reference genome (build hg19).

One reason for this relative scarcity of CpGs is that methylated cytosines can spontaneously deaminate to thymines ([@Scarano:1967ue]). Thus, over time, many mCpGs will become TpGs, leading to a genome-wide reduction in the proportion of CpGs and a genome-wide increase in the proportion of TpGs (see __FIGURE__).

There are many other evolutionary pressures on the distribution of nucleotides in a genome. One effect of this is that the distribution of CpGs is far from uniform. I discuss this distribution in the next section.

### CpG islands and other sandy metaphors

One way to explore the distribution of CpGs in the human genome is to look at the distribution of distances from one CpG to the next. __FIGURE__ is a kernel density plot of the distances between neighbouring CpGs  and __FIGURE__ is the empirical cumulative distribution function of the distances between neighbouring CpGs for the human reference genome (hg19)[^cpg_distance].

[^cpg_distance]: See __APPENDIX__ for the `R` code used to generate these figures.

As can be seen from __FIGURE__, there are a subset of CpGs that are within __HOW MANY__ bases of the next CpG.

The 20-40% of CpGs that are methylated in mammalian genomes are mostly found in clusters called _CpG islands_ (_CGIs_). Again, there are evolutionary arguments for why these CGIs form (__REFERENCE__).

CGIs are important regulatory elements in the genome and are where most differences in DNA methylation between different cell types are found (__REFERENCE__). [@GardinerGarden:1987we] gave the classical definition of a CGI:

> For the purpose of this survey, regions of DNA with a moving average of %G+C over 50 and Obs/Exp CpG over 0.6 have been classed as CpG-rich regions. CpG-rich regions over 200 bp in length are unlikely to have occurred by chance alone, so, as a working definition, have been labelled as CpG islands.

Other authors have since provided alternative definitions (__REFERENCES__ Adrian Bird, Rob Klose, Hao Wu, etc.)

[@Wu:2010do] developed a Hidden Markov Model to define CGIs. In my thesis I have used the predicted CGIs from this model. These predicted CGIs are available from [http://rafalab.jhsph.edu/CGI/](http://rafalab.jhsph.edu/CGI/).

The "sandy/beachy" metaphors have been continued (i.e. stretched to breaking point) with various authors defining CGI shores, CGI shelves, CpG canyons, CpG deserts and CpGs in the "open sea". CGI shores [@Irizarry:2008hg] are defined as regions within 2kb of CGIs. They regions have been demonstrated to have an increased variability of CpG methylation [@Irizarry:2008hg]. CGI shelves (__WHO__) are defined as __WHAT__ and __WHY ARE THEY USEFUL__. CpG deserts are __WHAT__, CpG canyons are __WHAT__ and the "open sea" is __WHAT AND HOW DOES THIS DIFFER FROM CpG DESERTS__. 

### Non-CpG methylation

In humans, cytosine methylation in most cell types is found almost exclusively at CpGs. There are, however, certain cell types with widespread non-CpG methylation. Non-CpG methylation is often classified as CHH methylation or CHG methylation, where H stands for any nucleotide except G (__SOURCE IUPAC__). The rule-of-thumb is that non-CpG methylation is rare in somatic cells but common in pluripotent cells. Of course, there are exceptions to every rule, especially in biology.

To give a few examples, [@Lister:2009hy] found that in a somatic cell line (_fibroblasts_ or skin cells) $99.98\%$ of methylcytosines occured at CpGs whereas as in an embryonic stem cell line $24.5\%$ of methylcytosines occured in a non-CpG context. A subsequent paper from the same group ([@Lister:2011kg]) extended this result when they reported that, more generally, non-CpG methylation accounts for $20-30\%$ of methylcytosines in _pluripotent_ cell lines (which includes embryonic stem cell lines along with induced pluripotent stem cell lines, see __SECTION___). An exception to the rule is provided by [@Lister:2013et], where they found that neurons, a somatic and not a pluripotent type of cell, also have non-CpG methylation, albeit at a lower level ($1.3-1.5\%$ of all non-CpG cytosines were methylated).

Non-CpG methylation is less well studied than CpG methylation. This is partly due to the design of popular assays for studying cytosine methylation; for example the popular Illumina 27k and 450k beadchips (see __SECTION__) only measure cytosine methylation at CpGs __CHECK__. However, recent technological advances (see __SECTION__) mean that cytosine methylation can be routinely assayed regardless of the sequence context. The biological role of non-CpG methylation is less well understood than that of CpG methylation (__SOURCE__).

### Modifications of a modification

Methylation is not the only chemical modification of cytosine nucleobase, although it is by far the most common. Further modifications, listed approximately from most frequent to least frequent, are 5-hydroxymethylcytosine (5hmC), 5-formylcytosine (5fC) and 5-carboxylcytosine (5caC). 

The biological significance of these marks is still being determined, in part because the assays for studying these are still in development (__SOURCE__) and because their relatively scarcity means that experiments to detect these modifications are more difficult and expensive. 

One genome-wide study of 5hmC found that less than $ 1\%$ of all assayed cytosines in mouse fetal cortex and adult cortex cells were hydroxymethylated [@Lister:2013et]. Most of the 5hmC was detected in the CpG context and, although the genome-wide level of 5hmC was low, the authors reported significant levels of 5hmC at particular cytosines in the genome.

[@Kriaucionis:2009bm] and [@Tahiliani:2009kl] discovered that the TET enzymes can convert 5mC to 5hmC, 5hmC to 5fC and 5fC to hcaC. This suggested a role for 5hmC, 5fC and 5caC in the process of removing 5mC marks, as is discussed in the next section.

### Writers, readers and erasers

A common analogy used in describing DNA methylation is that of "writers", "readers" and "erasers" [@Moore:2013in]. Writers catalyse the methyl group onto the DNA, readers recognise methylated DNA and erasers remove the methyl group from the DNA.

In mammalian cells, the writers are the DNA _methyltransferase_ (DNMT) enzymes. These DNMTs are commonly split into two groups, namely the maintainence methyltransferase (thought to be DNMT1) and the _de novo_ methyltransferases (thought to be DNMT3a and DNMT3b)[^DNMT]. DNA methylation is not preserved by the DNA replication machinery and so it is the role of DNMT1 to restore the methylation pattern on the daughter strand of DNA following DNA replication. In contrast, DNMT3a and DNMT3b lay down new methylation marks and are particularly active during development when there are widespread changes in DNA methylation (__SOURCE__). DNMT1 and DNMT3b appear to be essential for mammalian development as mouse knockouts[^KO] for either gene is embryonically lethal (__SOURCE: http://labs.genetics.ucla.edu/fan/papers/npp2012112a.pdf__) whereas mouse knockouts for DNMT3a are runted but survive for ~4 weeks after birth (__SOURCE: http://labs.genetics.ucla.edu/fan/papers/npp2012112a.pdf__).

[^DNMT]: DNMT2, now known as TRDMT1, was once thought to be a DNA methyltransferase, however, __http://www.sciencemag.org/content/311/5759/395__ showed that it in fact methylates a small RNA and not DNA. Another protein, DNMT3L, is homologous to DNMT3a and DNMT3b but does contain catalytic domain that is necessary for methyltransferase activity. Instead, DNMT3L is thought to stimulate the activity of DNMT3a and DNMT3b (__SOURCE__).

[^KO]: A knockout mouse for gene $X$ is a mouse that has been genetically engineered to remove or otherwise inactivate gene $X$. Mouse knockouts can be either heterozygous knockouts (one copy still of the gene is still present/active) or homozygous knockouts (both copies of gene absent/inactive). 

The readers of DNA methylation recognise methylated DNA. These readers can recruit additional proteins to the site of the methylated cytosine to perform a variety of functions related to gene expression. For example, the _MBD_ (methyl-CpG-binding domain) group of proteins bind to DNA containing a methylated CpG which then suppress gene expression by preventing transcription factor binding at that site (__SOURCE__). Another group, the _UHRF_ (ubiquitin-like, containing PHD and RING finger domain) proteins, help DNMT1 methylate _hemimethylated_ DNA, which is DNA where one strand is methylated and the other is not (e.g. the parent strand compared to the daughter strand following DNA replication).

The removal or erasure of DNA methylation, _demethylation_, is characterised as _passive_ loss or _active_ removal. Passive loss occurs when the maintainence methyltransferases do not efficiently perform their role of restoring DNA methylation following cell division. This leads to a gradual, stochastic and genome-wide loss of DNA methylation after multiple cell divisions. This form of passive demethylation, sometimes called replication-dependent demethylation, cannot explain observations of precise targeted demethylation (__EXAMPLE__) nor the two stages of rapid global demethylation that occur during development [@Wu:2014gw].

Active demethylation is currently a hot topic in epigenetics research. Multiple mechanisms have been proposed, and it is indeed likely that there are multiple ways to achieve active demethylation. These mechanisms were recently reviewed by [@Wu:2014gw], which I now briefly summarise: 

1. The direct removal of the methyl group from 5mC is considered unlikely due to the strong carbon-carbon bond between the methyl group and the cytosine. 
2. There is evidence that the DNA repair machinery can be co-opted to remove a methylated base or the surrounding region. The excised base or region is then repaired with unmethylated cytosines replacing 5mCs.
3. 5mC oxidation-dependent active DNA demethylation. This follows from the observation that the TET enzymes can iteratively oxidate a 5mC -> 5hmC -> 5fC -> 5caC reaction. The removal of 5hmC, 5fC or 5caC is biochemically "easier" than the removal of 5mC and could occur via a more efficient form of replication-dependent demethylation, direct removal of the oxidized methyl group or through the DNA repair machinery (__FIGURE OF TET-CYCLE__).

One question raised by point (3) is whether 5hmC, 5fC and 5caC are simply intermediate products in an active demethylation cycle or if they themselves are bona fide epigenetic marks.

### Function of DNA methylation: "normal", cancer and development


### DNA methylation in other organisms and kingdoms

#### Insects {-}

* Honeybees
* Fruit fly? DNA methylation's existance in other organisms, such as the common fruit fly (_Drosophila Melanogaster_), is an area of contention although recent studies suggest it is present albeit at very small quantities and only in very limited contexts [@Takayama:2014jp].
* General

#### Plants {-}

* Arabadopsis
* Non-CG methylation

#### Fungi {-}

#### Bacteria {-}

#### Archaea {-}


## Assays for studying DNA methylation
A challenge to measuring DNA methylation is that it is erased by standard molecular biology techniques, such as _polymerase chain reaction_ (PCR) and bacterial cloning, and it is not revealed by DNA hybridization assays [@Laird:2010iv]. Therefore, almost all assays of DNA methylation require one of the following _pre-treatments_ of the DNA:

1. _Enzyme digestion_
2. _Affinity enrichment_
3. _Sodium bisulfite conversion_

Each of these pre-treatments can then be combined with a second step of amplification or hybridisation. For example, 

1. Gel-based analysis
2. Sanger sequencing
3. Microarray hybridisation
4. Massively parallel sequencing

An exception to this classification scheme are a new class of assay that seek to directly "read" whether a position is methylated or unmethylated without requiring a pre-treatment of the DNA. For example, both [@Laszlo:2013kf, Schreiber:2013in] measure the change in current as a DNA molecule passes through a nanopore to infer whether a cytosine was methylated while __PACBIO AND ONT DO WHAT?__.

Almost all assays of DNA methylation measure a "population average" from a pool of hundreds or thousands of cells. For a diploid organism, this is an average over several different levels: the two DNA strands, the two homologous chromosomes within a diploid cell and the hundreds or thousand of cells used in the assay (see __FIGURE__). Hundreds or thousands of cells are required in order to have sufficient material as input for the assay. Assays that require only a single cell as input do exist but are still in development and not in widespread use (__SOURCE/EXAMPLES__).

The _resolution_ and _throughput_ are two key variables when choosing which assay to use for an experiment. The resolution of an assay is the scale on which DNA methylation can be measured[^resolution]. For example, a high resolution assay allows a researcher to quantify the level of DNA methylation at a single base whereas a low resolution assay might only allow for qualitative assessment (i.e. presence or absence) of DNA methylation at larger regions, such as CpG islands.

[^resolution]: Depending on the experiment and its aims, the resolution of an assay might instead be defined as the scale on which DNA methylation can be quantified or that which allows inference to address a specific hypothesis.

The throughput of an assay can be quantified in two ways. The first is per-sample throughput, which is how many measurements of DNA methylation are made per-sample[^throughput1]. This is typically what people mean when they describe an assay as being "high-throughput" or "low-throughput" and is the definition I use in the title of my thesis. Depending on your definition of "high", a high-througput assay will produce on the order of tens of thousands to billions of measurements per sample. The second definition of throughput is per-dollar, that is, "how many samples can I process within my budget?". To avoid confusion, __I should refer to the latter definition as something else__.

[^throughput1]: This might reasonably be argued as being a definition of resolution.

The choice of which assay to use for an experiment is a trade-off between resolution, per-sample throughput and per-dollar throughput. Experiments that use an assay with high-resolution and high per-sample throughput generally have fewer samples than experiments using a lower resolution assay or an assay with lower per-sample throughput. 

In this section I will describe each of the pre-treatments but will focus on the bisulfite-conversion assays. In particular, I will describe in detail the "gold standard" assay of DNA methylation, _whole-genome bisulfite-sequencing_ (WGBS), that combines the sodium bisulfite conversion pre-treatment with massively parallel sequencing to produce whole-genome maps of DNA methylation at base pair resolution.


### Enzyme digestion assays
_Restriction endonucleases_ are an important technique in molecular biology. These enzymes can preferrentially "cut/cleave/digest" DNA at particular sequence motifs. The motif at which a restriction enzyme cleaves DNA is called the _recognition motif_ or _restriction sequence_. The methylation of a position in the recognication motif can inhibit a restriction enzyme from cleaving the DNA. This can be used to design an assay to infer the methylation state of a DNA fragment.

For example, the recognition site of the restriction enzyme _HpaII_ is `CCGC`. However, _HpaII_ will only digest DNA when the second cytosine in the motif is unmethylated. The _HELP_ (_HpaII_ tiny fragment enrichment by ligation-mediated PCR) assay compares DNA digested by _HpaII_ to one digested with another restriction enzyme that has the same recognition motif but is methylation-insensitive (_MspI_) to identify _hypomethylated_ regions of a genome (__SOURCE__). 

Another popular restriction enzyme for studying DNA methylation is _McrBC_, "an enzyme with the unusual and desirable property of cutting methylated DNA promiscuously (recognition sequence $R^{m} C(N)_{55– 103} R^{m} C$)" [@Irizarry:2008hg]. 

Assays based on restriction enzyme were some of the first developed for studying DNA methylation. They were initially developed for studying particular loci although they have been extended to genome-scale analysis approaches [@Laird:2010iv]. Restriction enzyme assays have a relatively low resolution (__REASONS/SOURCE__) but may also be used as a first-pass to enrich a sample for methylated or unmethylated regions which are then assayed using a higher resolution technique.

### Affinity enrichment assays
Affinity enrichment assays compare measurements between an "enriched" version and an "input" (control) version of the same sample to infer the presence or absence of DNA methylation. This is very similar to the idea behind chromatin immunoprecipitation followed by microarray hybridisation (ChIP-chip) or sequencing (ChIP-seq). Enrichment may be done using a methylation-specific antibody, such as __EXAMPLE__, or a __unmethylation-specific(?)__ antibody, such as __EXAMPLE__.

Some examples of affinity enrichment assays for DNA methylation are the microarray-based MeDIP, mDIP and mCIP (all based on methylated DNA immunoprecipitation) and their sequencing-based relatives, MeDIP-seq and mDIP-seq. These are all low resolution assays since they are based on the enrichment of regional differences between the enriched and input samples. Furthermore, the bioinformatic analysis of data from these assays is complicated by the varying CpG density along the genome, which leads to different enrichment affinities for different regions of the genome. However, these assays are can provide a relatively cheap and efficient genome-wide assessment of DNA methylation [@Laird:2010iv}].

### Sodium bisulfite conversion assays
In the 1980s, __Hayatsu OR Wang et al.__ discovered that when _denatured_ (__DEFINE?__) DNA is treated with sodium bisulfite (__CHEMICAL SYMBOL__), unmethylated cytosines deaminate to uracils much faster than do methylated cytosines (__SOURCE__). The methylated cytosines are said to be "protected" from conversion to uracils. When the bisulfite treated is amplified by PCR the uracils are subsequently converted to thymines.

This discovery led to the development of assays for studying studying cytosine methylation based on the pre-treatment of DNA with sodium bisulfite, which are sometimes referred to as bisulfite-conversion assays (__FROMMER; CLARK__). These assays are all based on the idea of comparing the sequence of the untreated DNA to the sequence of the bisulfite treated DNA to infer the methylation state of all cytosines in the sequence by whether or not they were converted to uracil (thymine) following the bisulfite treatment (__see FIGURE__). 

#### Analysis platforms combined with bisulfite conversion {-}

Initial experiments based on the sodium bisulfite pre-treatment of DNA used Sanger sequencing of cloned PCR products, a very laborious task, which restricted experiments to studying a limited number of short segments of DNA (__HOW MANY & LONG?__). Although subsequent enhancements in the automation of Sanger sequencing improved the throughput of these assays, it was never going to be able to deliver a cost-effective, genome-scale assay of DNA methylation. The development of hybridisation microarrays provided an alternative way of analysing bisulfite-treated DNA that could provide cheap, genome-wide measurements of DNA methylation(__SOURCE__). 

Microarrays contain thousands, even millions, of short _oligonucleotide_ probes. Each probe is is designed to hybridise to a particular DNA sequence and emits a signal, such as flurousence, that can be measured to infer the strength of the hybridisation. Therefore, an (idealised) way to analyse DNA methylation with a microarray is to hybridise bisulfite-converted DNA to a microarray that contains probes for both the methylated and unmethylated versions of all sequences of interest. The relative methylation of each sequence can be inferred from the relative intensities of the "methylated probe" to the "unmethylated probe". Such an idealised experiment brushes over many complications including [@Laird:2010iv]: 

* The reduced complexity of bisulfite-converted DNA (from a 4-base alphabet to a largely 3-base alphabet) leads to decreased hybridisation specificity
* Sequences containing multiple cytosines require multiple probe versions in order to assay all possible methylation patterns
* This approach requires the design of organism-specific microarrays

The Illumina 450k (__INSERT FULL NAME__) provides a modern implementation of the idea of bisulfite pre-treatment of DNA followed by microarray hybridisation for studying human samples. This array assays __HOW MANY__ CpGs (__IN WHAT TYPE OF REGIONS__) and uses two different probe designs (__WHY?__). There are still substantial bioinformatics challenges in analysing data from these assays (__CITE: SNPs in probes, normalisation, etc.__).

As mentioned in __SECTION__, genomics research was revolutionised by the development of cheap high-throughput sequencing technology and the study of DNA methylation was no exception. In 2008, two papers were published describing methods for whole-genome shotgun sequencing of bisulfite-converted DNA using the nascent Illumina sequencing technology [@Cokus:2008fc, Lister:2008bh]. [@Cokus:2008fc] termed their approach _BS-seq_ while [@Lister:2008bh] called theirs _methylC-seq_. These techniques, along with methods that these with targeting of particular regions of the genome, are described further in __SECTION__.

One final approach to mention is the use of mass spectrometry to analyse bisulfite-treated DNA, as implemented in Sequenom's EpiTYPER (__SOURCE__). This platform can provide quantitative measurements of CpG methylation across hundreds of loci and multiple samples and is frequently used to validate findings discovered using other platforms (__SOURCE__).

### DNA kinetics assays
__TODO__

### Pros and cons of bisulfite-conversion assays
Bisulfite-conversion assays are considered the "gold standard" for studying DNA methylation. Currently, multiple cells are required as input although the minimum number required is being lowered and efforts are well produce single-cell, whole-genome DNA methylation maps using bisulfite-treated DNA (__SOURCE__).

A key advantage of bisulfite-conversion assays over enzyme digest and affinity enrichment assays is that it provide base pair resolution of DNA methylation. In fact, single-molecule, base pair resolution is even possible for short DNA sequences when bisulfite-treated DNA is analysed with sequencing[^microarrays]. 

[^microarrays]: Microarray hybridisation can provide base pair resolution but not single molecule resolution. The signal from a microarray-based experiment a sample-wide average since, for a given locus, the signal is relative to the proportion of DNA fragments in the sample methylated at that locus.

A recently discovered disadvantage is that bisulfite-conversion assays are unable to distinguish between 5mC and 5hmC because both modifications similarly protect cytosine from deamination to uracil (__SEE [@Yu:2012ej] for references__). In effect, the detection of 5mC, and all subsequent inference, is confounded with that of 5hmC. For most experiments this isn't much of a problem - most cells have very low levels of 5hmC and so there is little confounding - however, in certain experiments this needs a more careful approach. To address this issue, [@Booth:2012fh] developed _oxidative bisulfite-sequencing_ (oxBS-seq) and [@Yu:2012ej] developed _Tet-assisted bisulfite sequencing_ (TAB-seq) to separate 5mC and 5hmC detection. 

oxBS-seq specifically measures 5mC. The input DNA is oxidated by potassium perruthenate (KRuO$_{4}$), which converts 5hmC to 5fC, prior to bisulfite-treatment. Only 5mC is protected from conversion during the bisulfite-treatment, which effectively means that only 5mC remains to be detected at the sequencing stage. The level of 5hmC can be estimated by performing traditional bisulfite-sequencing and then "subtracting" the oxBS-seq signal (5mC) from the bisulfite-sequencing signal (5mC + 5hmC).

TAB-seq takes the opposite approach to oxBS-seq by specifically measuring 5hmC. The input DNA is treated with a $\beta$-glucosyltransferase, which converts 5hmC to $\beta$-glucosyl-5-hydroxymethylcytosine (5gmC), followed by TET oxidation. Only 5gmC is protected from TET oxidation, which effectively means that only 5hmC remains to be detected at the sequencing stage. The level of 5mC can be estimated by performing traditional bisulfite-sequencing and then "subtracting" the TAB-seq signal (5hmC) from the bisulfite-sequencing signal (5mC + 5hmC).

A disadvantage of bisulfite-conversion assays is that they require knowledge of the underyling DNA sequence in order to infer the methylation states of cytosines. This requires either a parallel experiment to sequence the target region(s) or reliance on a reference genome. When relying on a reference genome, the inference of the methylation state can be confounded by the DNA sequence of the sample (__SOURCE__). __FIGURE__ describes such an example.

The bisulfite-treatment of DNA can introduce biases and other problems. Three examples are _PCR-bias_, _incomplete bisulfite conversion_ and _DNA degradation_. _PCR-bias_ is the difference in amplification efficiency of methylated and unmethylated versions of the same DNA sequence [@Warnecke:1997eh]. _Incomplete bisulfite-conversion_ leads to cytosines being incorrectly inferred as 5mC (__SOURCE__). Finally, bisulfite-treatment can lead to _DNA degredation_, which, as the name suggets, leads to a loss in DNA quality or quantity and can make it particularly hard to sequence long DNA molecules (__SOURCE: Anecdotal__).

## Whole-genome bisulfite-sequencing

### Targeted BS-seq, e.g. RRBS, NimbleGen and Agilent's capture kits, MRE-seq

* "RRBS’s ability to interrogate a locus is dependent on its MspI-cut-site (CCGG) density and consequently measures 10%–15% of the CpGs in the human genome (Bock et al. 2010; Harris et al. 2010)" \citep{Stevens:2013hv}
* See \citep{Stevens:2013hv} for MRE-seq details and a good summary of DNA methylation assays



## Bioinformatics analysis of bisulfite-sequencing data
While there are many tools and methods for analysing bisulfite-sequencing data, there are four fundamental steps:

1. Data quality checking
2. Read mapping and post-processing of mapped reads
3. Methylation calling
4. Inference

These ideas will be familiar to anyone who analyses high-thoughput sequencing data, but each requires a "twist" to work with bisulfite-sequencing data. In this section I describe each step in broad terms. Most of my thesis is focused on step 4 of this process, with some work on step 3. A more detailed study, including comparisons of different software, is give in [@]Krueger:2012ks].

I only describe an analysis pipeline for data from the directional bisulfite-sequencing protocol, which is the standard and simpler protocol.

### Data quality checking
The first step in any analysis of high-throughput sequencing data is to perform a quality check of the data. Much of this is done visually by comparing summary graphs of the current sample(s) to previous "good" samples. As such, data quality checking relies on the judgement of, hopefully, experienced analyst.

The `FastQC` software (__CITE__) is a very useful tool for performing this first step. It produces summary graphs of many key measures such as base quality scores, read length distribution and sequence contamination. As `FastQC` is a general purpose tool, some of its output is not so useful, or at least must be cautiously interpreted, for bisulfite-sequencing data. For example, `FastQC` will report a __warning/error__ if the GC-percentage, the percentage of read bases that are guanine or cytosine, is __less than what value__. Due to the bisulfite-treatment, bisulfite-sequencing data will naturally have a very low percentage of cytosine bases (__around what%__) and therefore a low GC-percentage. So, a low GC-percentage is no cause for concern.

More important is the identification and removal of contaminating sequences. FastQC will screen a subset of the reads against a list of known, common contaminants. Examples of such contaminants include the sequencing adapters, PCR primers and __another example__.

Illumina adapter sequences are ligated to each DNA molecule in the library in order to perform Illumina sequencing. The sequencer can "read into" the adapter sequence, particularly when using paired-end sequencing of short DNA fragments such as those created in WGBS libraries. This means that some reads are a chimera containing sequence of interest (from the sample) and junk (from the adapters). This junk needs to be removed for two reasons:

1. Reads containing adapter contamination will generally not map to the reference genome, meaning these reads go to waste.
2. If they do map, then this will result incorrect inferences - the garbage in, garbage out maxim.

Using a tool such as `Trim Galore!` (__CITE__) or `trimmomatic` (__CITE__), the reads can be _trimmed_ to remove these contaminants. Reads might also be trimmed to remove low quality positions, which are common at the 3' end of reads, although this isn't as essential as trimming to remove contaminants.

### Read mapping and post-processing of mapped reads
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

### Methylation calling
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


## Outline of thesis

## General TODOs
* Nucleotide, residue, (nucleo)base; what's the difference and which term should I use?
	* A nucleotide is composed of a nucleobase (also termed a nitrogenous base), a five-carbon sugar (either ribose or 2-deoxyribose), and one or more phosphate groups \[Coghill, Anne M.; Garson, Lorrin R., ed. (2006). The ACS style guide: effective communication of scientific information (3rd ed.). Washington, D.C.: American Chemical Society. p. 244. ISBN 978-0-8412-3999-9.)\]
* Should I give more detail on Rafa's CGI HMM?
* Check how deep sub-sections should go and fix.
* Define eukaryotes, prokaryotes and viruses
* Hyphenation or not of "bisulfite-sequencing"?
* General term for methylC-seq, BS-seq, etc. "whole-genome bisulfite sequencing (WGBS)?"
* Is it most accurate to say that restriction enzymes "cleave", "cut" or "digest" DNA?
* Define/describe DNA sequencing and the genomics revolution in the section describing DNA (or reference genome section?)
	* Define "library"
* Fix up reference formating
* Note that when I refer to a read I mean either a single-end read or a paired-end read. I will use to read_1 and read_2 if I want to refer a particular read that makes up a paired-end read.
* \citet{Zhang:2013uu} credit Laurent et al. (2010) with BS-seq and Cokus with WGBS.
