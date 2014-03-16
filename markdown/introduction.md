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

#### The Central Dogma of molecular biology {-)
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
"Epigenetics" is a real Humpty Dumpty phrase: each author seems to believe that, "When I use a word, ... it means just what I choose it to mean - neither more nor less" ([@Carroll:1897uu]). I am sure to continue in this grand tradition.

Conrad Waddington coined the phrase in 1942 as a portmanteau of the words "epigenesis" and "genetics" (Source: [http://en.wikipedia.org/wiki/Epigenetics](http://en.wikipedia.org/wiki/Epigenetics)). He meant epigenetics as the study of how _genotypes_ give rise to phenotypes during development ([@Waddington:1957ub]).

A popular contemporary definition of epigenetics is attributed to the epigeneticist (__CHECK SOURCE AND QUOTE__) Art Riggs: "[epigenetics is] the study of mitotically and/or meiotically heritable changes in gene function that cannot be explained by changes in DNA sequence[^mitotic_meiotic]" (p1 [@russo1996epigeneti]). This is in line with the _epi_ prefix being derived from the Greek word for "upon, near to or in addition" (__INSERT OXFORD DICTIONARY REFERENCE__), hence the idea that epigenetics encodes information "on top of" the DNA sequence. However, it is quite different to Waddington's original definition.

[^mitotic_meiotic]: Mitotically heritable means heritable during cell division and meitotically heritable means heritable during sexual reproduction.

More recently, the definition of epigenetics has taken on a "more biochemical flavour" ([@Daxinger:2010bu]). It is used to desribe marks whose heritably is yet to be proven, such as histone modifications (__SOURCE__). Sir Adrian Bird, an esteemed British geneticist, attempts to unite these definitions in the following ([@Bird:2007em]):

> The following could be a unifying definition of epigenetic events: the structural adaptation of chromosomal regions so as to register, signal or perpetuate altered activity states. This definition is inclusive of chromosomal marks [e.g. histone marks], because transient modifications associated with both DNA repair or cell-cycle phases and stable changes maintained across multiple cell generations qualify.

Regardless of which definition you subscribe to, an _epigenetic mark_ is the modification that causes this "epigenetic change". Actually, "causes" may be too strong a phrase as much of current epigenetics research is in identifying associations rather than causations and the question of whether the epigenetic mark is a cause or consequence of the ascribed function is oft-debated.

The _epigenome_ of a cell is the set of epigenetics marks present on the cell's genome. In contrast to the genome, which is identical between cells within an individual, the epigenome is highly variable between cells within an individual. Indeed, we can identify variation in a single epigenetic mark within cells of the same cell type from the same individual (__EXAMPLE__).

Methylation, acetylation, ubiquination and phosphorylation of histones, which alter the chromatin structure, are examples of the "more biochemical flavour" of recent epigenetics. But the classic epigenetic mark is DNA methylation, which is described in the next section.

## DNA methylation
DNA methylation is a chemical modification of DNA that can impart information on top of the DNA sequence. It is heritable during mitotic cell division, which means that it is faithfully copied across to the daughter cell during cell division[^dnam_copying]. It therefore fits into Art Riggs' definition of an epigenetic modification ([@russo1996epigeneti]). So pervasive is DNA methylation that [@Lister:2009bo] dubbed it the "fifth base" of the DNA code. 

[^dnam_copying]: In practice this copying is not as faithful as, say, the copying of DNA from the parent to the daughter strand. Furthermore, the faithfulness of this copying will be different in different conditions, such as in a healthy, well-differentiated liver cell compared to a cancerous liver cell. Nevertheless, the copying of DNA methylation is faithful enough for most biologists to consider it as a mitotically heritable mark, most of the time (__TODO__: Is this just true of 5mC in mammals or more generally?)

What is referred to as DNA methylation is typically methylation of the cytosine nucleotide. Cytosine methylation is by far the most common form of DNA methylation in the animal and plant kingdoms (__ARE THERE OTHERS?__), however, I will continue to use the terms _DNA methylation_ and _cytosine methylation_ interchangeably as is standard in the literature. 

A German chemist, W.G. Ruppel, first identified a methylated nucleic acid in 1898. Ruppel was studying _tuberculinic acid_, the poison of _Mycobacterium tuberculosis_[^tuberculosis] , and discovered that it contained a methylated nucleotide (__SOURCE__). In 1925, Johnson and Coghill isolated this methylated nucleotide as a product of hydrolysis of tuberculinic acid (@{Johnson:1925js}). However, Johnson and Coghill's results were disputed for over twenty years by other researchers who were unable to replicate the original findings (@{VISCHER:1949ty}). 

In 1945, Hotchkiss ultimately proved Johnson and Coghill correct when he isolated 5-methylcytosine from nucleic acid prepared from cow thymus (@{HOTCHKISS:1948va}). Using paper chromatography, Hodgkiss demonstrated that methylated cytosine existed and was distinct from conventional cytosine and uracil.

[^tuberculosis]: _Mycobacterium tuberculosis_ was then known as _Tubercle bacillus_.

In the next section I focus on DNA methylation in mammals, in particular humans, but include brief descriptions of DNA methylation in other organisms. 

### DNA methylation in mammals
DNA methylation is essential for normal development in mammals, such as human and mouse (__REFERENCE__). My thesis focuses on DNA methylation in humans (_Homo sapiens_) but much of the material in this section is true for other mammals.

The typical site of DNA methylation is at the C5 carbon position of a cytosine nucleotide. This modified cytosine is called _5-Methylcytosine_ or _5mC_. [__FIGURE__](http://en.wikipedia.org/wiki/File:5-Methylcytosine.svg) shows the structure of 5mC. In mammalian genomes, most cytosines are unmethylated except for those at the _CpG dinucleotide_. I explain the importance of CpG dinucleotides in the next section. 

#### CpG dinucleotides {-}

A _CpG dinucleotide_, or more simply a _CpG_, is a cytosine followed by a guanine in the linear DNA sequence (__FIGURE__). The "p" stands for the phosphate backbone of DNA and some authors omit it in favour of simply calling it CG methylation. Approximately 60-80% of CpGs are methylated in mammals, meaning that the cytosine in the dinculeotide is a 5-Methylcytosine (__REFERENCE__). I will use the abbreviation _mCpG_ to refer to a methylated CpG, _uCpG_ to refer to an unmethylated CpG and simply _CpG_ when the methylation state is irrelevant.

CpGs are underrepresented in the human genome. The GC-content of the human genome, which is defined as the percentage of nucleotides that are either guanines or cytosines, is approximately __WHAT %__. Under a random sampling model of nucleotides we would expect approximately __WHAT %__ of the human genome to be comprised of CpGs. In fact, only __WHAT %__ of the dinucleotides are CpGs. __FIGURE__ shows the frequency of each dinucleotide in the human reference genome (hg19)[^cpg_frequency].

[^cpg_frequency]: See __APPENDIX__ for the `R` code used to compute these values for the human reference genome (build hg19).

One reason for this relative scarcity of CpGs is that methylated cytosines can spontaneously deaminate to thymines ([@Scarano:1967ue]). Thus, over time, many mCpGs will become TpGs, leading to a genome-wide reduction in the proportion of CpGs and a genome-wide increase in the proportion of TpGs (see __FIGURE__).

There are many other evolutionary pressures on the distribution of nucleotides in a genome. One effect of this is that the distribution of CpGs is far from uniform. I discuss this distribution in the next section.

#### CpG islands and other sandy metaphors {-}

One way to explore the distribution of CpGs in the human genome is to look at the distribution of distances from one CpG to the next. __FIGURE__ is a kernel density plot of the distances between neighbouring CpGs  and __FIGURE__ is the empirical cumulative distribution function of the distances between neighbouring CpGs for the human reference genome (hg19)[^cpg_distance].

[^cpg_distance]: See __APPENDIX__ for the `R` code used to generate these figures.

As can be seen from __FIGURE__, there are a subset of CpGs that are within __HOW MANY__ bases of the next CpG.

The 20-40% of CpGs that are methylated in mammalian genomes are mostly found in clusters called _CpG islands_ (_CGIs_). Again, there are evolutionary arguments for why these CGIs form (__REFERENCE__).

CGIs are important regulatory elements in the genome and are where most differences in DNA methylation between different cell types are found (__REFERENCE__). [@GardinerGarden:1987we] gave the classical definition of a CGI:

> For the purpose of this survey, regions of DNA with a moving average of %G+C over 50 and Obs/Exp CpG over 0.6 have been classed as CpG-rich regions. CpG-rich regions over 200 bp in length are unlikely to have occurred by chance alone, so, as a working definition, have been labelled as CpG islands.

Other authors have since provided alternative definitions (__REFERENCES__ Adrian Bird, Rob Klose, Hao Wu, etc.)

[@Wu:2010do] developed a Hidden Markov Model to define CGIs. In my thesis I have used the predicted CGIs from this model. These predicted CGIs are available from [http://rafalab.jhsph.edu/CGI/](http://rafalab.jhsph.edu/CGI/).

The "sandy/beachy" metaphors have been continued (i.e. stretched to breaking point) with various authors defining CGI shores, CGI shelves, CpG deserts and CpGs in the "open sea". CGI shores (@{Irizarry:2008hga}) are defined as regions within 2kb of CGIs. They regions have been demonstrated to have an increased variability of CpG methylation (@{Irizarry:2008hga}). CGI shelves (__WHO__) are defined as __WHAT__ and __WHY ARE THEY USEFUL__. CpG deserts are __WHAT__ and the "open sea" is __WHAT AND HOW DOES THIS DIFFER FROM CpG DESERTS__

#### Non-CpG methylation {-}

In humans, cytosine methylation in most cell types is found almost exclusively at CpGs. However, there are certain cell types with widespread non-CpG methylation. Non-CpG methylation is often classified as CHH methylation or CHG methylation, where H stands for any nucleotide except G (__SOURCE IUPAC__). The rule-of-thumb is that non-CpG methylation is rare in somatic cells but common in pluripotent cells. Of course, there are exceptions to every rule, especially in biology.

To give a few examples, [@Lister:2009hy] found that in a somatic cell line (_fibroblasts_ or skin cells) $99.98\%$ of methylcytosines occured at CpGs whereas as in an embryonic stem cell line $24.5\%$ of methylcytosines occured in a non-CpG context. A subsequent paper from the same group ([@Lister:2011kg]) extended this result when they reported that, more generally, non-CpG methylation accounts for $20-30\%$ of methylcytosines in _pluripotent_ cell lines (which includes embryonic stem cell lines along with induced pluripotent stem cell lines, see __SECTION___). An exception to the rule is provided by [@{Lister:2013et}], where they found that neurons, a somatic and not a pluripotent type of cell, also have non-CpG methylation, albeit at a lower level ($1.3-1.5\%$ of all non-CpG cytosines were methylated).

Non-CpG methylation is less well studied than CpG methylation. This is partly due to the design of popular assays for studying cytosine methylation; for example the popular Illumina 27k and 450k beadchips (see __SECTION__) only measure cytosine methylation at CpGs __CHECK__. However, recent technological advances (see __SECTION__) mean that cytosine methylation can be routinely assayed regardless of the sequence context. The biological role of non-CpG methylation is less well understood than that of CpG methylation (__SOURCE__).

#### {-} 5hmC, 5acC, etc.

Methylation is not the only chemical modification of cytosine nucleobase, although it is by far the most common. Further modifications, listed approximately from most frequent to least frequent, are 5-hydroxymethylcytosine (5hmC), 5-formylcytosine (5fC) and 5-carboxylcytosine (5caC). 

The biological significance of these marks is still being determined, in part because the assays for studying these are still in development (__SOURCE__) and because their relatively scarcity means that experiments to detect these modifications are more difficult and expensive. 

One genome-wide study of 5hmC found that less than $ 1\%$ of all assayed cytosines in mouse fetal cortex and adult cortex cells were hydroxymethylated (@{Lister:2013et}). Most of the 5hmC was detected in the CpG context and, although the genome-wide level of 5hmC was low, the authors reported significant levels of 5hmC at particular cytosines in the genome.

__Kriaucionis and Heintz, 2009; Tahiliani et al., 2009__ discovered that the TET enzymes can convert 5mC to 5hmC, 5hmC to 5fC and 5fC to hcaC. This suggested a role for 5hmC, 5fC and 5caC in the process of removing 5mC marks, as is discussed in the next section.

#### Readers, writers and erasers {-}

#### Function of DNA methylation {-}

#### DNA methylation in other organisms and kingdoms {-}

##### Insects {-}

* Honeybees
* Fruit fly? DNA methylation's existance in other organisms, such as the common fruit fly (_Drosophila Melanogaster_), is an area of contention although recent studies suggest it is present albeit at very small quantities and only in very limited contexts [@Takayama:2014jp].
* General

##### Plants {-}

* Arabadopsis
* Non-CG methylation

##### Fungi {-}

##### Bacteria {-}

##### Archaea {-}


## Assays for studying DNA methylation

## Bioinformatics analysis of bisulfite-sequencing data

## Outline of thesis

## General TODOs
* Nucleotide, residue, (nucleo)base; what's the difference and which term should I use?
	* A nucleotide is composed of a nucleobase (also termed a nitrogenous base), a five-carbon sugar (either ribose or 2-deoxyribose), and one or more phosphate groups \[Coghill, Anne M.; Garson, Lorrin R., ed. (2006). The ACS style guide: effective communication of scientific information (3rd ed.). Washington, D.C.: American Chemical Society. p. 244. ISBN 978-0-8412-3999-9.)\]
* Should I give more detail on Rafa's CGI HMM?
* Check how deep sub-sections should go and fix.
* Define eukaryotes, prokaryotes and viruses
