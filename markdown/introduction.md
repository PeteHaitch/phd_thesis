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

#### Genes and DNA transcription {-}
DNA would not be quite so interesting if all it did was self-replicate. __IT ALSO ENCODES INFO FOR NEW MOLECULES__
* DNA is _transcribed_ to produce _RNA_ (ribonucleicacid). 
* RNA is _translated_ to produce _proteins_. 




#### The Central Dogma of molecular biology {-)
Francis Crick, who along with James Watson and Rosalind Franklin discovered the double helix structure of DNA, proposed the Central Dogma of molecular biology [__REFERENCE__(http://sandwalk.blogspot.com.au/2007/01/central-dogma-of-molecular-biology.html)]:
> ... once (sequential) information has passed into protein it cannot get out again (F.H.C. Crick, 1958)

Crick re-stated this in 1970:
> The central dogma of molecular biology deals with the detailed residue-by-residue transfer of sequential information. It states that such information cannot be transferred from protein to either protein or nucleic acid. (F.H.C. Crick, 1970)

Crick further added (1970):
> It is not the same, as is commonly assumed, as the sequence hypothesis, which was clearly distinguished from it in the same article (Crick, 1958). In particular, the sequence hypothesis was a positive statement, saying that the (overall) transfer nucleic acid → protein did exist, whereas the central dogma was a negative statement saying that transfers from protein did not exist.

Lewin’s GENES VIII defines the Central Dogma of Molecular Biology as,
> The central dogma states that information in nucleic acid can be perpetuated or transferred but the transfer of information into protein is irreversible. (B. Lewin, 2004)

#### Genomes {-}
* Most DNA in humans is non-coding
* Much of it is junk
* junk != non-coding
* DNA -> chromosomes -> genomes
* Some data are mapped to hg18 and some to hg19.

### Genetics

### Epigenetics
* Waddington and others' definitions
* 

## DNA methylation
DNA methylation is a chemical modification of DNA that can impart information on top of the DNA's genetic code. This notion of DNA methylation "encoding information on top of DNA" gives rise to its definition as an _epigenetic_ mark, _epi_ being derived from the Greek word for "upon, near to or in addition"__INSERT OXFORD DICTIONARY REFERENCE__. The broader role of epigenetics is discussed below in the section on [Epigenetics].

DNA methylation is heritable during mitotic cell division, which means that it is faithfully copied across to the daughter cell during cell division[^mitotic].

[^mitotic]: In practice this copying is not as faithful as, say, the copying of DNA from the parent to the daughter strand. Furthermore, the faithfulness of this copying will be different in different conditions, such as in a healthy, well-differentiated liver cell compared to a cancerous liver cell. Nevertheless, the copying of DNA methylation is faithful enough for most biologists to consider it as a mitotically heritable mark, most of the time (__TODO__: Is this just true of 5mC in mammals or more generally?)

So pervasive is DNA methylation that [@Lister:2009bo] dubbed it the "fifth base" of the DNA code. In fact, what is referred to as DNA methylation is typically referring to methylation of the cytosine nucleotide. Cytosine methylation is by far the most common form of DNA methylation in the animal and plant kingdoms, however, I will continue to use the terms _DNA methylation_ and _cytosine methylation_ interchangeably as is standard in the literature. In the next section I focus on DNA methylation in mammals, in particular humans, but include brief descriptions of DNA methylation in other organisms. 

### DNA methylation in mammals
DNA methylation is essential for normal development in mammals, such as human and mouse (__REFERENCE__). My research has focused on DNA methylation in modern humans (_Homo sapiens sapiens_) but much of the material in this section is true for other mammals.

The typical site of DNA methylation is at the C5 carbon position of a cytosine nucleotide. This modified cytosine is called _5-Methylcytosine_ or _5mC_. [__FIGURE__](http://en.wikipedia.org/wiki/File:5-Methylcytosine.svg) shows the structure of 5mC. Most cytosines are unmethylated except for those at the _CpG dinucleotide_. I explain the importance of CpG dinucleotides in the next section. 

#### CpG dinucleotides {-}
A _CpG dinucleotide_, or more simply a _CpG_, is a cytosine followed by a guanine in the linear DNA sequence (__FIGURE__). Approximately 60-80% of CpGs are methylated in mammals, meaning that the cytosine in the dinculeotide is a 5-Methylcytosine (__REFERENCE__). I will use the abbreviation _mCpG_ to specifically refer to a methylated CpG.

CpGs are underrepresented in the human genome. The GC-content of the human genome, which is defined as the percentage of nucleotides that are either guanines or cytosines, is approximately __WHAT %__. Under a random sampling model of nucleotides we would expect approximately __WHAT %__ of the human genome to be comprised of CpGs. In fact, only __WHAT %__ of the dinucleotides are CpGs. __FIGURE__ shows the frequency of each dinucleotide in the human reference genome (hg19)[^cpg_frequency].

[^cpg_frequency]: See __APPENDIX__ for the `R` code used to compute these values for the human reference genome (build hg19).

One reason for this relative scarcity of CpGs is that methylated cytosines can spontaneously deaminate to thymines [@Scarano:1967ue]. Thus, over time, many mCpGs will become TpGs, leading to a genome-wide reduction in the proportion of CpGs and a genome-wide increase in the proportion of TpGs (see __FIGURE__).

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

__TODO__: CpG shores, shelves, deserts, etc.


#### Non-CG methylation {-}

#### DNA methyltransferases {-}

#### Function of DNA methylation {-}

#### DNA methylation in other organisms and kingdoms

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
