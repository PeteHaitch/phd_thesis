# Thesis outline

## Thesis title

Statistical analysis of high-throughput assays for studying DNA methylation

## Research questions

* What DNA methylation is
* What bisulfite-sequencing is and how people use it
* Why statistics is needed
* What questions I have addressed
* How my work is useful

# Chapter outlines

For each chapter I outline my plans and add checkboxes for what is and is not drafted.

## Preamble

Do last.

## Introduction

- [x] Chapter overview
- [ ] DNA
  - [ ] basic structure, mitotic replication, DNA -> nucleosomes -> chromatin -> chromosomes -> genome, DNA -> RNA -> protein.
- [ ] Epigenetics
- [x] DNA methylation
  - Mostly done.
- [x] Assays for studying DNA methylation
  - Mostly done.
- [ ] WGBS
- [ ] Outline of thesis

## Bioinformatics analysis of bisulfite-sequencing data

This chapter is mostly complete. Need to re-work the intro a little to be more general than differential methylation experiments.

Put this in context with Krueger et al.'s Nature Methods paper and any other papers that explain the "pipeline" for a bioinformatics analysis of BS-seq data.

- [x] Chapter overview
- [x] Data quality control checks
- [x] Read mapping and post-processing of mapped reads
- [x] Methylation calling
- [ ] Downstream analyses

## A statistical framework for analysing bisulfite-sequencing data

I should emphasise that this chapter is also useful because this framework is lacking, which makes it difficult to compare methods.

This chapter should focus on the statistical framework and not the bioinformatics problems (these should go in the previous chapter) and so this chapter requires some re-working.

- [x] Chapter overview
- [ ] One sample
- [ ] _n_ samples
- [ ] Parameter estimation
- [ ] Statistical properties of $\beta$
- [ ] Downstream analyses

## Datasets used in thesis

This chapter is a good candidate to write during the bootcamp.


## Differential methylation in Agouti viable yellow mice

This chapter is a good candidate to write during the bootcamp.

- [ ] Chapter overview
- [ ] Introduction
- [ ] Background + literature review
- [ ] Methodology
- [ ] Discussion of findings
- [ ] Conclusion

## Co-methylation

I want to combine the "Co-methylation review" and current "Co-methylation" chapter into a single chapter.

- [ ] Chapter overview
- [ ] Introduction
- [ ] Background
- [ ] Literature review
  - [ ] Correlations of aggregate methylation
  - [ ] Within-fragment co-methylation
  - [ ] Between sample-comethylation
- [ ] Methodology
  - [ ] Correlations of $\beta$-values
  - [ ] Within-fragme co-methylation
- [ ] Implementation
- [ ] Case studies
  - [ ] Lister data
  - [ ] EPISCOPE
  - [ ] Ziller
- [ ] Discussion of findings
  - [ ] Biology
  - [ ] Limitations of using 2-tuples
  - [ ] Accounting for and exploiting co-methylation in analyses of differential methylation
- [ ] Conclusions

## Co-methylation review

This chapter is mostly complete. Just need to review some old papers and decide whether this should really be a stand-alone chapter or a part of the next chapter.

- [x] Chapter overview
- [x] What is co-methylation and why study it
- [ ] Correlations of aggregate methylation
    - Mostly done, still need to review Landan et al.
- [x] Within-fragment co-methylation
- [ ] Between sample co-methylation
  - Still to review Capra and Kotska (2014) and Liu et al. (2014)
- [ ] Exploiting co-methylation in analyses of differential methylation

## Co-methylation

- [x] Chapter overview
- [ ] Correlations of beta-values
- [ ] Within-fragment co-methylation
  - Partially written
- [ ] Using higher-order m-tuples
- [ ] Results of co-methylation analyses

## A simulation model of DNA methylation data

This chapter is a good candidate to write during the bootcamp.

- [ ] Chapter overview
- [ ] Introduction
- [ ] Background + literature review
- [ ] Methodology
- [ ] Implementation
- [ ] Case study
- [ ] Discussion of findings
- [ ] Conclusion

## Concluding remarks

## Bibliography

## Appendices
