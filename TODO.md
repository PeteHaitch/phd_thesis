
# Short term

- [x] Avoid first person in thesis? Ask Terry for opinion.
  - Modify if I grimace while reading the sentence.
- [x] Give each chapter a summary.
  - Checked each chapter
- [x] Put table captions at the top, figure captions at the bottom.
  - Check in proofreading
- [x] Add genome-wide line to legend of MH plots
  - __DONE__
- [x] Address Saskia's comments
  - __DONE__
- [x] Address Goknur's comments
  - __DONE__
- [x] Address Terry's comments
  - __DONE__
- [ ] Rename `methsim` plots without the `_1` and `_2` suffixes
- [ ] Plots of beta-correlations separated by strand in a clearer way.
  - Use loess to summarise trend
  - In general, plots of beta-correlations are noisy. Just show trend?
- [ ] "The conditional estimator, works better than the unconditional estimator when the sample size, n, is small [Agresti 2007]." re conditional odds ratio estimators
  - What does "better" mean.
- [ ] Add a representative plot of the within-fragment data for a single sample
- [ ] Fix formatting of $2 \times 2$ contingency tables to be consistent.
- [ ] Write abstract
- [ ] Complete `methsim` chapter.
- [ ] Write concluding remarks
- [ ] Email communications about printing
- [ ] Find submission documents
- [ ] Add blank page between title and abstract
  - Dineika's solution gives the first blank page a number, which isn't what Belinda's thesis looks like.
- [ ] Add analysis folder to github
  - First upload `.Rmd` analysis scripts.
  - Upload <sample>.CG.<m>.tsv.gz Lister, Seisenberger, Ziller samples (m = 1, 2, 2ac) to figshare. All 1-tuples are under Figshare's 250 MB limit but the same is not true of 2-tuples or 2ac-tuples.
- [ ] Fix MethylationTuples::mantelhaen(). Works if copy-pasted into session but not as part of MethylationTuples package.
- [x] Give every figure and table a brief caption.
  - __DONE__
- [ ] Look at some scatterplots of beta-value pairs from simulated data
- "dataset" vs. "data set"
- [ ] \texttt{--all-combinations} prints as '-all-combinations'

# Discussing co-methylation chapter

- 2-tuples vs. pairs. Use pairs instead of 2-tuples (e.g. p153)
- Spatial correlation is not like traditional autocorrelation.
- Plots of beta-correlations separated by strand in a clearer way.
  - Use loess to summarise trend
- "Figures 7.11, 7.12, 7.13 and 7.14 are plots of Spearman correlations of β-values for pairs of CpGs stratified by whether the pair is inside or outside of a CpG island" (p163) __Figure references are wrong; increment by four__
- In general, plots of beta-correlations are noisy. Just show trend?
- p175 "The odds ratio is not perfece, however. For one, the possible values of ψ are highly". The odds ratio is indepndent of the margins (Terry's favourite point about the odds ratio; emphasise).
- top of p176: latex error
- "The conditional estimator, works better than the unconditional estimator when the sample size, n, is small [Agresti 2007]."
  - What does "better" mean.
- p177: s/convert/converge/
- p178: "we analysed a sim. study" -> "we carried out a simulation study"
- When $K$ is known, use the actual value. Improve legends.
- Increase axis labels in simulation figures
- Include figure number when talking about Simulation 3 or 4, etc.
- chrM results are worth commenting on
- chrX with NIL \geq 0 look similar to CGI in that they increase separation.
- Add SE and PE information to plots of within-fragment co-methylation
- Add a representative plot of the within-fragment data for a single sample
- Could remove spearman vs. pearson figures by showing proper "versus" plot and then skipping ahead.
- p197 Something weird has happened.
- Comment on other genomic features that might be interesting to stratify by.
- Comment on how you might reduce heterogeneity by subsampling.



# Long term

- MethylomeParam and simulate,MethylomeParam-method are going to need _another_ re-write :(
  - Slots will be something like MethLevelDT, ComethDT (using Mantel-Haenszel estimates), MethLevelCorDT and all will have to be with respect to the levels of the PartitionedMethylome.
- Once switch from "2" to "1" classes/methods made in methsim, fix other TODOs
- Suggest BPPARAM argument for bsapply(). Suggest bpreplicate() added to BiocParallel.
- Figure out why an assay in a SummarizedExperiment can be a DataFrame but cannot be inserted using `assay<-`.
- Read Hadley's recommendations on capitalising parameter descriptions and implement them for all packages.
- Finish ASD2 report.
- Update excludes for coveralls packages
