Tonight and on Saturday I want to work on my simulation chapter. On Sunday I will work on the Emma Whitelaw chapter. I want to get as much of the "non-results" writing done as possible and to set out the structure for the results so that it is easily incorporated. If time permits, or I get bored/frsutrated, then I will move to the Introduction or the Statistical Analysis of WGBS data chapters. I won't work on the Co-methylation chapter(s) because these are mostly editing and not writing.

## Simulation chapter

These are the main points:

* Realistic simulations are vitally important for methods development.
* Ultimately the real test of a method is how it works on real data but there is no "gold standard" WGBS datasets, hence the need to simulate.
* Current simulation methods suck. Like, really suck.
* Simulation method requirements are:
  1. Realistic
  2. Fast/cheap
  3. Output data in a useful format
  4. Easily set and interpreted parameters
  5. Software implementation
* How the simulation software will be useful
* A simple case study, e.g. 3 vs. 3 DMR experiment
* How to determine if a simulation produces "realistic" data?

## Emma Whitelaw chapter

These are the main points:

* Avy is a cool mouse and what it has taught us about DNA methylation
* WGBS lets us see what else is out there in terms of methylation variability
* The experimental design is unusual, even lacking, which necessitated the development of new methods. But isogenic mice are awesome for looking for DNA methlyation differences.
* It's essentially a one-group experiment where we're looking for outlier sample(s). Moreover, these outliers need to be consistent in order to identify DMRs.
* My role in this project was to get something written quickly that gave reasonable result and supply a software implementation. I started with a table of M + U counts, so relying on Harry's pre-processing of data.
* Basic idea is chi-square test of beta values and then finding runs. Not a new or very sophisticated idea, but pragmatic.
* Do we have any control regions?
* Visualisation of data is __critical__.
* How can this one-group experiment be framed as a GLM?

## Other stuff

I'm a little worried that I only have 4 chapters.
