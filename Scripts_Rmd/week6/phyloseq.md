Phyloseq
================
Nicholas Baetge
10/31/2020

# Intro

We explore the processed ACIDD 16S sequences using
[phyloseq](https://joey711.github.io/phyloseq/)

# Install phyloseq

``` r
# BiocManager::install("phyloseq")
```

``` r
library(tidyverse) 
library(phyloseq)
library(RColorBrewer)
```

# Import Data

``` r
count.tab <- read_rds("~/GITHUB/eemb144l/Input_Data/week6/seqtab-nochimtaxa.rds") #table of counts for each sequence in each sample
tax.tab <- read_rds("~/GITHUB/eemb144l/Input_Data/week6/taxa.rds") #table that matches ASV to sequence
sample.tab <- read_rds("~/GITHUB/eemb144l/Output_Data/week4/ACIDD_Exp_Processed_DOC_BGE.rds") %>% 
  drop_na(DNA_SampleID) %>% 
  column_to_rownames(var = "DNA_SampleID") 
```

# Phyloseq Object

We need to create a phyloseq object that merges all three datasets.
Sometimes this doesn’t work beacuse of the format of the data files.
Make sure all the sample names between the sampleinfo.txt and
seqtab-nochimtaxa.txt are the same

``` r
OTU = otu_table(count.tab, taxa_are_rows = TRUE) 
TAX = tax_table(tax.tab)
SAM = sample_data(sample.tab)
ps = phyloseq(OTU,TAX,SAM) 
```

# Filter sequences

We will filter out chloroplasts and mitochondria, because we only
intended to amplify bacterial sequences. It’s good to check you don’t
have anything lurking in the taxonomy table.

``` r
sub_ps <- ps %>%
  # subset_samples(Experiment == "ASH172") %>%  #use this function if you want to only include some subset of your sample set in the subsequent analysis
  subset_taxa(Family  != "mitochondria" & Order  != "Chloroplast")
```

# Sample Summary

As a first analysis, we will look at the distribution of read counts
from our
samples

<img src="phyloseq_files/figure-gfm/unnamed-chunk-6-1.png" style="display: block; margin: auto;" />

``` r
# mean, max and min of sample read counts
summary(sample_sum_df)
```

    ##       sum      
    ##  Min.   :6048  
    ##  1st Qu.:7152  
    ##  Median :7633  
    ##  Mean   :7687  
    ##  3rd Qu.:8454  
    ##  Max.   :9253

# Beta Diversity

Beta diversity involves calculating metrics such as distances or
dissimilarities based on pairwise comparisons of samples – they don’t
exist for a single sample, but rather only as metrics that relate
samples to each other. i.e. beta diversity = patterns in community
structure between samples

Since differences in sampling depths between samples can influence
distance/dissimilarity metrics, we first need to somehow normalize the
read depth across our samples.

## Subsample

We will rarefy (random subsample with replacement) the read depth of the
samples first (scale to the smallest library size).

[Case for not
subsampling](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003531)

[Response blog for
subsampling](https://www.polarmicrobes.org/how-i-learned-to-stop-worrying-and-love-subsampling-rarifying/)

Read depth is an artifact of a machine made by a company in San Diego,
not anything about your samples or their biology. It is totally
artifactual, and controlling for artifacts is critical in science.
Subsampling randomly is the simplest way to control for this, and the
question is whether this is the “best” way of controlling for it. See
links above for alternative arguments about what the best way of
controlling for this artifact is.

A strong reason to subsample is to standardize effort. The bottom line
is that in all experimental design you should not be comparing things to
which you devote different effort in resolution. For instance, you don’t
sample one site once a week and another once a month if you want to
compare the dynamics between the sites. You standardize effort.

With that said, the bigger your differential in mean (or median) read
depth (reads/sample) between pre- and post-subsampling, the greater the
“effect” on beta diversity.

Examples:

  - means reads before = 40k, mean reads after = 1k, big effect.
  - mean reads before = 40k, mean reads after = 20k, small effect.
  - mean reads before = 2k, mean reads after = 1k, small effect.

We will subsample to the minimum read depth of all samples and not
subsample. We’ll then compare the mean reads pre- and post-subsampling
and also compare beta diversity
patterns

``` r
ps_min <-  rarefy_even_depth(sub_ps, sample.size = min(sample_sums(sub_ps)))
```

    ## You set `rngseed` to FALSE. Make sure you've set & recorded
    ##  the random seed of your session for reproducibility.
    ## See `?set.seed`

    ## ...

``` r
mean(sample_sums(sub_ps)) #7686
```

    ## [1] 7686.875

``` r
mean(sample_sums(ps_min)) #6048 this is also the same as min(sample_sums(sub)ps) 
```

    ## [1] 6048

Based on the mean reads pre- and post-subsampling, subsampling here
shouldn’t have a major effect on our beta diversity patterns

## NMDS

One of the best exploratory analyses for amplicon data is unconstrained
ordinations. Here we will look at non-metric multidimensional scaling
(NMDS) ordinations of our full community samples. For NMDS plots it’s
important to set a seed since the starting positions of samples in the
alogrithm is random.

``` r
set.seed(1)
# Ordinate
nmds <- ordinate(sub_ps, method = "NMDS",  distance = "bray") # stress = 0.04
```

    ## Square root transformation
    ## Wisconsin double standardization
    ## Run 0 stress 0.04221621 
    ## Run 1 stress 0.17179 
    ## Run 2 stress 0.04221623 
    ## ... Procrustes: rmse 8.273279e-05  max resid 0.0001418332 
    ## ... Similar to previous best
    ## Run 3 stress 0.0422162 
    ## ... New best solution
    ## ... Procrustes: rmse 2.301456e-05  max resid 3.890594e-05 
    ## ... Similar to previous best
    ## Run 4 stress 0.0422162 
    ## ... New best solution
    ## ... Procrustes: rmse 2.836402e-05  max resid 4.43833e-05 
    ## ... Similar to previous best
    ## Run 5 stress 0.04221622 
    ## ... Procrustes: rmse 4.5349e-05  max resid 7.537188e-05 
    ## ... Similar to previous best
    ## Run 6 stress 0.0422162 
    ## ... Procrustes: rmse 3.680198e-05  max resid 5.67041e-05 
    ## ... Similar to previous best
    ## Run 7 stress 0.04221622 
    ## ... Procrustes: rmse 2.862554e-05  max resid 5.607923e-05 
    ## ... Similar to previous best
    ## Run 8 stress 0.04221621 
    ## ... Procrustes: rmse 3.203284e-05  max resid 5.156382e-05 
    ## ... Similar to previous best
    ## Run 9 stress 0.04221626 
    ## ... Procrustes: rmse 7.681552e-05  max resid 0.0001262773 
    ## ... Similar to previous best
    ## Run 10 stress 0.04221621 
    ## ... Procrustes: rmse 5.109692e-05  max resid 8.127412e-05 
    ## ... Similar to previous best
    ## Run 11 stress 0.04221621 
    ## ... Procrustes: rmse 4.87695e-05  max resid 7.580248e-05 
    ## ... Similar to previous best
    ## Run 12 stress 0.227169 
    ## Run 13 stress 0.04221619 
    ## ... New best solution
    ## ... Procrustes: rmse 1.32125e-06  max resid 2.475487e-06 
    ## ... Similar to previous best
    ## Run 14 stress 0.0422162 
    ## ... Procrustes: rmse 2.608597e-05  max resid 4.349033e-05 
    ## ... Similar to previous best
    ## Run 15 stress 0.04221621 
    ## ... Procrustes: rmse 4.857337e-05  max resid 7.662114e-05 
    ## ... Similar to previous best
    ## Run 16 stress 0.0422162 
    ## ... Procrustes: rmse 2.000833e-05  max resid 3.188521e-05 
    ## ... Similar to previous best
    ## Run 17 stress 0.0422162 
    ## ... Procrustes: rmse 3.929749e-05  max resid 6.287313e-05 
    ## ... Similar to previous best
    ## Run 18 stress 0.0422162 
    ## ... Procrustes: rmse 7.720879e-06  max resid 1.303811e-05 
    ## ... Similar to previous best
    ## Run 19 stress 0.0422162 
    ## ... Procrustes: rmse 2.031009e-05  max resid 3.256571e-05 
    ## ... Similar to previous best
    ## Run 20 stress 0.0422162 
    ## ... Procrustes: rmse 3.396999e-05  max resid 5.407073e-05 
    ## ... Similar to previous best
    ## *** Solution reached

``` r
set.seed(1)
# Ordinate
nmds_min <- ordinate(ps_min, method = "NMDS",  distance = "bray") # stress = 0.04
```

    ## Square root transformation
    ## Wisconsin double standardization
    ## Run 0 stress 0.0450746 
    ## Run 1 stress 0.1672519 
    ## Run 2 stress 0.04507461 
    ## ... Procrustes: rmse 3.382258e-05  max resid 6.200966e-05 
    ## ... Similar to previous best
    ## Run 3 stress 0.2118775 
    ## Run 4 stress 0.1662421 
    ## Run 5 stress 0.0450746 
    ## ... Procrustes: rmse 1.178065e-05  max resid 1.994252e-05 
    ## ... Similar to previous best
    ## Run 6 stress 0.0450746 
    ## ... Procrustes: rmse 7.886302e-06  max resid 1.464252e-05 
    ## ... Similar to previous best
    ## Run 7 stress 0.0450746 
    ## ... Procrustes: rmse 1.267304e-05  max resid 2.314787e-05 
    ## ... Similar to previous best
    ## Run 8 stress 0.0450746 
    ## ... New best solution
    ## ... Procrustes: rmse 6.722885e-07  max resid 1.745986e-06 
    ## ... Similar to previous best
    ## Run 9 stress 0.0450746 
    ## ... Procrustes: rmse 5.929639e-06  max resid 9.927679e-06 
    ## ... Similar to previous best
    ## Run 10 stress 0.0450746 
    ## ... New best solution
    ## ... Procrustes: rmse 3.0888e-06  max resid 5.47317e-06 
    ## ... Similar to previous best
    ## Run 11 stress 0.0450746 
    ## ... Procrustes: rmse 1.169706e-06  max resid 2.451393e-06 
    ## ... Similar to previous best
    ## Run 12 stress 0.1985986 
    ## Run 13 stress 0.0450746 
    ## ... Procrustes: rmse 7.881676e-06  max resid 1.398886e-05 
    ## ... Similar to previous best
    ## Run 14 stress 0.0450746 
    ## ... Procrustes: rmse 4.53905e-06  max resid 7.800073e-06 
    ## ... Similar to previous best
    ## Run 15 stress 0.0450746 
    ## ... Procrustes: rmse 3.839374e-06  max resid 9.302879e-06 
    ## ... Similar to previous best
    ## Run 16 stress 0.0450746 
    ## ... New best solution
    ## ... Procrustes: rmse 1.436219e-06  max resid 2.63293e-06 
    ## ... Similar to previous best
    ## Run 17 stress 0.0450746 
    ## ... Procrustes: rmse 7.539374e-06  max resid 1.370446e-05 
    ## ... Similar to previous best
    ## Run 18 stress 0.0450746 
    ## ... Procrustes: rmse 4.447174e-06  max resid 8.32403e-06 
    ## ... Similar to previous best
    ## Run 19 stress 0.0450746 
    ## ... Procrustes: rmse 5.592817e-06  max resid 1.012632e-05 
    ## ... Similar to previous best
    ## Run 20 stress 0.0450746 
    ## ... Procrustes: rmse 2.189923e-06  max resid 4.387174e-06 
    ## ... Similar to previous best
    ## *** Solution reached

<img src="phyloseq_files/figure-gfm/unnamed-chunk-11-1.png" style="display: block; margin: auto;" />

<img src="phyloseq_files/figure-gfm/unnamed-chunk-12-1.png" style="display: block; margin: auto;" />

NMDS plots attempt to show ordinal distances between samples as
accurately as possible in two dimensions. It is important to report the
stress of these plots, because a high stress value means that the
algorithm had a hard time representing the distances between samples in
2 dimensions. The stress of this plot was good - it was .04 (generally
anything below .2 is considered acceptable).

Subsampling doesn’t appear to affect the patterns we see in beta
diversity, so moving forward, we will focus on the subsampled dataset.

# Alpha Diversity

Estimating alpha diversity of microbial communities is
[problematic](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC93182/) no
matter what you do.

We are going to calculate the Chao1 index for richness and the Shannon
diversity index.

**it is important to note that the alpha diversity values are not
interpretable as “real” numbers of anything (due to the nature of
amplicon data), but they can still be useful as relative metrics of
comparison. If Chao1 richness goes up, but Shannon diversity goes down,
it indicates that the sample may have more ASVs but is dominated by a
few of them.**

We will use the subsampled library, which retains estimates of the
species abundance of the real population while standardizing sampling
effort.

[subsampling and alpha diversity
paper](https://www.frontiersin.org/articles/10.3389/fmicb.2019.02407/full)

[Chao1: nonparametric estimation of minimum community
richness](https://www.jstor.org/stable/4615964?seq=1#metadata_info_tab_contents)

``` r
richness <- estimate_richness(ps_min, measures = c("Chao1", "Shannon")) %>% 
  rownames_to_column(., var = "DNA_ID") %>% 
  mutate_at(vars(DNA_ID), str_replace_all, pattern = "171.", "171-") %>% 
   mutate_at(vars(DNA_ID), str_replace_all, pattern = "172.", "172-")
```

Let’s add the sample metadata into this
dataframe

``` r
alphadiv <- left_join(richness, sample.tab %>% rownames_to_column(., var = "DNA_ID")) 
```

    ## Joining, by = "DNA_ID"

<img src="phyloseq_files/figure-gfm/unnamed-chunk-15-1.png" style="display: block; margin: auto;" />

Boxes represent the 1.5 interquartile range, with the internal solid
line representing the median. Circles represent data points. p-values
are reported the non-parametric two sample Wilcoxon test, which tests
whether the means between two groups are equal (ns: p \> 0.05, \* : p≤
0.05, \*\* : p ≤ 0.01).

Difference in the alpha diversity indexes among conditions were tested
using pairwise Wilcoxon tests; p \< 0.05 was considered the threshold
significance for a difference between conditions.

From this plot we can see within the treatments that the richness (via
Chao index) of our samples significantly changed, while overall
diversity (via Shannon index) did not change. This suggest that while
richness decreased in both the control and ash leachate treatments, the
eveness was similar between the initial and final
conditions.

<img src="phyloseq_files/figure-gfm/unnamed-chunk-16-1.png" style="display: block; margin: auto;" />

From this plot we can see between the treatments that the richness of
the control samples were higher at the initial condition than the ash
leachate, suggesting that there may have been some quality control
issues as we would expect the initial samples to all have the same
richness. By timepoint 6, it looks like the richness was about the same
between the control and the ash leachate. Overall diversity was similar
between the treatments at the initial condition, but not by the end of
the experiment. The ash leachate samples at timepoint 6 may have been
less even.

**In sum, we observe that richness similarly decreased within the
treatments over time, but overall diversity was lower in the ash
leachate experiments, suggesting that those communities became less rich
and less even relative to the control. i.e. relative to the control, ash
leachate samples had less ASVs and were dominated by fewer of them.**

# Who??

Which taxa were important? Which taxa were contributing to the change in
community compositon?

**Note: Recovered 16S rRNA gene copy numbers do not equal organism
abundance.**

That said, we can generate a heat map of our samples showing us how the
relative abundance of different taxonomic groups change…potentially
giving us a visual of which taxa are most important to the alpha and
beta diversity patterns we observed. First, we’re going to generate a
custom table that will be easier to work with than a phyloseq object.

## Generate relative abundances

Our data currently shows number gene copies recovered, so we’ll convert
to percentages (relative abundances)

``` r
ps_std <- transform_sample_counts(ps_min, function(x) x/sum(x))
#extract the relative abundance table and coerce into dataframe
ps_std.tab <- as(otu_table(ps_std), "matrix")
ps_std.df = as.data.frame(ps_std.tab) 
```

## Make table

``` r
#first coerce the taxa table into a data frame
tax.df <-  as.data.frame(tax.tab) 
#then combine the data frames
custom.tab <- tax.df %>% 
  rownames_to_column(., var = "asv") %>% 
  left_join(., ps_std.df %>% rownames_to_column(., var = "asv")) %>% 
  #create a new index of that combines the  class, order, family, and genus values, you can play around here!!
  mutate(#pcofg = paste(Phylum, "_", Class, "_", Order,"_", Family, "_", Genus),
         # pcof = paste(Phylum, "_", Class, "_", Order,"_", Family,),
         pco = paste(Phylum, "_", Class, "_", Order)) %>% 
  select(-c(asv:Genus)) %>% 
  # select(pcof,everything()) %>% 
  # group_by(pcof) %>% 
  select(pco,everything()) %>% 
  group_by(pco) %>% 
  #here we are combining the relative abundances based on our grouping
  summarise_at(vars(contains(c("ASH171", "ASH172"))), sum, na.rm = T) %>% 
  ungroup()
```

    ## Joining, by = "asv"

``` r
#save the row names and then make them into the column names
colnames <- custom.tab[,1] 

#transpose the dataframe so we can merge with the sample info table
t_custom.tab <-  as.data.frame(t(custom.tab[,-1]))
# colnames(t_custom.tab) <- colnames$pcof
colnames(t_custom.tab) <- colnames$pco

#merge
sweet.tab <- t_custom.tab %>% 
  rownames_to_column(., var = "sample") %>% 
  left_join(., sample.tab %>% rownames_to_column(., var = "sample") %>% select(sample, Experiment, Location, Bottle, Treatment, Timepoint, days, cells)) %>% 
  select(sample, Experiment:cells, everything())
```

    ## Joining, by = "sample"

``` r
relabund <- sweet.tab %>% 
  select(-c(sample:cells)) %>% 
  #remove groups that are completely absent
  .[ , colSums(.) > 0] %>% 
  #arrange by biggest contributors
  .[, order(colSums(-.))] %>% 
  bind_cols(sweet.tab %>% select(sample:cells), .)
```

## Heatmap

<img src="phyloseq_files/figure-gfm/unnamed-chunk-19-1.png" style="display: block; margin: auto;" />

From this kind of plot, you can get an idea of who responded to the ash
leachate…that is, who uniquely increased in relative abundance in the
ash leachate treatment that did not in the control? what about decrease?

It looks like ASVs belonging to Oceanospirialles, Campylobacterales, and
Alteromonadales all increased in relative abundance in the ash leachate
experiments, but did not in the control. It looks any decreases in ASVs
were similar between the treatments\! Pretty cool data\!\!

Everything shown here is just a snapshot of what you can look at with
your community composition data. There are many other resources you can
use to get ideas of how to look at different aspects of your data,
including the [phyloseq tutorial](https://joey711.github.io/phyloseq/)
and [happy belly bioinformatics](https://astrobiomike.github.io). It’s
up to you and your
questions\!\!

# Save and knit

``` r
saveRDS(sweet.tab, "~/GITHUB/144l_students/Output_Data/week6/Custom_ASV_Table.rds")
saveRDS(sub_ps, "~/GITHUB/144l_students/Output_Data/week6/phyloseq_obj.rds")
saveRDS(ps_min, "~/GITHUB/144l_students/Output_Data/week6/subsampled_phyloseq_obj.rds")
saveRDS(alphadiv, "~/GITHUB/144l_students/Output_Data/week6/alphadiv.rds")
```
