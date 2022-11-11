---
bibliography: references.bib
output:
  pdf_document:
    keep_tex: true
csl: asm.csl
geometry: margin=1.0in
header-includes:
 - \usepackage{upgreek}
 - \usepackage{booktabs}
 - \usepackage{longtable}
 - \usepackage{graphicx}
 - \usepackage{array}
 - \usepackage{multirow}
 - \usepackage{wrapfig}
 - \usepackage{float}
 - \usepackage{colortbl}
 - \usepackage{pdflscape}
 - \usepackage{tabu}
 - \usepackage{threeparttable}
 - \usepackage{threeparttablex}
 - \usepackage[normalem]{ulem}
 - \usepackage{makecell}
 - \usepackage{setspace}
 - \doublespacing
 - \usepackage[left]{lineno}
 - \linenumbers
 - \modulolinenumbers
 - \usepackage{helvet} % Helvetica font
 - \renewcommand*\familydefault{\sfdefault} % Use the sans serif version of the font
 - \usepackage[T1]{fontenc}
 - \usepackage[shortcuts]{extdash}
---




# Rarefaction is currently the best approach to control for uneven sampling effort in amplicon sequencing analyses

\vspace{20mm}

**Running title:** Rarefy your data

\vspace{20mm}

Patrick D. Schloss${^\dagger}$

\vspace{40mm}

${\dagger}$ To whom corresponsdence should be addressed:


\href{mailto:pschloss@umich.edu}{pschloss@umich.edu}

Department of Microbiology & Immunology

University of Michigan

Ann Arbor, MI 48109

\vspace{20mm}

**Research article**

\newpage



## Abstract

How to control for uneven sampling effort in amplicon sequencing analyses has become a contentious question in the microbial ecology field since it is common to find as much as 100-fold variation in the number of 16S rRNA gene sequences found across samples in a study. Some have argued that rarefaction is inadmissible because it omits valid data and reduces the statistical power to detect differences between treatment groups. A number of alternative approaches have been developed to normalize and rescale the data that purport to be invariant to the number of observations. I generated community distributions based on 12 published datasets where I was able to assess the sensitivity of different methods to uneven sampling effort. Rarefaction was the only method that could control for variation in uneven sampling effort when measuring commonly used alpha and beta diversity metrics. Next, I compared the false detection rate and power to detect true differences between simulated communities with a known effect size using various alpha and beta diversity metrics. Although all methods of controlling for uneven sampling effort had an acceptable false detection rate when samples were randomly assigned to two treatment groups, rarefaction was consistantly able to contorol for differences in sampling effort when sequencing depth was confounded with treatment group. Finally, the statistical power to detect differences in alpha and beta diversity metrics was consistently the highest when the data were rarefied. These simulations underscore the importance of rarefying the number of sequences across samples in amplicon sequencing analyses.

## Importance

Sequencing of 16S rRNA gene sequences has become a fundamental tool for understanding the diversity of microbial communities and the factors that affect that diversity. Because of technical challenges, it is common to observe wide variation in the number of sequences that are collected from different samples. Because most of the diversity metrics used by microbial ecologists are affected by sampling effort, tools are needed to control for the uneven levels of sampling. This simulation-based analysis shows that despite ongoing controversy, rarefaction is the most robust approach to controlling for uneven sampling effort. I contend that the controversy was started because of confusion over how rarefaction is performed and because methods have been borrowed from other fields that make important assumptions that are not correct with microbial community data.

\newpage

## Introduction

The ability to generate millions sequence reads from 16S rRNA gene amplicons and metagenomic datasets has allowed researchers to multiplex multiple samples on the same sequencing run by pooling separate PCRs that can be deconvoluted later based on index (aka barcode) sequences that embeded into the primer sequence [@Kozich2013; @Caporaso2010]. Unfortunately, it is common to observe variation in the number of sequence reads from each sample vary by as much as 100-fold (e.g, see Figure S1). This occurs because pooling of DNA from multiple PCRs is frought with numerous opportunities for technichal errors to compound leading to a skewed distribution. Aside from developing better methods of pooling DNA, the question of how to control for uneven sampling effort in microbial ecology studies has become controversial.

The practice of rarefaction has been commonly used in ecology for more than 50 years as a tool to control for uneven sampling effort across experimental replicates [@Sanders1968; @Hurlbert1971]. Microbial ecologists have used it to compare 16S rRNA gene sequence data for the past 25 years [@Dunbar1999; @Hughes2001; @Moyer1998]. When data are rarefied, the investigator selects a desired threshold of sampling effort and removes any samples below that threshold. They then randomly select that many sequences, with replacement from each sample. Based on the obseved sequence counts, the researcher can then calcuate alpha diversity metrics including richness and diversity indices or beta divesity metrics such as a Jaccard or Bray-Curtis dissimilarity index. I refer to this single sampling as a subsample; this method is implemented as the sub.sample function in mothur [@Schloss2009] and the rrarefy function in the vegan R package [@Oksanen2022]. Rarefaction repeats the subsampling a large number of times (e.g., 100 or 1,000) and calculates the mean of the alpha or beta diersity metric over those subsamplings; rarefaction is implemented in mothur using the summary.single and dist.shared functions [@Schloss2009] and with the vegan R package using the rarefy or avgdist fuctions [@Oksanen2022]. Rarefaction effectively tells a researcher what an alpha or beta diversity metric would have been for a collection of samples if they were all evenly sampled. Although a closed form equation exists to calcualte the expected richness [@Hurlbert1971], it is computationally easier to empirically calculate rarefied richness and other alpha and beta diversity metrics.

In 2014, McMurdie and Holmes [@McMurdie2014] announced that rarefaction of microbial community data is "statistically inadmissible" because it omits valid data. In their simulations, they observed that the use of rarefaction reduced the statistical power to compare communities based on beta diversity metrics. The communities used in this study were simulated by mixing operational taxonomic unit (OTU) frequencies from marine and human fecal communities, which were then sampled to generate 40 replicate communities each sampled to different depths. The simulation with the greatest depth of coverage obtained an average of 10,000 sequences per sample. Analysis of the R code that they used indicates that the simulations were actually only based on a single subsample and not truly rarefied data. Furthermore, the authors removed any OTU from the simulated communities that was not observed in three or more samples. This is a frequently used approach for analyzing microbiome data [@Bokulich2013]. Unfortunately, it removes valid data, skews the distribution of the community, and removes true patchiness from the representation of the communities [@Schloss2020]. Finally, the analysis did not include analysis of alpha diersity metrics and chose to assess beta diversity metrics based on clustering quality of samples rather than using the more frequently used approach of employing non-parametric multivariate analysis of variance [@McArdle2001]. Regardless of these caveats, McMurdie and Holmes recommended the use of variance stabilizing transformations such as those applied for gene expression analysis such as that found in the DeSeq R package. 

Beyond the critiques from McMurdie and Holmes others have developed alternative approaches to controlling for uneven sampling effort in amplicon sequencing studies. For alpha diversity metrics, Willis used toy datasets to demonstrate that one could estimate the richness for each sample in a dataset and use those values for statistical comparisons [@Willis2019]. Non-parametric estimators of richness and diversity [@Chao2003; @Chao2016] and parametric estimators of richness [Willis2015] that have been used in microbial ecology studies. For beta diversity metrics at least four approaches have been pursued. First, one could use relative abundance values where the observed number of sequences in an OTU is divided by the total number of sequences in the sample [@Lin2020]. Second, normalization strategies have been developed where the relative abundance is mutiplied by the size of the minimum desired sampling effort and fractional values are reapportioned to the OTUs to obtain integer values [@Beule2020; @Paulson2013]. Third, a variety of center log-ratio methods have been developed where the compositional nature of the OTU counts is removed and used to calculated Euclidean distances (aka Aitchison distances) [@Quinn2019; @Paulson2013; @Lin2020; @Martino2019; @Costea2014. This strategy is purported to control for uneven sampling effort [@Martino2019; @Quinn2018]; however, some have noticed that this feature breaks down under certain conditions [@teBeest2021]. Finally, as mentioned above, variance stabilization transformations have been recommended to generate values that can be used to calcualte Euclidean distances [@McMurdie2014].

The ongoing controversy over the use of rarefaction and the development of numerous strategies to control for uneven sampling effort caused me to question how these methods compared to each other. My analysis included 16S rRNA gene sequence data from from 12 studies that characterized the variation in bacterial communities from diverse environments (Table 1 and Figure S1). The sequences were assigned to OTUs using a standard pipeline and their frequencies and the number of sequences found in each sample were used to generate simulated communities and treatment effects. For each dataset and simulation, 100 replicate datasets were generated and used as input to each of the strategies for controlling for uneven sampling effort. My overall conclusion was that rarefaction outperformed the alternative strategies.

## Results

***Without rarefaction, metrics of alpha diversity are sensitive to sampling effort.*** To test the sesitivity of various approaches of measuring alpha diversity to sampling effort, I generated null models for each dataset. Under a null model, each community from the same dataset would be expected to have the same alpha diversity regardless of the sampling effort. I measured the richness of the communities in each dataset without any correction, using scaled ranked subsampling (SRS) normalized OTU counts, with estimates based on non-parametric and parametric approaches, and using rarefaction (e.g. Figure S2). For each dataset, all of the approaches, except for rarefaction, showed a strong correlation between richness and the number of sequences in the sample (Figure 1A). Next, I assessed diversity using the Shannon diversity index and the inverse Simpson diversity index without any correction, using normalized OTU counts, and rarefaction; I also used a non-parametric estimator of Shannon diversity. The correlation between sampling depth and the diversity metric was not as strong as it was for richness and the inverse Simpson diversity values were less sensitive than the Shannon diversity values; however, the correlation to the rarefied diversity metrics were the lowest for all of the metrics and studies (Figure 1A). The rarefied alpha-diversity metrics consistently demonstrated a lack of sensitivity to sampling depth.

***Without rarefaction, metrics of beta diversity are sensitive to sampling effort.*** To test the sesitivity of various approaches of measuring beta diversity to sampling effort, I used the same null models used for studying the sensitivity of alpha diversity. Under a null model, the ecological distance between any pair of samples would be the same regardless of the difference in the number of sequences observed in each sample (e.g., Figure S3). First, I calculated the Jaccard distance coefficient between all pairs of communities within a dataset. The Jaccard distance coefficient is the fraction of OTUs that are unique to either community and does not account for the abundance of the OTUs. Jaccard distances were calcualted using the uncorrected OTU counts, with rarefaction, relative abundances, and following normalization using cumulative sum scaling (CSS) and SRS. Only the rarefied distances showed a lack of sensitivity to sampling effort (Figure 1B). Second, I analyzed the sensitivity of the Bray-Curtis distance coefficient, which is a popular metric that incorporates the abundance of each OTU. Similar to what I observed with the Jaccard coefficient, only the rarefied data showed a lack of sensitivity to sampling effort (Figure 1B). Third, I calcualted the Euclidean distance on raw OTU counts where the central log-ratio (CLR) was calculated (i.e., Aitchison distances) by ignoring OTUs in samples with zero counts (Robust CLR), adding a pseudocount of one to all OTU counts prior to calculating the CLR (One CLR), adding a pseudocount of one divided by the total number of sequences obtained for the community (Nudge CLR), and imputing the value of zero counts (Zero CLR). The Aitchison distances were all strongly sensitive to sampling effort (Figure 1B). Finally, I used the variance-stabilization transformation (VST) from DeSeq2 prior to calculating Euclidean distances. Again, there was a strong sensitivity to sampling effort (Figure 1B). Although Euclidean distances are not typically used on raw or rarefied count data in ecology, rarefied Euclidean distances were not sensitive to sampling effort. Across each of the beta diversity metrics and approaches used to account for uneven sampling effort and sparsity, rarefaction was the least sensitive approach to differences in sampling effort.

***Rarefaction limits the detection of false positives when sampling effort and treatment group are confounded.*** Next, I investigated the impact of the various strategies and metrics on falsely detecting a significant difference using the the same communities generated from the null model in the analysis of alpha and beta diversity metrics. To test for differences in alpha and beta diversity I used the non-parametric Wilcoxon test and non-parametric permutation-based multivariate analysis of variance (PERMANOVA). First, I employed an unbiased null treatment model to measure the false detection rate, which should not have meaningfully differed from 5%. Indeed, for each dataset and alpha and beta diversity metric and strategy for accounting for uneven sampling, approximately 5% of the tests yielded a significant result (Figure 2). Second, I employed a biased null treatment model where the treatment group was completely confounded with the number of sequences in each sample. Under this model, only the rarefied data consistently resulted in a 5% false positive rate for alpha and beta diversity metrics (Figure 2). These results aligned with the observed sensitivity of alpha and beta diversity metrics to sampling effort and underscore the value of rarefaction.

***Rarefaction preserves the statistical power to detect differences between treatment groups.*** To assess the impact of different approaches to control for uneven sampling effort I performed two additional simulations. In the first simulation, I implemented a skewed abundance distribution model to create two treatment groups for each datasets that were each populated with half of the samples each with the same number of sequences as the samples had in the observed data. The power to detect differences in richness between the two simulated treatment groups by all approaches was low (Figure 4A). This was likely because the approach for generating the perturbed community did not necessarily change the number of OTUs in each treatment group. Regardless, the simulations testing differencse in richness using the Rice and Stream datasets had the greatest power when the richness data were rarefied. To explore this further, a richness-adjusted community model was created by removing 3% of the OTUs from a null model model. As suggested by the first simulation, the rarefied richness data had a higher statistical power than the other approaches when measuring richness (Figure 5). The simulations testing the power to detect differences in Shannon diversity also showed that rarefied data performed other methods (Figure 4A). When testing for differences in the Inverse Simpson diversity index the the difference between rarefaction and the other methods was negligible (Figure 4A). For tests of beta diversity I found that rarefaction was the most reliable approach to maintain statistical power to detect differences between two communities. Among the tests using the Jaccard and Bray-Curtis metrics, raw count data and CSS normalized data had little power relative to rarefied, relative abundance, and SRS normalized data. The differences in power between rarefied, relative abundance, and SRS normalized data was small, but if there were differences, the power obtained using rarefied data was greater than the other methods. Among the tests using Euclidean distances, using raw counts and CLR and DeSeq2 transformed data had little power relative to the distances calcualted using rarefied and relative abundance data. This power-based analysis of the simulated communities using different methods of handling uneven sample sizes demonstrated the value of rarefaction for preserving the statistical power to detect differences between treatment groups for measures of alpha and beta diversity.



***Increased rarefaction depth reduces intra-sample variation in alpha and beta diversity.*** Once concern with rarefying communities is the perceived loss of sequencing information when more a large fraction of data appears to be removed when the community with the greatest sequencing depth is rarefied to the size of the community with the least (e.g., 99% with the Bioethanol dataset). To assess the sensitivity of alpha and beta diversity metrics to rarefaction depth, I again used the dataset generated using the null models, but rarefied each community to varying sampling depths (Figure 6). The richness values increased with sampling effort as rare OTUs would continue be detected. In contrast, the Shannon diversity and Bray-Curtis values plateaued with increased sampling effort. This result was expected since increased sampling would lead to increased precision in the measured abundance of OTUs. Next, I measured the coefficient of variation (i.e., the mean divided by the standard deviation) between samples for richness, Shannon diversity, and Bray-Curtis distances. Although the richness values appeared to increased unbounded with smapling effort, the coefficient of variation for each dataset was relatively stable. In general, the coefficient of variation increased slightly with sampling depth only to decline once smaller samples were removed from the analysis at higher sampling depths. Interestingly, the coefficient of variation between Shannon diveristy values decreased towards zero with increased sampling effort and the coefficient of variation between Bray-Curtis distances tended to increased. Regardless, the coefficients of variation were relatively small.

***Most studies have a high level of sequencing coverage.*** To explore the concern over loss of sequencing depth further, I calculated the Good's coverage for the observed data. The median coverage for each dataset ranged between 89.4 and 99.8% for the Seagrass and Human datasets, respectively (Figure 7). When I rarefied each dataset to the size of the smallest community in the dataset, with the exception of the Seagrass, Rice, and Stream datasets, the median coverage for the rarefied communities was still greater than 90%. These results suggest that most studies had a level of sequencing coverage that aligned with the diversity of the communities. Next, I used the null model for each dataset to ask how much sequencing effort was required to obtain higher levels of coverage. To obtain 95 and 99% coverage, an average of 2.70 and 101.20-fold more sequence data was estimated to be required than was required to obtain 90% coverage, respectively (Figure 7). As suggested by the simulated coverages curve in Figure 7, these estimates are conservative. Regardless, the sequencing effort required to acheive higher sequencing depth would likely limit the number of samples that could be sequenced when controlling for costs. Although it may be disconcerting to rarefy to a sequencing depth that is considerably lower than that obtained for the best sequenced community in a dataset, sequencing coverage for many environments is probably adequate even at the lower sequencing depth. Of course, the results above have demonstrated that rarefaction is necessary to avoid problems with making inferences.

## Discussion

Over the past decade, the question of whether to rarefy microbial community sequence data has become controversial. The analyses I presented here strongly indicate that rarefaction is necessary to control for uneven sampling effort when comparing communities based on alpha and beta diversity indices. Compared to all other methods, rarefaction removed the correlation between sequencing depth and alpha or beta diversity metrics when the sequencing depth varied by as much as 97-fold across samples. I showed that this correlation could lead falsely detecting differences between treatment groups if sampling depth and sequencing effort are confounded [e.g., @Morris2013]. The correlation with sampling effort leads to an artificial increase in the variation between samples and a reduced power to detect true differences in alpha and beta diversity. For these reasons, rarefaction is a valuable tool to control for uneven sampling effort until improved statistical procedures are developed or it becomes possible to more evenly distribute sequencing effort across samples.

The primary alternative to rarefaction that many advise is to estimate alpha diversity has been to use non-parametric or parametric methods using raw data and to then compare the estimates [@Willis2019]. My data demonstrates that such approaches are limited for several reasons. First, it has been widely observed that non-parametric richness esimates such as ACE and Chao1 are sensitive to sampling effort. Therefore, these strategies do not, in practice, remove the effects of sampling effort. Second, parametric approaches, such as those implemented in the breakaway R package, generate confidence intervals that are likely to include the true richness and theoretically shrink with increased sampling effort. For most samples, the confidence intervals around the estimates are too wide to be meaningful, again leading to an inability to remvoe the effects of sampling effort. Third, it has become an increasingly common practice for researchers to remove sequences that are rare in a sample (e.g., those that appear once). Although that approach was not taken in this study, removing rare sequences would skew the distribtion of sequences and OTUs leading to a distortion of the measurement of alpha diversity [@Schloss2020]. The effects of removing rare sequences would vary across samples depending on the number of sequences in each sample. One interesting result of this analysis was the demonstration that as metrics that depend less on rare taxa are used, the effect of uneven sampling effort was reduced. For example, richness counts a sequence appearing once as much as sequence appearing 1,000 times, while Shannon diversityindex places more emphasis on the more abundant sequence, and the inverse Simpson index even more. Although normalizing communities to a common number of sequences is also suggested (e.g., SRS normalization) to control for uneven sampling effort, the current analysis indicates that its performance does not meet that of rarefaction. For analysis of alpha diversity metrics, rarefaction is the most effective and consistent approach to control fo uneven sampling effort.

Use of relative abundances, normalized counts, variance stabilizing transformations, and centered log ratios have each been recommended as more robust alternatives to rarefaction. Again, the only approach that consistently removed the effects of uneven sampling effort was rarefaction. Most of the recommendations borrow techniques from methods for identifying differential gene expression. Unfortunately, there is an important but underappreciated difference between gene expression and community data. This is the interpretation of unobserved data. For gene expression analysis in a single organism the lack of any sequencing data for a gene would indicate that its expression was below the limit of detection. Sequencing the same organism under multiple conditions would not lead to a seemingly unbouned number of genes in the organims. Rather, the number of genes has a definite limit that is knowable from the genome sequence. With microbiome data, an unobserved sequence could mean that the organism was present, but below the limit of detection or that the organism was missing. Because we have yet to exhaustively sample any community in the same way we have sequenced a single genome it is unreliable to impute the presence of all organisms. Yet, this is exactly what the variance stabilization transformation and most CLR techniques do. This analysis has demonstrated a clear correlation between distances calculated by these methods and sequencing effort. This result is at odds with the claims by others that the distances are scale invariant [@Martino2019; @Quinn2018]. Again, rarefaction is the most effective and consistent approach to control for uneven sampling effort when calculating beta diversity metrics.

Two common critique of rarefaction is that the approach ommits valid data and that the selection of the value to rarefy the sequence data to is arbitrary [@McMurdie2014]. I disagree with both critiques. All of the data are used to calculate the mean value of the rarefied metrics; none of it is excluded. When the dataset is subsampled, every sequence has a random chance of being included in the calculated metric. When that subsampling is repeated a large number of times (e.g., 100 or 1000) the risks of ignoring or oversampling rare taxa are mititgated. Furthermore, it is curious that the study making the original critique removed sequences that did not appear in at least three samples or had an abundance of one in any sample. A parallel analysis to this study has demonstrated that many of these sequences are likely valid and that removal of rare sequences can bias alpha and beta diversity metrics and reduce statistical power [@Schloss2020]. As for the second critism, I would resist the claim that the selection of the rarefaction threshold is arbitrary. In practice, there is a tradeoff between sampling breadth and depth. Greater breadth will increase the statistical power to compare treatment groups and greater depth will increase the resolution to describe the communities. My process for picking a break involves looking for a natural break in the distribution of the number of sequences. For example, the Lake dataset used in this study had a clear break around 10,000 sequences. I would also consider what samples are below any break that I select. If there were critical samples below the break I would either reduce the break point or I would generate more sequences from those samples. Furthermore, regardless of the break point one selects, they could repeat an analysis at multiple thresholds to assess the tradeoff between statistical power and resolution. As shown in my analysis of Good's coverage values, most studies obtain an ample level of coverage and would need to increase their sampling depth by 10 fold to increase the coverage by several percentage points. I contend that for most studies increased sequencing effort should be applied to improve sampling breadth and statistical power over depth and resolution.

The up to 100-fold difference in sample sizes is an unfortunate byproduct of how seqencing libraries are constructed. Researchers perform separate PCRs for each sample with unique index (aka, barcode) sequencs that all them to later identify which sample each sequence belongs. When the PCRs are pooled efforts are taken to pool the fragments in equimolar ratios. Researchers use one of two approaches. First, they often will quantify the concentration of DNA from each PCR and then pool DNA in the desired amounts. Alternatively, they may use normalization plates where each well can hybridize a uniform amount of DNA that is then eluted and pooled. Clearly, both approaches have limitations that reduce the ability to truly achieve equimolar mixture. In addition, for some samples it is common to co-amplify non-specific DNAs which may add to the challenges of achieving equimolar mixtures of the desired amplicons [@Morris2013]. Regardless, it is clear that better strategies are needed to reduce the variation in the number of sequences generated for each sample.

All simulations have weaknesses and should be interpreted with caution. However, the simulated communites generated and analyzed in this study had the advantage of being designed with known properties including the alpha and beta diversity and the their differences between treatment groups. For now, is is perfectly admissible and proper to rarefy microbial community data or risk reaching unwarrented conclusions.


## Materials and Methods

**Choice of datasets.** The specific studies used in this study were selected because their data was publicly accessible through the Sequence Read Archive, the original investigators sequenced the V4 region of the 16S rRNA gene using paired 250 nt reads, and my previous familiarity with the data. The use of paired 250 nt reads to sequence the V4 region resulted in a near complete two-fold overap of the V4 region resulting in high quality contigs with a low sequencing error rate [@Kozich2013]. These data were processed through a standardized sequence curation pipeline to generate operational taxonomic units (OTUs) using the mothur software package [@Kozich2013; @Schloss2009]. OTUs were identified using the OptiClust algorithm to cluster sequences together that were not more than 3% different from each other [@Westcott2017].

**Null community model.** Null community models were generated such that within a dataset the number of sequences per sample and the number of sequences per OTU across all samples within the dataset were the same as was observed in original. This model effectively generated statistical samples of a single community so that there should not have been a statistical difference between the samples. This model implemented by randomly assigning each sequence in the dataset to an OTU and sample while keeping constant the number of sequences per sample and the total number of sequences in each OTU. This is a similar approach to that of the IS algorithm described by Ulrich and Gotelli [@Ulrich2010]. Because the construction of the null models was a stochastic process, 100 replicates were geneated for each dataset.

**Null treatment models.** I created an unbiased and biased treatment model. In the unbiased model, samples were randomly assigned to one of two treatment groups. In the biased treatment model, samples that had more than the median number of seqeunces for a dataset were assigned to one treatment group and the rest were assigned to a second treatment group. Comparison of any diversity metric across the two treatment groups should have only yielded a significant result in 5% of the simulations when testing under a Type I error (i.e., $\upalpha$) of 0.05. 

**Skewed abundance community model.** In the skewed abundance community model, communities were randomly assigned to one of two simulated treatment groups. Communities in the first treatment group were generated by calculating the relative abundance of each OTU across all samples and using those values as the probability of sampling each OTU. This probability distribution was sampled until each sample had the same number of sequences that it did in the observed data. Samples in the second treatment group were generated by adjusting the relative abundances of the OTUs in the firs treatment group by increasing the relative abundance of 10% of the OTUs by 5%. These values were determined after empirically searching for conditions that resulted in a large fraction of the randomizations yieleding a significant result across most of the studies. Sequences were sampled from the skewed community community until each sample had the same number of sequences that it did in the observed data. Under the skewed abundance community model each sample represented a statistical sampling of two communities such that there should not have been a statistically significant difference within a treatment group, but there was between the treatment groups. Because the construction of the skewed abundance community model was a stochastic process, 100 replicates were geneated for each dataset.

**Richness-adjusted community model.** In the richness-adjusted community model, communities were randomly assigned to one of two simulated treatment groups. Communities in the first treatment group were generated by calculating the relative abundance of each OTU across all samples and using those values as the probability of sampling each OTU. This probability distribution was sampled until each sample had the same number of sequences that it did in the observed data. Samples in the second treatment group were generated by removing 3% of the OTUs from the dataset and recalculating the relative abundance of the remaining OTUs. Sequences were sampled from the richness-adjusted community distrubtion until each sample had the same number of sequences that it did in the observed data. Under the richness-adjusted community model each sample represented a statistical sampling of two communities such that there should not have been a statistically significant difference within a treatment group, but there was between the treatment groups. Because the construction of the richness-adjusted community model was a stochastic process, 100 replicates were geneated for each dataset.

**Test of statistical significance.** Statistical comparisons of alpha diversity metrics across the simulated treatment groups were performed using the non-parametric two-sample Wilcoxon test as implemented in R. This test was selected because the alpha-diversity metrics tended to not be normally distributed and each dataset required a different transformation to normalize the data. Comparisons of beta diversity metrics were performed using the adonis2 function from the vegan (v.2.6.2) R package [@Oksanen2022]. The adonis2 function implements a non-parametric multivariate analysis of variance using distances matrices as described by McArdle and Anderson [@McArdle2001]. Throughout this study I used 0.05 as the threshold for asssessing statistical significance.

**Power analysis.** The parameters used to design the skewed abundance and richness-adjusted community models were set to impose a known effect size when using rarefied community data. The statistical power to detect these differences was determined by calculating the p-value for each of 100 replicate simulated set of samples from each dataset using the various alpha and beta diversity metrics. The percentage of tests that yielded a significant P-value was considered the statistical power (i.e., 1 minus the Type II error) to detect the difference.

**Alpha diversity calculations.** Various strategies for handling uneven sampling effort were evaluated to identify the best approach for calculating community richness and Shannon and inverse Simpson diversity indeices. Raw OTU counts were used as input to calculate sample richness and Shannon and inverse Simpson diversity using mothur [@Schloss2009; @Magurran2004]. Shannon diversity was calculated as

$$H_{shannon} = - \sum_{i=1}^{S_{obs}} \frac{n_i}{N} ln \frac{n_i}{N}$$

The Simpson diversity was calculated as

$$D_{simpson} = \frac {\sum_{i=1}^{S_{obs}} {n_i \left ( n_i - 1 \right )}}{N \left( N-1 \right )}$$

The inverse Simpson diversity was calculated as $1/D_{simpson}$. In both formulae, $n_i$ is the number of sequences in OTU $i$ and $N$ is the number of sequences in the sample. Rarefaction of richness, Shannon diversity and Inverse Simpson diversity values were carried out in mothur. Briefly, mothur calculates each value on a random draw of the same number of sequences from each sample and obtains a mean value based on 1,000 random draws. Scaled ranked subsampling (SRS) was used to normalize OTU counts to the size of the smallest sample in each dataset using the SRS R package (v.0.2.3)[@Beule2020]. Normalized OTU counts were used to calculate sample richness and Shannon and inverse Simpson diversity values using mothur. Data normalized by cumulative sum scaling (CSS) were not reported for alpha-diversity values since the relative abundances of the features do not change with the normalization procedure [@Paulson2013]. The non-parametric bias-corrected Chao1 and ACE richness estimators [@Chao2016] and a non-parametric estimator of the Shannon diversity [@Chao2003] were calculated using raw OTU counts with mothur. Parametric estimates of sample richness were calculated using the breakaway (BA) R package (v.4.7.9)[@Willis2015]. The current analysis reports both the results from running default model selection procedure and the Poisson model. The default model selection returned either the Kemp, Negative Binomial, or Poisson models. Relative abundance data were not used to calculate alpha diversity metrics since the richness and evenness does not change from the raw data when dividing each sample by the total number of sequences in the sample. 

**Beta diversity calculations.** Similar to the alpha diversity calculations, multiple approaches were used to control for uneven sampling effort and calculate beta diversity. Raw and OTU counts were used for input to calculate the Jaccard, Bray-Curtis, and Euclidean dissimilarity indices using the vegdist function from the vegan R package (v.2.6.2)[@Oksanen2022]. The Jaccard index was calculated as

$$D_{Jaccard}=1-\frac{S_{AB}}{S_A+S_B-S_{AB}}$$

where $S_A$ and $S_B$ were the number of OTUs in samples $A$ and $B$ and $S_{AB}$ was the number of OTUs shared between the two samples. The Bray-Curtis index was calculated as 

$$D_{Bray-Curtis}=1-\frac{\sum_{i=1}^{S_T} \left| n_{A,i} - n_{B,i}\right| }{ N_A + N_B}$$

where $n_{A,i}$ and $n_{B,i}$ are the number of sequences observed in OTU $i$ from samples $A$ and $B$, respectively. $N_A$ and $N_B$ are the total number of sequences in samples $A$ and $B$, respectively. $S_T$ is the total number of OTUs observed between the two samples. The Euclidean distance was calculated as

$$D_{Euclidean}=1-\sqrt{\sum_{i=1}^{S_T}\left(n_{A,i} - n_{B,i}\right)^2}$$

These metrics were calculated using the relative abundance of each OTU using the vegdist function from vegan. The relative abundance was calculated as the number of sequences in the OTU (e.g., $n_{A,i}$) divided by the total number of sequences in the sample (e.g., $N_A$).

Rarefied beta-diversity values were calculated using the avgdist function in vegan. Briefly, vegan's avgdist function calculates each pairwise dissimilarity index after obtaining a random draw of the same number of sequences from each sample. After obtaining 100 random draws it returns the mean value.

Three approaches were taken to normalize the number of sequences across samples within a dataset. Scaled ranked subsampling (SRS) and cumulative sum scaling (CSS) were used to normalize raw OTU counts using the SRS (v.0.2.3) and metagenomeSeq (v.1.36.0) R packages, respectively [@Beule2020; @Paulson2013]. The normalized counts were then used to calculate Jaccard and Bray-Curtis dissimilarity indices. Finally, the variance-stabilization transformation (VST) as implemented in the DESeq2 (v.1.34.0) R package was used to normalize the data as described by McMurdie and Holmes [@Love2014; @McMurdie2014]. Because the VST approach generates negative values, which are incompatible with calculating Jaccard and Bray-Curtis dissimilarity values, Euclidean distances were calcualted instead.

Raw OTU counts were used to calculate centered log ratio values for each OTU, which were then used to calcualte Euclidean distances; such distances are commonly refered to as Aitchison distances. Centered log-ratio (CLR) abundances are calculated as:

$$
clr\left(n_j\right) = \left[ \ln\frac{x_{ij}}{g(x_j)}, ..., \ln\frac{x_{S_Tj}}{g(x_j)}\right]
$$

where $x_{ij} was the number of sequences observed for OTU $i$ in sample $j$ and $g()$ was the geometric mean $x_j$ was the counts of the $S_T$ OTUs in sample $j$. Because the geometric mean is zero if any OTU is absent from a sample, the CLR is undefined when there are unobserved OTUs in a sample. To overcome this problem, I attempted a four approaches. The first, Zero CLR, imputed the value of the zeroes based on the observed data using the zCompositions (v.1.4.0.1) R package [@Palarea2015]. The second, One CLR, added a pseudocount of 1 to the abundance of all OTUs [@Paulson2013; @Lin2020]. The third, Nudge CLR, added a pseudocount of 1 divided by the total number of sequences in a sample to each OTU in the sample [@Costea2014; @Lin2020]. The final approach, Robust CLR, removed unobserved OTUs prior to calcualting the CLR [@Martino2019].

**Analysis of sequencing coverage.** To assess the level of sequencing coverage I calcualted Good's coverage ($C_{Good}$) using mothur:

$$C_{Good} = 100\% \times \left(1-\frac{n_1}{N_T} \right)$$

where $n_1$ is the number of OTUs with only one sequence in the sample and $N_T$ is the total number of sequences in the sample. Good's coverage was calculated using (i) the observed OTU counts for each sample and dataset, (ii) following rarefaction (1,000 iterations) of the observed OTU counts to the size of the smallest sample in each dataset, and (iii) after rarefying (1,000 iterations) the null community distribution.

**Reproducible data analysis.** A complete reproducible workflow written in Snakemake (v.7.15.2) and conda (v.4.12.0) computational environment can be obtained from the GitHub hosted git repository for this project (https://github.com/SchlossLab/Schloss_Rarefaction_XXXXX_2022). This paper was written in R markdown (v.2.16) with the kableExtra (v.1.3.4) package. The mothur (v.1.47.0) and R (4.1.3) software packages were used for all analyses with extensive use of functions in the tidyverse metapackage (v.1.3.1). A preliminary version of this analysis was presented as the Rarefaction video series on the Riffomonas Project YouTube channel (https://www.youtube.com/playlist?list=PLmNrK_nkqBpJuhS93PYC-Xr5oqur7IIWf).

\vspace{10mm}

**Acknowledgements.** 

I am grateful to the researchers who generated the datasets used in this study. I also thank the individuals who asked questions and commented on the preliminary version of this project, which was released as a YouTube playlist on the Riffomonas channel. This work was supported in part by funding from the National Institutes of Health (U01AI124255, P30DK034933, R01CA215574).

\newpage

## References

\setlength{\parindent}{-0.25in}
\setlength{\leftskip}{0.25in}
\noindent

<div id="refs"></div>
\bibliography{ref}
\setlength{\parindent}{0in}
\setlength{\leftskip}{0in}

\newpage
## Tables

**Table 1. Summary of studies used in the analysis.** For all studies, the number of sequences used from each dataset was rarefied to the smallest sample size. A graphical represenation of the distribution of sample sizes for each dataset and the samples that were removed from each dataset are provided in Figure S1.

\footnotesize


|\textbf{Dataset\nobreakspace{}(Ref)} | \textbf{Samples}| \makecell[c]{\textbf{Total}\\\textbf{sequences}}| \makecell[c]{\textbf{Median}\\\textbf{sample size}}| \makecell[c]{\textbf{Mean}\\\textbf{sample size}}| \makecell[c]{\textbf{Range of}\\\textbf{sample sizes}}| \makecell[c]{\textbf{SRA study}\\\textbf{accession}}|
|:------------------------------------|----------------:|------------------------------------------------:|---------------------------------------------------:|-------------------------------------------------:|------------------------------------------------------:|----------------------------------------------------:|
|Bioethanol\ [@Li2015]                |               95|                                        3,970,972|                                              16,014|                                            41,799|                                 3,690\Hyphdash*356,027|                                            SRP055545|
|Human\ [@Baxter2016]                 |              490|                                       20,828,275|                                              32,452|                                            42,506|                                10,439\Hyphdash*422,904|                                            SRP062005|
|Lake\ [@Beall2015]                   |               52|                                        3,145,486|                                              69,205|                                            60,490|                                15,135\Hyphdash*110,993|                                            SRP050963|
|Marine\ [@Henson2016]                |                7|                                        1,484,068|                                             213,091|                                           212,009|                               132,895\Hyphdash*256,758|                                            SRP068101|
|Mice\ [@Kozich2013]                  |              348|                                        2,785,641|                                               6,426|                                             8,004|                                  1,804\Hyphdash*30,311|                                            SRP192323|
|Peromyscus\ [@Baxter2014]            |              111|                                        1,545,288|                                              12,393|                                            13,921|                                  4,454\Hyphdash*33,502|                                            SRP044050|
|Rainforest\ [@LevyBooth2018]         |               69|                                          936,666|                                              11,464|                                            13,574|                                  4,880\Hyphdash*37,403|                                            ERP023747|
|Rice\ [@Edwards2015]                 |              490|                                       22,623,937|                                              43,399|                                            46,171|                                 2,777\Hyphdash*192,200|                                            SRP044745|
|Seagrass\ [@Ettinger2017]            |              286|                                        4,135,440|                                              13,538|                                            14,459|                                  1,830\Hyphdash*45,076|                                            SRP092441|
|Sediment\ [@Graw2018]                |               58|                                        1,151,389|                                              17,606|                                            19,851|                                  7,686\Hyphdash*67,763|                                            SRP097192|
|Soil\ [@Johnston2016]                |               18|                                          932,563|                                              50,487|                                            51,809|                                 46,622\Hyphdash*58,935|                                            ERP012016|
|Stream\ [@Hassell2018]               |              201|                                       21,017,610|                                              90,621|                                           104,565|                                 8,931\Hyphdash*394,419|                                            SRP075852|
\normalsize

\newpage

## Figures

\includegraphics[height=17cm]{figure_1.png}

**Figure 1. Rarefaction eliminates the correlation between sequencing depth and alpha diversity (A) and between differences in sampling depth and beta (B) diversity metrics when using null community models.** Examples of the relationship between different metrics and methods for controlling for uneven sequencing effort are provided in Figures S2 and S3 for alpha and beta diversity metrics, respectively. Each point represents the mean of 100 random null community models; the standard deviation was smaller than the size of the plotting symbol.

\newpage

\includegraphics[height=17cm]{figure_2.png}

**Figure 2. The risk of falsely detecting a difference between treatment groups drawn from a null model does not meaningfully vary from 5%, regardless of approach for controlling for uneven sequencing depth**. Samples were randomly assigned to different treatment groups. To calculate the false detection rate, datasets were regenerated 100 times and differences in alpha diversity were tested using a Wilcoxon test (A) and differences in beta diversity were tested using PERMANOVA (B) at a 5% threshold. The false positive rate was the number of times a dataset yeilded a significant result.

\newpage

\includegraphics[height=17cm]{figure_3.png}

**Figure 3.The risk of falsely detecting a difference between treatment groups drawn from a null model does not meaningfully vary from 5% when data are rarefied when sequencing depth is confounded with treatement group**. Samples were assigned to different treatment groups based on whether they were above the median number of sequences for each dataset. To calculate the false detection rate, datasets were regenerated 100 times and differences in alpha diversity were tested using a Wilcoxon test (A) and differences in beta diversity were tested using PERMANOVA (B) at a 5% threshold. The false positive rate was the number of times a dataset yeilded a significant result.


\newpage

\includegraphics[height=17cm]{figure_4.png}

**Figure 4. The ability to detect true differences in treatment groups for alpha (A) and beta (B) diversity metrics is greatest when communities differing in the relative abundance of their OTUs are rarefied.** For each dataset samples were randomly assigned to one of two community distributions where the abundance of OTUs differed. To calculate the power for each datasets, datasets were regenerated 100 times and differences in alpha diversity were tested using a Wilcoxon test (A) and differences in beta diversity were tested using PERMANOVA (B) at a 5% threshold. The power was the number of times a dataset yielded a significant result.


\newpage

\includegraphics{figure_5.png}

**Figure 5. The ability to detect true differences in treatment groups for alpha diversity metrics is greatest when communities differing in richness are rarefied.** For each dataset samples were randomly assigned to one of two community distributions where one distribution contained a subset of OTUs found in the other. To calculate the power for each dataset, datasets were regenerated 100 times and differences in alpha diversity were tested using a Wilcoxon test (A) and differences in beta diversity were tested using PERMANOVA (B) at a 5% threshold. The power was the number of times a dataset yielded a significant result. 


\newpage

\includegraphics{figure_6.png}

**Figure 6. The mean and coefficient of variation for rarefied richness, shannon diversity, and Bray-Curtis dissimilarity vary with sequencing depth.** For each dataset, a null community distribution was created and samples were created to have the same sequencing depth as they did originally. The placement of the plotting symbol indicates the size of the smallest sample. Results are only shown for sequencing depths where a dataset had 5 or more samples.

\newpage

\includegraphics{figure_7.png}

**Figure 7. Most datasets are sequenced to a level that provides a high level of coverage.** Each plotting symbol represents the observed Good's coverage for a different sample in each dataset. The smoothed line indicates the simulated coverage for varying levels of sampling effort when a null community is generated from the observed data. The box and whisker plot indicates the range of coverage values when the observed commmunity data were rarefied to the size of the least sequenced sample.

\newpage

\includegraphics{figure_s1.png}

**Figure S1. The number of sequences observed in each sample for each dataset included in this analysis generally varied by 10 to 100-fold.** The threshold for specifying the number of sequences per sample varied by dataset and was determined based on identifying natural breaks in the data.

\newpage

\includegraphics{figure_s2.png}

**Figure S2. Examples of the richness in each of the 490 samples that were generated for one randomization of the null model using the human dataset.** The x-axis indicates the number of sequences in each of the samples prior to each method's appraoch of controlling for uneven sampling effort. The Spearman correlation coefficient ($\uprho$) and test of whether the coefficient was significantly different from zero are indicated for each panel.

\newpage

\includegraphics{figure_s3.png}

**Figure S3. Examples of differences in beta diversity in each of the 490 samples that were generated for one randomization of the null model using the human dataset.** The x-axis indicates the difference in the number of sequences in each of the samples that went into calcualting the pairwise distance prior to each method's appraoch of controlling for uneven sampling effort. The Spearman correlation coefficient ($\uprho$) and test of whether the coefficient was significantly different from zero are indicated for each panel.
