# Social Support and Network Formation in a Small-Scale Horticulturalist Population (Simpson, Accepted, Scientific Data)


## Abstract
Evolutionary studies of cooperation in traditional human societies suggest that helping family and responding in kind when helped are the primary mechanisms for informally distributing resources vital to day-to-day survival (e.g., food, knowledge, money, childcare). However, these studies generally rely on forms of regression analysis that disregard complex interdependences between aid, resulting in the implicit assumption that kinship and reciprocity drive the emergence of entire networks of supportive social bonds. Here I evaluate this assumption with actor-based simulations of network formation (i.e., Stochastic Actor-Oriented Models). Specifically, I test standard predictions of cooperation derived from the evolutionary theories of kin selection and reciprocal altruism alongside well-established sociological predictions around the self-organisation of asymmetric relationships. Simulations are calibrated to exceptional public data on genetic relatedness and the provision of tangible aid amongst all 108 adult residents of a village of indigenous horticulturalists in Nicaragua (11,556 ordered dyads). Results indicate that relatedness and reciprocity are markedly less important to whom one helps compared to the supra-dyadic arrangement of the tangible aid network itself.

![](https://github.com/cohensimpson/smallnet_ScientificData/blob/main/Figure%201%20SciData%20Social%20Support%20and%20Network%20Formation%20(2022-03-02).svg "Kinship and tangle aid in Arang Dak. Each arc (i.e., directed relationship) indicates the provision of tangible aid by some villager i to some other villager j (108 villagers; 1,485 arcs, 1,422 of which are verified by the source and target of aid). Arcs are coloured to reflect whether the source and target of aid are “close kin” (i.e., consanguineal relatedness or affinal relatedness ≥ 0.125). Red arcs emanate from spouses (i.e., affinal relatedness = 1.0) and dark-blue arcs emanate from primary kin — i.e., relatives with a consanguineal or affinal relatedness equal to 0.5 (e.g., full-siblings, parent and child, wife and brother-in-law, husband and mother-in-law). Lighter-blue arcs and green arcs emanate from near kin (i.e., 0.375 ≥ consanguineal/affinal relatedness ≥ 0.125). And yellow arcs emanate from kin who are not close, either genetically or through marriage (i.e., consanguineal/affinal relatedness ≤ 0.125), or emanate from residents of Arang Dak who are not related (i.e., consanguineal/affinal relatedness = 0).") _Kinship and tangle aid in Arang Dak. Each arc (i.e., directed relationship) indicates the provision of tangible aid by some villager i to some other villager j (108 villagers; 1,485 arcs, 1,422 of which are verified by the source and target of aid). Hover over network visualisation for full caption._



## R Code
Here, you will find a single R Script and two ".csv" data files. I have treated the R script a bit like a running "notebook". Accordingly, throughout the script, you will find code to carry out the analyses reported in my paper alongside commands used to produce useful print out (e.g., descriptive statistics, small tables, plots, etc.) and comments that (hopefully) give you insight into the thinking behind the decisions I take.

**_After_** you have placed the data files and the R script in the same R working directory, installed the necessary packages, and set the number of available computing cores for your machine (see circa Line 68 of the R script), you should be able to simply hit the "source" button in RStudio or run "source("RI_Nicaragua_Replication.R")" to redo my analyses. This will also carry out all of the goodness-of-fit tests and generate the text files containing the numbers used to produce Table 1, Table 2, Table 3, and Table 4 in the manuscript (N.B. Online-Only Table 1 is made by hand).



## Executables & Packages
In addition to the R scripts, I have included in the repository the installation files for the version of the two R packages integral to my analyses — i.e, the packages "RSiena" (https://github.com/snlab-nl/rsiena/wiki) and "groundhog" (https://groundhogr.com). Note that you may need to first install GCC (https://formulae.brew.sh/formula/gcc) — i.e., the GNU Compiler Collection — before attempting to install RSiena from source. Also, see the list of loaded packages at the very beginning of the R scripts for other short notes on dependencies that you may need to address.

For the unfamiliar, the groundhog package is a fabulous innovation that is designed to make package installation in the name of reproducible research very, very easy. Specifically, it uses the CRAN/MRAN database and date-based version control to load the packages necessary for a given set of analyses, as well as their dependencies. Please see the first 50-ish lines of the R scripts for details. Versions of RSiena are inconsistently pushed to CRAN/MRAN and thus RSiena will need to be installed from source using the file I have included in the repository.

Finally, when re-running my analyses, some numerical results may differ slightly from those reported in the paper due to stochastic perturbations. I have used the same random seed (20180709) to ensure exact reproducibility wherever possible. However, this is not always an option depending on the function.



## Summary of Files in Repository
 1) RI_Nicaragua_Replication.R (Script for Data Preparation, Transformation, Analyses, and Goodness-of-Fit)

 2) RSOS_corrected_data (14 May 2018).csv (Monadic Covariates + Dyadic Covariates + Network Data Collected by Koster (2018) for his Royal Society Open Science Paper — Excluding the Ethnicity Variable and Including Koster’s [1] Corrected Kinship Measures) 

 3) RSOS_corrected_data_attributes (17 May 2018).csv (Monadic Data Collected by Koster [1] Including the Ethnicity Variable) 

 4) RI_Arang_Dak_2021_Estimation History.txt" (Detailed Output from Each Iteration of the RSIENA Estimation Algorithm for the Eight Models in my Paper [N.B., Each Fully-Converged Model Requires Multiple Iterations of the RSIENA Algorithm]) 

 5) Koster_2018.pdf [1]
 
 6) groundhog-1.5.0.tar.gz ("groundhog" source code)

 7) Ripley et al. - 2022 - Manual for RSiena (v. 1.3.6).pdf" [2]

 8) rsiena-1.3.6.tar.gz ("RSiena" source code)

 9) R-4.1.2-arm64.pkg (Apple M1 Mac OSX version of R)

 10) SAOM_Model_1_Conventional_Model_Lambda_36.RData (SAOM/RSIENA Fitted Model Object - M1)
 11) SAOM_Model_2_Extended_Model_Lambda_36.RData (SAOM/RSIENA Fitted Model Object - M2)
 12) SAOM_Model_3_Network_Aid_Model_Restricted_Lambda_36.RData (SAOM/RSIENA Fitted Model Object - M3)
 13) SAOM_Model_4_Network_Aid_Model_Full_Lambda_36.RData (SAOM/RSIENA Fitted Model Object - M4)
 14) SAOM_Model_5_Conventional_Model_Lambda_108.RData (SAOM/RSIENA Fitted Model Object - M5)
 15) SAOM_Model_6_Extended_Model_Lambda_108.RData (SAOM/RSIENA Fitted Model Object - M6)
 16) SAOM_Model_7_Network_Aid_Model_Restricted_Lambda_108.RData (SAOM/RSIENA Fitted Model Object - M7)
 17) SAOM_Model_8_Network_Aid_Model_Full_Lambda_108.RData (SAOM/RSIENA Fitted Model Object - M8)
 
 18) Figure 1 SciData Social Support and Network Formation (2022-03-02).pdf
 19) Figure 1 SciData Social Support and Network Formation (2022-03-02).svg
 20) Reproducible Visone Visualisation Directions.txt (Directions to recreate Figure 1)
 21) Relatedness Custom Colour Palette KARPFENBLAU_GOLD.pdf (Colour palette used to create Figure 1)
 22) Arang.Dak.Tangible.Support.Intersection.graphml (Network file used to plot Figure 1 with Visone [https://visone.ethz.ch])



## Key Citations for Replicators
[1] Simpson, C.R. Accepted. "Social Support and Network Formation in a Small-Scale Horticulturalist Population". Scientific Data.

[2] Koster, Jeremy. (2018). "Family Ties: The Multilevel Effects of Households and Kinship on the Networks of Individuals." Royal Society Open Science 5(4):172159. https://royalsocietypublishing.org/doi/10.1098/rsos.172159

[3] Ripley, R.M., Snijders, T.A.B., Boda, Z., Vörös, A., Preciado, P., 2021. Manual for RSiena (v. 1.3.3). University of Oxford and University of Groningen. Available at: http://www.stats.ox.ac.uk/~snijders/siena/RSiena_Manual.pdf



## Notes
1) Thank you for your interest in my work! Please do let me know if something goes wrong. I am more than happy to help and you can always email me.

2) The estimation of the SAOMs takes a very, very long time. Although one can estimate these models more quickly (i.e., hours) using multiple CPU cores, doing so ignores the random seed and thus prevents exact reproducibility. Depending on the strength of your CPU, redoing my analyses is expected to take at least/around four days — primarily owing to the SAOMs that allow a larger number of average changes to network members' portfolios of outgoing ties during the simulation (for details, see the manuscript for discussion of the SAOM parameter "Lambda"). Note that the four-day expectation is optimistic as it is based on how long SAOMs take to estimate and completely converge using the Apple M1 Max Chip which, at the time of writing, has an exceptional single-core speed (https://browser.geekbench.com/macs/macbook-pro-16-inch-2021-apple-m1-max). Slower CPUs will see longer estimation times.


