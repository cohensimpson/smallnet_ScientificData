################################# Social Support and Network Formation in a Small-Scale Horticulturalist Population
################################# Replication Code


## To help ensure the reproducibility of findings from this project, I use the R package "groundhog"
library(groundhog)


## The groundhog.library() command requires two values. Like library(), you indicate which package
## you want to load. And then you must enter a date — any date (formatted as "yyyy-mm-dd"). 
## Groundhog will load the most recent version of the named package from CRAN corresponding to the  
## entered date. Groundhog will also load all dependencies of that package that are current on the
## entered date. Put simply, and borrowing from the package creator (see: https://groundhogr.com ):
   
## 1) Groundhog Makes R scripts reproducible by "replacing library("pkg")" with "groundhog.library("pkg", "date")"
## 2) groundhog.library() loads a package and its dependencies, as available on the chosen date.
## 3) Packages and their dependencies get automatically installed on their initial loading via groundhog if needed.
## 4) You will need to create a dedicated folder for groundhog package installs during the initial executing of "library(groundhog)".
## 4) Installation keeps, rather than replaces, existing other versions of a package of interest (e.g., versions from other dates).
## 5) If you want to stop using groundhog at any time, simply: replace "groundhog.library("pkg", "date")" with library("pkg") when loading a package of interest.

## Note that versions of RSiena are inconsistently pushed to CRAN (which groundhog pulls from). 
## Accordingly, I have included the source code for the version of RSiena used for this project
## with the other replication materials on the Open Science Framework. You will have to
## install RSiena manually (see: https://github.com/snlab-nl/rsiena/releases ). 
## Similarly, I have included the source code for the version of groundhog used for this 
## project with the replication materials and it should also be installed manually first.

## Note that Groundhog does not always play nicely with one's normal R library.
## Specifically, you may need to uninstall packages prior to use groundhog.library()
## Accordingly, please pay attention to the prompnts/warnings provided by groundhog

## N.B., you may need to give Groundhog permission to create a folder on your hard drive (follow the prompts)
## N.B., RStudio can preload packages in the background when using the package operators "::" and ":::" (https://community.rstudio.com/t/when-does-r-studio-autoload-packages/96649/2). This can break groundhog loading via package clashes.
## N.B., you may need to install GFortran which is used to install the Matrix package
## N.B., you may need to install the package ‘tcltk’ which RSiena uses to display model progress when estimation is not set to silent.

groundhog.library("network", "2022-01-31", quiet.install = TRUE, force.source = TRUE)
groundhog.library("sna", "2022-01-31", quiet.install = TRUE, force.source = TRUE)

groundhog.library("abind", "2022-01-31", quiet.install = TRUE, force.source = TRUE)
groundhog.library("dplyr", "2022-01-31", quiet.install = TRUE, force.source = TRUE) 
groundhog.library("purrr", "2022-01-31", quiet.install = TRUE, force.source = TRUE)

groundhog.library("stargazer", "2022-01-31", quiet.install = TRUE, force.source = TRUE)
groundhog.library("pastecs", "2022-01-31", quiet.install = TRUE, force.source = TRUE)

library(RSiena) 

library(parallel) ## Not on CRAN, automatically installed with Base R.



set.seed(20180709)
options(scipen = 8)
options(digits = 5)
options(max.print = 5000000)



#################################### SET NUMBER OF AVAILABLE COMPUTING CORES ####################################  
## Unfortunately, RSiena ignores the "seed" argument for the purposes of random number generation
## when estimating SAOMs using multiple CPU cores. In order for the "seed" argument below in sienaAlgorithmCreate()
## to not be ignored, siena07() — i.e., the function used to estimate SAOMs — must be run using a single core/not in 
## parallel. This slows estimation considerably. However, adherence to the random seed (here, 20180709) is necessary 
## to exactly reproduce results. 

cores <- 1



####################################  LOAD NETWORK AND ATTRIBUTE DATA  ####################################
###### N.B. The data are arranged to reflect resource flows. Accordingly, given the arc (i.e., asymmetric relation) X_ij, actor i is the resident that provides tangible support to j!
###### More specifically, and quoting Koster (2018, p. 6 and 8), the dataset documents responses to two sociometric questions: 
###### 1) "Who provides tangible support to you at least once per month?"  
###### 2) "To whom do you provide tangible support to at least once per month?"
###### Responses to the first question, which are given in the dataset in the column labelled “y”, is actor j’s report on whether actor i provides assistance to j (i.e., the "support seeking" ties). 
###### Responses to the the second, which are given in the dataset in the column labelled “y.donor.oriented,” is actor i's report on whether actor i provides assistance to j (i.e., the "support giving" ties).
###### Koster, J.M., 2018. Family ties: The multilevel effects of households and kinship on the networks of individuals. Royal Society Open Science, 5(4):172159. https://doi.org/10.1098/rsos.172159

#### Load the complete dataset. N.B.: Cases/rows are for *all* possible asymmetric relationships between the 108 residents of Arang Dak.
arang.dak.data <- read.csv("RSOS_corrected_data (14 May 2018).csv", header = TRUE, stringsAsFactors = FALSE)
arang.dak.data$i_ID_number <- arang.dak.data$i_ID ## Re-assign the original resident ID number of i which will be used to create a more informative unique label
arang.dak.data$j_ID_number <- arang.dak.data$j_ID ## Re-assign the original resident ID number of j which will be used to create a more informative unique label

#### Rearrange the data frame
#### N.B.: Each actor i belongs to a household k and each actor j belongs to a household l
arang.dak.data <- arang.dak.data[c("y","y.donor.oriented","i_ID","j_ID","i_ID_number","j_ID_number","ij_dyad_ID","k_ID","l_ID","kl_dyad_ID","il_ID","jk_ID",
                                   "i_age","j_age","i_sex","j_sex","i_bmi","j_bmi","i_skin","j_skin","k_wealth","l_wealth","ij_god_relation","ij_deg_r","ij_affinal_r",
                                   "kl_avg_interhouse_R","kl_sd_interhouse_R","kl_distance","kl_affair","i_avg_community_r","j_avg_community_r")]
arang.dak.data$i_ID <- paste0(  paste0("R", arang.dak.data$i_ID ), paste0("H", arang.dak.data$k_ID) ) ## Using a resident's ID number and their house ID number, create a unique ID to match their ID in the attribute data
arang.dak.data$j_ID <- paste0(  paste0("R", arang.dak.data$j_ID ), paste0("H", arang.dak.data$l_ID) ) ## Using a resident's ID number and their house ID number, create a unique ID to match their ID in the attribute data


#### Load the data frame containing all attribute variables for the 108 respondents
attributes <- read.csv("RSOS_corrected_data_attributes (17 May 2018).csv", header = TRUE, stringsAsFactors = FALSE)
rownames(attributes) <- attributes$ID ## These match the i_ID/j_ID created above

attributes$sex <- ifelse(attributes$sex == "F", 1, 0) ## 1 == Female
attributes$ethnicity <- ifelse(attributes$ethnicity == "Mis", 1, 0) ## 1 == Miskito
attributes$out_of_town <- ifelse(attributes$ID %in% c("R85H27", "R94H9"), 1, 0) ## These two residents had absconded from the community to have an affair
attributes$HH.size <- table(attributes$houseID)[attributes$houseID]

## Intra-village ranking of households by their wealth. 
## N.B. — the ranking needs to be done for the 32 households — not for all 108 villagers — as villagers who live together will of course have the same wealth and thus need to have the same rank.
## Rankings are then assigned to the 108 villagers based on their household ID.
wealth_hh_rank <- attributes[,c("houseID", "wealth_hh")][!duplicated(attributes[,c("houseID", "wealth_hh")]),]
wealth_hh_rank$wealth_hh_rank <- rank(wealth_hh_rank$wealth_hh)

HH.wealth.rank <- wealth_hh_rank$wealth_hh_rank
names(HH.wealth.rank) <- wealth_hh_rank$houseID

attributes$wealth_hh_rank <- as.numeric(HH.wealth.rank[attributes$houseID]) ## Index HH.wealth.rank for each houseID in attributes

rm(wealth_hh_rank, HH.wealth.rank)


#### Extract edgelist for the "support seeking" ties amongst the 108 respondents from the complete dataset
#### Support Seeking Tie == "y" == 1 == i provides tangible support to j **according to the seeker of aid j**
tangible_support.edge.list.seeking <- arang.dak.data[c("y", "i_ID", "j_ID" ,"i_ID_number", "j_ID_number")]
tangible_support.edge.list.seeking <- subset(tangible_support.edge.list.seeking, tangible_support.edge.list.seeking$y == 1)


#### Extract edgelist for "support giving" ties amongst the 108 respondents from the complete dataset
#### Support Giving Tie == "y.donor.oriented" == 1 == i provides tangible support to j **according to the giver of aid i**
tangible_support.edge.list.giving <- arang.dak.data[c("y.donor.oriented", "i_ID", "j_ID" ,"i_ID_number", "j_ID_number")]
tangible_support.edge.list.giving <- subset(tangible_support.edge.list.giving, tangible_support.edge.list.giving$y.donor.oriented == 1)



################################ Tangible Support Seeking Network ################################ 
tangible.support.matrix.seeking <- matrix(data = 0, nrow = nrow(attributes), ncol = nrow(attributes))
colnames(tangible.support.matrix.seeking) <- attributes$ID
rownames(tangible.support.matrix.seeking) <- attributes$ID

## Add support seeking ties to tangible.support.matrix.seeking
for(i in 1:nrow(tangible_support.edge.list.seeking)){
  
  source <- as.character( tangible_support.edge.list.seeking$i_ID[i] ) ## "as.character()" is a sanity check
  target <- as.character( tangible_support.edge.list.seeking$j_ID[i] )
  
  ## Index the matrix of zeros by its row names and column names (i.e., the resident IDs) to add the support seeking tie (i.e., a one) in the appropriate cell
  tangible.support.matrix.seeking[which( rownames(tangible.support.matrix.seeking) == source ), which( colnames(tangible.support.matrix.seeking) == target )] <- 1
  
  rm(i, source, target)
}

sum(tangible.support.matrix.seeking)



################################ Tangible Support Giving Network ################################ 
tangible.support.matrix.giving <- matrix(data = 0, nrow = nrow(attributes), ncol = nrow(attributes))
colnames(tangible.support.matrix.giving) <- attributes$ID
rownames(tangible.support.matrix.giving) <- attributes$ID

## Add support giving ties to tangible.support.matrix.giving
for(i in 1:nrow(tangible_support.edge.list.giving)){
  
  source <- as.character( tangible_support.edge.list.giving$i_ID[i] )
  target <- as.character( tangible_support.edge.list.giving$j_ID[i] )
  
  ## Index the matrix of zeros by its row names and column names (i.e., the resident IDs) to add the support giving tie (i.e., a one) in the appropriate cell
  tangible.support.matrix.giving[which( rownames(tangible.support.matrix.giving) == source ), which( colnames(tangible.support.matrix.giving) == target )] <- 1
  
  rm(i, source, target)
}

sum(tangible.support.matrix.giving)



################################ Inter-household Distances ################################ 
geodist.edge.list <- arang.dak.data[c("i_ID", "j_ID", "kl_distance")] ## Grab the inter-household distances for all possible arcs amongst the 108 residents

geographic.distance <- matrix(data = 0, nrow = nrow(attributes), ncol = nrow(attributes))
colnames(geographic.distance) <- attributes$ID
rownames(geographic.distance) <- attributes$ID


for(i in 1:nrow(geodist.edge.list)){
  
  source <- as.character( geodist.edge.list$i_ID[i] ) ## i_ID is a member of household k
  target <- as.character( geodist.edge.list$j_ID[i] ) ## j_ID is a member of household l
  
  inter_household_distance <- geodist.edge.list$kl_distance[i]
  
  ## Index the matrix of zeros by its row names and column names (i.e., the resident IDs) to add the inter-household distance in the appropriate cell
  geographic.distance[which( rownames(geographic.distance) == source ), which( colnames(geographic.distance) == target )] <- inter_household_distance
  
  rm(i, source, target, inter_household_distance)
}

table(geographic.distance == t(geographic.distance)) ## Sanity Check: Compare constructed matrix to its transposed version to ensure that the matrix is symmetric



################################ Consanguineal Degree of Relatedness ################################ 
relatedness.edge.list <- arang.dak.data[c("i_ID", "j_ID", "ij_deg_r")] ## Grab Wright's coefficient of relatedness (i.e., consanguineal relatedness) for all possible arcs amongst the 108 residents

degree.of.relatedness <- matrix(data = 0, nrow = nrow(attributes), ncol = nrow(attributes))
colnames(degree.of.relatedness) <- attributes$ID
rownames(degree.of.relatedness) <- attributes$ID


for(i in 1:nrow(relatedness.edge.list)){
  
  source <- as.character( relatedness.edge.list$i_ID[i] )
  target <- as.character( relatedness.edge.list$j_ID[i] )
  
  relatedness <- relatedness.edge.list$ij_deg_r[i]
  
  ## Index the matrix of zeros by its row names and column names (i.e., the resident IDs) to add the consanguineal degrees of relatedness in the appropriate cell
  degree.of.relatedness[which( rownames(degree.of.relatedness) == source ), which( colnames(degree.of.relatedness) == target )] <- relatedness
  
  rm(i, source, target, relatedness)
}

table(degree.of.relatedness == t(degree.of.relatedness)) ## Sanity Check: Compare constructed matrix to its transposed version to ensure that the matrix is symmetric



################################ Affinal Degree of Relatedness ################################
affinal.relatedness.edge.list <- arang.dak.data[c("i_ID", "j_ID", "ij_affinal_r")] ## Grab Wright's coefficient of relatedness between i's spouse s and the alter j (i.e., affinal relatedness) for all possible arcs amongst the 108 residents

affinal.degree.of.relatedness <- matrix(data = 0, nrow = nrow(attributes), ncol = nrow(attributes))
colnames(affinal.degree.of.relatedness) <- attributes$ID
rownames(affinal.degree.of.relatedness) <- attributes$ID


for(i in 1:nrow(affinal.relatedness.edge.list)){
  
  source <- as.character( affinal.relatedness.edge.list$i_ID[i] )
  target <- as.character( affinal.relatedness.edge.list$j_ID[i] )
  
  affinal.relatedness <- affinal.relatedness.edge.list$ij_affinal_r[i]
  
  ## Index the matrix of zeros by its row names and column names (i.e., the resident IDs) to add the affinal degree of relatedness in the appropriate cell
  affinal.degree.of.relatedness[which( rownames(affinal.degree.of.relatedness) == source ), which( colnames(affinal.degree.of.relatedness) == target )] <- affinal.relatedness
  
  rm(i, source, target, affinal.relatedness)
}

table(affinal.degree.of.relatedness == t(affinal.degree.of.relatedness)) ## Sanity Check: Compare constructed matrix to its transposed version to ensure that the matrix is symmetric



################################ Fictive Relatedness: Godparental Ties ################################
fictive.relatedness.edge.list <- arang.dak.data[c("i_ID", "j_ID", "ij_god_relation")] ## Grab godparental ties for all possible arcs amongst the 108 residents

fictive.relatedness <- matrix(data = 0, nrow = nrow(attributes), ncol = nrow(attributes))
colnames(fictive.relatedness) <- attributes$ID
rownames(fictive.relatedness) <- attributes$ID


for(i in 1:nrow(fictive.relatedness.edge.list)){
  
  source <- as.character( fictive.relatedness.edge.list$i_ID[i] )
  target <- as.character( fictive.relatedness.edge.list$j_ID[i] )
  
  godparental.tie <- fictive.relatedness.edge.list$ij_god_relation[i]
  
  ## Index the matrix of zeros by its row names and column names (i.e., the resident IDs) to add the godparental tie in the appropriate cell
  fictive.relatedness[which( rownames(fictive.relatedness) == source ), which( colnames(fictive.relatedness) == target )] <- godparental.tie
  
  rm(i, source, target, godparental.tie)
}

table(fictive.relatedness == t(fictive.relatedness)) ## Sanity Check: Compare constructed matrix to its transposed version to ensure that the matrix is symmetric



################################ Infidelity Relations: Ties via Illegitimate Children ################################
infidelity.relatedness.edge.list <- arang.dak.data[c("i_ID", "j_ID", "kl_affair")] ## Grab ties via illegitimate children for all possible arcs amongst the 108 residents

infidelity.relatedness <- matrix(data = 0, nrow = nrow(attributes), ncol = nrow(attributes))
colnames(infidelity.relatedness) <- attributes$ID
rownames(infidelity.relatedness) <- attributes$ID


for(i in 1:nrow(infidelity.relatedness.edge.list)){
  
  source <- as.character( infidelity.relatedness.edge.list$i_ID[i] )
  target <- as.character( infidelity.relatedness.edge.list$j_ID[i] )
  
  infidelity.tie <- infidelity.relatedness.edge.list$kl_affair[i]
  infidelity.tie <- ifelse(infidelity.tie == "affair", 1, 0) ## Recode the character entries for fitting the models; "affair" == Yes == 1; "ref" == No == 0
  
  ## Index the matrix of zeros by its row names and column names (i.e., the resident IDs) to add the affair-related tie in the appropriate cell
  infidelity.relatedness[which( rownames(infidelity.relatedness) == source ), which( colnames(infidelity.relatedness) == target )] <- infidelity.tie
  
  rm(i, source, target, infidelity.tie)
}

table(infidelity.relatedness == t(infidelity.relatedness)) ## Sanity Check: Compare constructed matrix to its transpose to ensure that the matrix is symmetric



######################################################################## STOCHASTIC ACTOR-ORIENTED MODELS (SAOMs) ########################################################################
############################################## FUNCTION FOR THE ITERATIVE RUNNING OF siena07() ############################################## 
### The following function will repeatedly execute siena07() until a SAOM has converged in line with the convergence criteria below.
### Here, prevAns (i.e., an earlier existing "on track" estimation result) is repeatedly used to determine the initial values for estimation at each successive iteration of the algorithm.
### Because each run of siena07() can take many hours, each iteration is saved as an .RData file so no progress is lost.
siena07RunToConvergence <- function(alg, dat, eff, ans0, modelName, ...){
  
  numr <- 0
  
  ans <- siena07(alg, data = dat, effects = eff, prevAns = ans0, returnDeps = TRUE, ...)
  
  repeat{
    numr <- numr + 1  ## Count the number of repeated runs
    tconv.max <- ans$tconv.max  ## Extract the overall maximum convergence ratio
    tratio.max <- max( abs( ans$tstat[eff$type[eff$include] != "rate"] ) ) ## Extract the maximum absolute value of the convergence t-ratios. Don't include the t-ratio for the rate parameter as it is fixed!
    
    if (tconv.max > 100) {  ## Divergence without much hope of returning to good parameter values
      print(ans)
      cat("WARNING: Extreme Divergence. Terminating run.\n")
      return("WARNING: Extreme Divergence. Terminating run")
    }
    ## These are the convergence criteria used for the study. Convergence is excellent when the overall maximum convergence ratio (tconv.max) 
    ## is less than 0.15, and, for all the individual parameters, the t-ratios for convergence (tratio.max) are all less than 0.1 in absolute value
    else if (tconv.max < 0.15 & tratio.max < 0.10) {
      print(ans)
      cat(paste0("Maximum Absolute Value Amongst Convergence t-Ratios: ", tratio.max, "\n"))
      cat(paste0("Model Has Converged After ", numr, " iterations. \n"))
      
      save(ans, file =  paste0(modelName,"_SIENA_Iteration_Number_", numr,"_CONVERGED.RData") )
      
      return(ans)
      
    }
    else {
      print(ans)
      cat("WARNING: Convergence Inadequate.\n")
      cat(paste0("Overall Maximum Convergence Ratio: ", tconv.max, "\n"))
      cat(paste0("Iteration Number: ", numr), "\n") ## Report how far along we are
      
      save(ans, file =  paste0(modelName,"_SIENA_Iteration_Number_", numr,"_NOT_CONVERGED.RData") )
      
      ans <- siena07(alg, data = dat, effects = eff, prevAns = ans, returnDeps = TRUE, ...)
      
    }
  }
}



################################# Network Descriptive Statistics #################################
gden(list(tangible.support.matrix.seeking, tangible.support.matrix.giving), mode = "digraph") ## Network Density
gcor(tangible.support.matrix.seeking, tangible.support.matrix.giving)
table(tangible.support.matrix.seeking, tangible.support.matrix.giving)


## For my analysis, I follow Lee and Butts (2018) as opposed to Kasper and Borgerhoff Mulder (2015) and Nolin (2010) who assume that people are "honest but forgetful"
## Lee, F., Butts, C.T., 2018. Mutual assent or unilateral nomination? A performance comparison of intersection and union rules for integrating self-reports of social relationships. Social Networks, 55:55–62. https://doi.org/10.1016/j.socnet.2018.05.005
## Kasper, C., Borgerhoff Mulder, M., 2015. Who helps and why? Cooperative Networks in Mpimbwe. Current Anthropology, 56(5):701–732. https://doi.org/10.1086/683024
## Nolin, D.A., 2010. Food-sharing networks in Lamalera, Indonesia. Human Nature, 21(3):243–268. https://doi.org/10.1007/s12110-010-9091-3
tangible.support.matrix.intersection <- tangible.support.matrix.seeking + tangible.support.matrix.giving 
tangible.support.matrix.intersection[tangible.support.matrix.intersection < 2] <- 0
tangible.support.matrix.intersection[tangible.support.matrix.intersection > 0] <- 1 ## tangible.support.matrix.intersection == "source-recipient/target-varied network" mentioned in the paper


## The two out_of_town individuals provided no data on who gave them tangible support (i.e., missing column/receiver info). Accordingly "impute" their incoming ties using the unilateral reports of the available aid donors.
tangible.support.matrix.intersection[, which(attributes$out_of_town == 1)] <- tangible.support.matrix.giving[, which(attributes$out_of_town == 1)] 


## The two out_of_town individuals provided no data on to whom they gave tangible support (i.e., missing row/sender info). Accordingly, "impute" their outgoing ties using the unilateral reports of the available aid seekers.
tangible.support.matrix.intersection[which(attributes$out_of_town == 1), ] <- tangible.support.matrix.seeking[which(attributes$out_of_town == 1), ]  


gden(tangible.support.matrix.intersection, mode = "digraph") ## Network Density
grecip(tangible.support.matrix.intersection, measure = "edgewise") ## Tie Reciprocity (i.e., proportion of edges/arcs which are reciprocated)
gtrans(tangible.support.matrix.intersection, mode = "digraph", measure = "weak") ## Graph-Level Transitivity

dyad.census(tangible.support.matrix.intersection)
triad.census(tangible.support.matrix.intersection)

quantile(degree(tangible.support.matrix.intersection, gmode = "digraph", cmode = "outdegree"), probs = seq(0, 1 ,0.05))
quantile(degree(tangible.support.matrix.intersection, gmode = "digraph", cmode = "indegree"), probs = seq(0, 1, 0.05))

table(degree(tangible.support.matrix.intersection, gmode = "digraph", cmode = "outdegree"))
table(degree(tangible.support.matrix.intersection, gmode = "digraph", cmode = "indegree"))

stat.desc(degree(tangible.support.matrix.intersection, gmode = "digraph", cmode = "outdegree"))
stat.desc(degree(tangible.support.matrix.intersection, gmode = "digraph", cmode = "indegree"))



################################ PREPARATION OF OBJECTS FOR MODEL FITTING #################################
villagers <- rownames(tangible.support.matrix.intersection)
villagers.size <- length(villagers)
villagers <- sienaNodeSet(villagers.size, nodeSetName = "villagers", names = villagers)



#################################### THE DEPENDENT NETWORK
### RSiena does not run when the observed waves of the analysed network are identical. 
### Therefore, the second wave is slightly modified by randomly changing one tie in the adjacency matrix.
### This directive for estimating the cross-sectional SAOM with RSiena comes directly from Snijders, T. A. B., & Steglich, C. E. G. (2015). 
### "Representing Micro–Macro Linkages by Actor-based Dynamic Network Models." Sociological Methods & Research, 44(2), 222–271. http://doi.org/10.1177/0049124113494573

### First we need a 3D array of network snapshots/waves where the SAME wave is used along the third dimension 
support_array <- array(data = c( tangible.support.matrix.intersection, tangible.support.matrix.intersection ), 
                          dim = c(villagers.size, villagers.size, 2) )


#### And second, we need to add a random tie to the "second" network wave but don't add a self-loop!
villagers.to.which.to.add.a.random.tie <- sample( c(1:villagers.size), 2, replace = FALSE)

support_array[villagers.to.which.to.add.a.random.tie[1], villagers.to.which.to.add.a.random.tie[2], 2] <- 1


#### Finally, formally create the SIENA dependent network object
support_net <- sienaNet(support_array, type = "oneMode", nodeSet = "villagers"
                           , allowOnly = FALSE)



#################################### MONADIC COVARIATES — CENTRED AT THEIR WITHIN-VILLAGE MEANS
log_wealth_hh_Z <- coCovar(scale( log( attributes$wealth_hh ), center = TRUE, scale = TRUE)[,1],
                           centered = FALSE, nodeSet = "villagers")


HH.size_Z <- coCovar(scale( attributes$HH.size, center = TRUE, scale = TRUE)[,1],
                     centered = FALSE, nodeSet = "villagers")


age_Z <- coCovar(scale( attributes$age, center = TRUE, scale = TRUE)[,1],
                 centered = FALSE, nodeSet = "villagers")


gender <- coCovar(attributes$sex, ## Female = 1; Male = 0
                  centered = FALSE, nodeSet = "villagers") 


ethnicity <- coCovar(attributes$ethnicity, ## Miskito = 1; Mayangna = 0
                     centered = FALSE, nodeSet = "villagers")


melanin_index_Z <- coCovar(scale( attributes$melanin_index, center = TRUE, scale = TRUE)[,1],
                           centered = FALSE, nodeSet = "villagers")


bmi_Z <- coCovar(scale( attributes$bmi, center = TRUE, scale = TRUE)[,1],
                           centered = FALSE, nodeSet = "villagers")


HH.ID <- coCovar(attributes$houseID, ## Household Number
                 centered = FALSE, nodeSet = "villagers")



#################################### DYADIC COVARIATES — NON-CENTRED
consanguineal.relatedness <- coDyadCovar( degree.of.relatedness, ## Wright's coefficient of relatedness
                         centered = FALSE, nodeSets = c("villagers", "villagers"),
                         type = "oneMode")


affinal.relatedness <- coDyadCovar( affinal.degree.of.relatedness, ## Wright's coefficient of relatedness through marriage
                                 centered = FALSE, nodeSets = c("villagers", "villagers"),
                                 type = "oneMode")


godparental.relation <- coDyadCovar( fictive.relatedness, ## Godparental Tie
                                 centered = FALSE, nodeSets = c("villagers", "villagers"),
                                 type = "oneMode")


infidelity.relation <- coDyadCovar( infidelity.relatedness, ## Inter-household illegitimate child tie
                                    centered = FALSE, nodeSets = c("villagers", "villagers"),
                                    type = "oneMode")


## log of the geographic distance following: Preciado, P., Snijders, T.A., Burk, W.J., Stattin, H., Kerr, M., 2012. Does proximity matter? Distance dependence of adolescent friendships. Social Networks, 34(1):18-31. https://doi.org/10.1016/j.socnet.2011.01.002
geodist.dyad <- coDyadCovar( log(geographic.distance + 1) ,
                             centered = FALSE, nodeSets = c("villagers", "villagers"),
                             type = "oneMode")


relative.wealth_hh.rank <- data.frame()
for(i in 1:nrow(attributes)){
  wealth_hh.rank <- attributes$wealth_hh_rank[i] ## What is the ranking of the household wealth of i in Arang Dak?
  relative.rank <- wealth_hh.rank - attributes$wealth_hh_rank ## What is the difference between i's wealth rank and the wealth ranks of all other residents in Arang Dak?

  relative.wealth_hh.rank <- rbind(relative.wealth_hh.rank, relative.rank)

  rm(i, wealth_hh.rank, relative.rank)
}
relative.wealth_hh.rank <- as.matrix(relative.wealth_hh.rank)
colnames(relative.wealth_hh.rank) <- attributes$ID
rownames(relative.wealth_hh.rank) <- attributes$ID
diag(relative.wealth_hh.rank) <- 0

relative.wealth_hh.rank <- coDyadCovar( relative.wealth_hh.rank, ## N.B. Do not take the absolute value as this will obfuscate the directionality of need! See commentary by Koster in response to Kasper and Borgerhoff Mulder (2015, p. 720).
                             centered = FALSE, nodeSets = c("villagers", "villagers"),
                             type = "oneMode")


## Cannot assess internal SIENA-created interaction effects with sienaRI(). Accordingly, do the dyadic interaction by hand.
consanguineal.relatedness_x_relative.wealth_hh.rank <- coDyadCovar( degree.of.relatedness*relative.wealth_hh.rank, ## Wright's coefficient of relatedness X relative wealth rank 
                                          centered = FALSE, nodeSets = c("villagers", "villagers"),
                                          type = "oneMode")


## Cannot assess internal SIENA-created interaction effects with sienaRI(). Accordingly, do the dyadic interaction by hand.
consanguineal.relatedness_x_geodist.dyad <- coDyadCovar( degree.of.relatedness*log(geographic.distance + 1), ## Wright's coefficient of relatedness X geographic distance
                                          centered = FALSE, nodeSets = c("villagers", "villagers"),
                                          type = "oneMode")



#################################### CREATE THE RSIENA DATA OBJECT FOR MODEL FITTING
multidata <- sienaDataCreate(support_net, 
                             relative.wealth_hh.rank,
                             geodist.dyad,
                             consanguineal.relatedness,
                             consanguineal.relatedness_x_relative.wealth_hh.rank,
                             consanguineal.relatedness_x_geodist.dyad,
                             affinal.relatedness, 
                             godparental.relation, 
                             infidelity.relation,
                             log_wealth_hh_Z,
                             HH.size_Z,
                             age_Z, 
                             gender,
                             ethnicity,
                             melanin_index_Z,
                             bmi_Z,
                             HH.ID,
                             nodeSets = list(villagers)) 

#### RSIENA-provided summary statistics about the dependent network and covariates
# print01Report(multidata, modelname = "RI_Arang_Dak_2021_sienaDataSummary", getDocumentation = FALSE)



############################################# DEFINE SIENA ALGORITHM AND RUN CROSS-SECTIONAL SAOMS
modelparams <- sienaAlgorithmCreate(projname = "RI_Arang_Dak_2022_Estimation_History", cond = FALSE, maxlike = FALSE
                                    
                                    # Number of subphases in phase 2.
                                    , nsub = 4
                                    
                                    # Number of iterations in phase 3. For regular use with the Method of Moments, n3 = 1000
                                    # mostly suffices. For use in publications and for Maximum Likelihood, at least n3 = 3000
                                    # is advised. Sometimes much higher values are required for stable estimation of standard errors.
                                    , n3 = 20000
                                    
                                    # This determines the step sizes in the estimation algorithm. If the algorithm is unstable 
                                    # (e.g., oscillating between wild parameter estimates and convergence from run to run),
                                    # use a smaller value (but greater than 0). The default value is 0.2. Sometimes
                                    # for difficult data-model combinations, the algorithm diverges very quickly, and this
                                    # may be countered by smaller values of firstg, e.g., 0.01 or 0.05.
                                    , firstg = 0.2
                                    
                                    # Number between 0 and 1 (bounds included), values outside this interval will be truncated; 
                                    # for diagonalize = 0 the complete estimated derivative matrix will be used for updates in the Robbins-Monro procedure; 
                                    # for diagonalize = 1 only the diagonal entries will be used; for values between 0 and 1, the weighted average will be 
                                    # used with weight diagonalize for the diagonalized matrix. Has no effect for Maximum Likelihood estimation.
                                    # Higher values are more stable, lower values potentially more efficient. 
                                    # Default for Method of Moments estimation is diagonalize = 0.2.
                                    , diagonalize = 0.2
                                    
                                    # The random seed will NOT be adhered to if one runs siena07RunToConvergence()/siena07() with multiple cores
                                    , seed = 20180709
)



##### Model 1: The Standard/"Conventional"/"Human Behavioural Ecology Model" with Lambda = 36 (i.e., the maximum out-degree in support_net)
##### Typically converges after one iteration of siena07RunToConvergence
fit1.modeffects <- getEffects(multidata)
fit1.modeffects <- setEffect(fit1.modeffects, Rate, initialValue = 36, name = "support_net", type = "rate", fix = TRUE, verbose = FALSE) ## Explicitly fix the rate parameter lambda at the maximum observed out-degree in the intersection/source-target verified network

fit1.modeffects <- includeEffects(fit1.modeffects, recip, name = "support_net", type = "eval", fix = FALSE, verbose = FALSE, include = TRUE) ## Contingent giving; N.B., the reciprocity effect is automatically included in SAOMs. Make it Explicit.
fit1.modeffects <- includeEffects(fit1.modeffects, X, name = "support_net", interaction1 = "geodist.dyad", type = "eval", fix = FALSE, verbose = FALSE)
fit1.modeffects <- includeEffects(fit1.modeffects, X, name = "support_net", interaction1 = "consanguineal.relatedness", type = "eval", fix = FALSE, verbose = FALSE) ## Un-restricted giving to genetic kin
fit1.modeffects <- includeEffects(fit1.modeffects, X, name = "support_net", interaction1 = "affinal.relatedness", type = "eval", fix = FALSE, verbose = FALSE) ## Un-restricted giving to marriage-based kin


fit1.ans <- siena07RunToConvergence(alg = modelparams, dat = multidata, eff = fit1.modeffects, ans0 = NULL, modelName = "fit1.ans", batch = TRUE, verbose = FALSE, silent = FALSE, nbrNodes = cores, useCluster = FALSE)



##### Model 2: The "Extended" Anthropology Model with Lambda = 36 
##### Typically converges after two iterations of siena07RunToConvergence
fit2.modeffects <- includeEffects(fit1.modeffects, X, name = "support_net", interaction1 = "relative.wealth_hh.rank", type = "eval", fix = FALSE, verbose = FALSE) ## Need-Based Transfer following: Thomas, M.G., Ji, T., Wu, J., He, Q.Q., Tao, Y., Mace, R., 2018. Kinship underlies costly cooperation in Mosuo villages. Royal Society Open Science, 5(2):171535. https://doi.org/10.1098/rsos.171535 

fit2.modeffects <- includeEffects(fit2.modeffects, XRecip, name = "support_net", interaction1 = "consanguineal.relatedness", type = "eval", fix = FALSE, verbose = FALSE) ## Kin-favoured reciprocity (Kasper and Borgerhoff Mulder 2015)
fit2.modeffects <- includeEffects(fit2.modeffects, X, name = "support_net", interaction1 = "consanguineal.relatedness_x_relative.wealth_hh.rank", type = "eval", fix = FALSE, verbose = FALSE) ## Kin-directed altruism (Kasper and Borgerhoff Mulder 2015)
fit2.modeffects <- includeEffects(fit2.modeffects, X, name = "support_net", interaction1 = "consanguineal.relatedness_x_geodist.dyad", type = "eval", fix = FALSE, verbose = FALSE) ## Travel greater distances to help kin (Thomas et al. 2018)


fit2.modeffects <- includeEffects(fit2.modeffects, X, name = "support_net", interaction1 = "godparental.relation", type = "eval", fix = FALSE, verbose = FALSE) ## Un-restricted giving to fictive kin
fit2.modeffects <- includeEffects(fit2.modeffects, X, name = "support_net", interaction1 = "infidelity.relation", type = "eval", fix = FALSE, verbose = FALSE) ## Un-restricted giving to illicit genetic kin


fit2.modeffects <- includeEffects(fit2.modeffects, egoX, name = "support_net", interaction1 = "log_wealth_hh_Z", type = "eval", fix = FALSE, verbose = FALSE)
fit2.modeffects <- includeEffects(fit2.modeffects, altX, name = "support_net", interaction1 = "log_wealth_hh_Z", type = "eval", fix = FALSE, verbose = FALSE)

fit2.modeffects <- includeEffects(fit2.modeffects, egoX, name = "support_net", interaction1 = "HH.size_Z", type = "eval", fix = FALSE, verbose = FALSE) ## Adults Only
fit2.modeffects <- includeEffects(fit2.modeffects, altX, name = "support_net", interaction1 = "HH.size_Z", type = "eval", fix = FALSE, verbose = FALSE) ## Adults Only

fit2.modeffects <- includeEffects(fit2.modeffects, egoX, name = "support_net", interaction1 = "age_Z", type = "eval", fix = FALSE, verbose = FALSE)
fit2.modeffects <- includeEffects(fit2.modeffects, altX, name = "support_net", interaction1 = "age_Z", type = "eval", fix = FALSE, verbose = FALSE)
fit2.modeffects <- includeEffects(fit2.modeffects, simX, name = "support_net", interaction1 = "age_Z", type = "eval", fix = FALSE, verbose = FALSE)

fit2.modeffects <- includeEffects(fit2.modeffects, egoX, name = "support_net", interaction1 = "gender", type = "eval", fix = FALSE, verbose = FALSE)
fit2.modeffects <- includeEffects(fit2.modeffects, altX, name = "support_net", interaction1 = "gender", type = "eval", fix = FALSE, verbose = FALSE)
fit2.modeffects <- includeEffects(fit2.modeffects, sameX, name = "support_net", interaction1 = "gender", type = "eval", fix = FALSE, verbose = FALSE)

fit2.modeffects <- includeEffects(fit2.modeffects, egoX, name = "support_net", interaction1 = "ethnicity", type = "eval", fix = FALSE, verbose = FALSE)
fit2.modeffects <- includeEffects(fit2.modeffects, altX, name = "support_net", interaction1 = "ethnicity", type = "eval", fix = FALSE, verbose = FALSE)
fit2.modeffects <- includeEffects(fit2.modeffects, sameX, name = "support_net", interaction1 = "ethnicity", type = "eval", fix = FALSE, verbose = FALSE)

fit2.modeffects <- includeEffects(fit2.modeffects, egoX, name = "support_net", interaction1 = "melanin_index_Z", type = "eval", fix = FALSE, verbose = FALSE)
fit2.modeffects <- includeEffects(fit2.modeffects, altX, name = "support_net", interaction1 = "melanin_index_Z", type = "eval", fix = FALSE, verbose = FALSE)
fit2.modeffects <- includeEffects(fit2.modeffects, simX, name = "support_net", interaction1 = "melanin_index_Z", type = "eval", fix = FALSE, verbose = FALSE)

fit2.modeffects <- includeEffects(fit2.modeffects, egoX, name = "support_net", interaction1 = "bmi_Z", type = "eval", fix = FALSE, verbose = FALSE)
fit2.modeffects <- includeEffects(fit2.modeffects, altX, name = "support_net", interaction1 = "bmi_Z", type = "eval", fix = FALSE, verbose = FALSE)
fit2.modeffects <- includeEffects(fit2.modeffects, simX, name = "support_net", interaction1 = "bmi_Z", type = "eval", fix = FALSE, verbose = FALSE)


fit2.ans <- siena07RunToConvergence(alg = modelparams, dat = multidata, eff = fit2.modeffects, ans0 = NULL, modelName = "fit2.ans", batch = TRUE, verbose = FALSE, silent = FALSE, nbrNodes = cores, useCluster = FALSE)



##### Model 3: The Networked Aid Model (Limited) with Lambda = 36
##### Typically converges after three iterations of siena07RunToConvergence
fit3.modeffects <- includeEffects(fit2.modeffects, inPop, name = "support_net", type = "eval", fix = FALSE, verbose = FALSE) 
fit3.modeffects <- includeEffects(fit3.modeffects, outPop, name = "support_net", type = "eval", fix = FALSE, verbose = FALSE) 

fit3.modeffects <- includeEffects(fit3.modeffects, transTrip, name = "support_net", type = "eval", fix = FALSE, verbose = FALSE) 

fit3.ans <- siena07RunToConvergence(alg = modelparams, dat = multidata, eff = fit3.modeffects, ans0 = NULL, modelName = "fit3.ans", batch = TRUE, verbose = FALSE, silent = FALSE, nbrNodes = cores, useCluster = FALSE)



##### Model 4: The Networked Aid Model (Comprehensive) with Lambda = 36
##### Typically converges after three iterations of siena07RunToConvergence
fit4.modeffects <- includeEffects(fit3.modeffects, outAct, name = "support_net", type = "eval", fix = FALSE, verbose = FALSE) 

fit4.modeffects <- includeEffects(fit4.modeffects, transRecTrip, name = "support_net", type = "eval", fix = FALSE, verbose = FALSE) 
fit4.modeffects <- includeEffects(fit4.modeffects, cycle3, name = "support_net", type = "eval", fix = FALSE, verbose = FALSE) 

fit4.modeffects <- setEffect(fit4.modeffects, denseTriads, parameter = 6, name = "support_net", type = "eval", fix = FALSE, verbose = FALSE)

fit4.modeffects <- includeEffects(fit4.modeffects, jumpXTransTrip, name = "support_net", interaction1 = "HH.ID", type = "eval", fix = FALSE, verbose = FALSE) ## The organising role of households (Koster 2018)

fit4.modeffects <- includeEffects(fit4.modeffects, sharedPop, name = "support_net", type = "eval", fix = FALSE, verbose = FALSE)


fit4.ans <- siena07RunToConvergence(alg = modelparams, dat = multidata, eff = fit4.modeffects, ans0 = NULL, modelName = "fit4.ans", batch = TRUE, verbose = FALSE, silent = FALSE, nbrNodes = cores, useCluster = FALSE)



################################ MULTI-PARAMETER WALD TESTS (MODELS 1-4) ################################ 
## RUN: ?Multipar.RSiena
fit2.ans.Walt.test <- Multipar.RSiena(ans = fit2.ans, 4, 7:9, 11:31) ## Positive integers specify the tested effects (as numbered in "print(ans)")
print(fit2.ans.Walt.test)

fit3.ans.Walt.test <- Multipar.RSiena(ans = fit3.ans, 4:6) 
print(fit3.ans.Walt.test)

fit4.ans.Walt.test <- Multipar.RSiena(ans = fit4.ans, 5:8, 11, 40)
print(fit4.ans.Walt.test)


################################# MODEL ESTIMATION: ROBUSTNESS CHECK USING A TRIPLED (i.e., 36*3) LAMBDA VALUE ################################# 
##### Model 5: The Standard/"Conventional"/"Human Behavioural Ecology Model" with Lambda = 108
##### Typically converges after one iteration of siena07RunToConvergence
fit5.modeffects <- setEffect(fit1.modeffects, Rate, initialValue = 108, name = "support_net", type = "rate", fix = TRUE, verbose = FALSE) ## Explicitly fix the rate parameter lambda at the maximum observed out-degree (in the intersection/source-target verified network)

fit5.ans <- siena07RunToConvergence(alg = modelparams, dat = multidata, eff = fit5.modeffects, ans0 = NULL, modelName = "fit5.ans", batch = TRUE, verbose = FALSE, silent = FALSE, nbrNodes = cores, useCluster = FALSE)


##### Model 6: The "Extended" Anthropology Model with Lambda = 108
##### Typically converges after two iterations of siena07RunToConvergence
fit6.modeffects <- setEffect(fit2.modeffects, Rate, initialValue = 108, name = "support_net", type = "rate", fix = TRUE, verbose = FALSE) ## Explicitly fix the rate parameter lambda at the maximum observed out-degree (in the intersection/source-target verified network)

fit6.ans <- siena07RunToConvergence(alg = modelparams, dat = multidata, eff = fit6.modeffects, ans0 = NULL, modelName = "fit6.ans", batch = TRUE, verbose = FALSE, silent = FALSE, nbrNodes = cores, useCluster = FALSE)


##### Model 7: The Networked Aid Model (Limited) with Lambda = 108
##### Typically converges after three iterations of siena07RunToConvergence
fit7.modeffects <- setEffect(fit3.modeffects, Rate, initialValue = 108, name = "support_net", type = "rate", fix = TRUE, verbose = FALSE) ## Explicitly fix the rate parameter lambda at the maximum observed out-degree (in the intersection/source-target verified network)

fit7.ans <- siena07RunToConvergence(alg = modelparams, dat = multidata, eff = fit7.modeffects, ans0 = NULL, modelName = "fit7.ans", batch = TRUE, verbose = FALSE, silent = FALSE, nbrNodes = cores, useCluster = FALSE)


##### Model 8: The Networked Aid Model (Comprehensive) with Lambda = 108
##### Typically converges after five iterations of siena07RunToConvergence
fit8.modeffects <- setEffect(fit4.modeffects, Rate, initialValue = 108, name = "support_net", type = "rate", fix = TRUE, verbose = FALSE) ## Explicitly fix the rate parameter lambda at the maximum observed out-degree (in the intersection/source-target verified network)

fit8.ans <- siena07RunToConvergence(alg = modelparams, dat = multidata, eff = fit8.modeffects, ans0 = NULL, modelName = "fit8.ans", batch = TRUE, verbose = FALSE, silent = FALSE, nbrNodes = cores, useCluster = FALSE)



################################# Save All Converged SAOMs for GitHub ################################# 
## As of March 2022, GitHub limits the size of individuals files to 100MB or less (https://docs.github.com/en/repositories/working-with-files/managing-large-files/about-large-files-on-github)
## Accordingly, I cannot upload the entire R workspace post model fitting and goodness-of-fit which is over 1GB
## Still, to save time for replicators, and for posterity, files for all eight fitted models reported in the paper are included on GitHub
## Note that siena07RunToConvergence() will automatically save each fitted SAOM until the model converges. 
## Accordingly, the files saved here are equivalent to the converged model returned by siena07RunToConvergence()
## For example, "SAOM_Model_1_Conventional_Model_Lambda_36.RData" is the same as "fit1.ans_SIENA_Iteration_Number_1_CONVERGED.RData"

save(fit1.ans, file = "SAOM_Model_1_Conventional_Model_Lambda_36.RData")
save(fit2.ans, file = "SAOM_Model_2_Extended_Model_Lambda_36.RData")
save(fit3.ans, file = "SAOM_Model_3_Network_Aid_Model_Restricted_Lambda_36.RData")
save(fit4.ans, file = "SAOM_Model_4_Network_Aid_Model_Full_Lambda_36.RData")

save(fit5.ans, file = "SAOM_Model_5_Conventional_Model_Lambda_108.RData")
save(fit6.ans, file = "SAOM_Model_6_Extended_Model_Lambda_108.RData")
save(fit7.ans, file = "SAOM_Model_7_Network_Aid_Model_Restricted_Lambda_108.RData")
save(fit8.ans, file = "SAOM_Model_8_Network_Aid_Model_Full_Lambda_108.RData")



################################# Combine All SAOM Fit Objects for Post-Processing ################################# 
ans <- list(fit1.ans, fit2.ans, fit3.ans, fit4.ans, ## Lambda = 36 (Main Models)
            fit5.ans, fit6.ans, fit7.ans, fit8.ans) ## Lambda = 108 (Robustness Check Models)

names(ans) <- c("Model_1", "Model_2", "Model_3", "Model_4",
                "Model_5", "Model_6", "Model_7", "Model_8")

intersection.Nicaragua.sienaFits <- ans


closeAllConnections()



################################# TABLE 2 (PART 1), TABLE 3 (PART 1), AND TABLE 4: PARAMETER ESTIMATES AND THE RELATIVE IMPORTANCE OF EFFECTS #################################
all.pretty.effects.of.interest <- c( ## Arranged based on the RSiena internal ordering of effects as they appear in Model 4 and Model 8; RUN: print(fit4.ans)
  "Out-degree",
  "Reciprocity",
  "Transitive Triplets",
  "Transitive Reciprocated Triplets",
  "Three Cycles",
  "Dense Triads",
  "Shared Popularity",
  "In-degree Popularity",
  "Out-degree Popularity",
  "Out-degree Activity",
  "Relative Wealth Rank",
  "Geographic Distance",
  "Consanguineal Relatedness",
  "Consanguineal Relatedness x Reciprocity",
  "Consanguineal Relatedness x Relative Wealth Rank",
  "Consanguineal Relatedness x Geographic Distance",
  "Affinal Relatedness",
  "Godparental Relation",
  "Infidelity Relation",
  "HH Wealth (Alter)",
  "HH Wealth (Ego)",
  "HH Size (Alter)",
  "HH Size (Ego)",
  "Age (Alter)",
  "Age (Ego)",
  "Age Similarity",
  "Gender: Female (Alter)",
  "Gender: Female (Ego)",
  "Same Gender",
  "Ethnicity: Miskito (Alter)",
  "Ethnicity: Miskito (Ego)",
  "Same Ethnicity",
  "Melanin Index (Alter)",
  "Melanin Index (Ego)",
  "Melanin Index Similarity",
  "Body Mass Index (Alter)",
  "Body Mass Index (Ego)",
  "Body Mass Index Similarity",
  "Transitive Triplets Jumping HHs"
)


reorder.all.pretty.effects.of.interest <- c( ## Arranged in the preferred order for tabular presentation in the paper
  "Out-degree",
  
  "Reciprocity",
  "Relative Wealth Rank",
  "Geographic Distance",
  "Consanguineal Relatedness",
  "Consanguineal Relatedness x Reciprocity",
  "Consanguineal Relatedness x Relative Wealth Rank",
  "Consanguineal Relatedness x Geographic Distance",
  "Affinal Relatedness",
  "Godparental Relation",
  "Infidelity Relation",

  "HH Wealth (Alter)",
  "HH Size (Alter)",
  "Age (Alter)",
  "Gender: Female (Alter)",
  "Ethnicity: Miskito (Alter)",
  "Melanin Index (Alter)",
  "Body Mass Index (Alter)",
  
  "HH Wealth (Ego)",
  "HH Size (Ego)",
  "Age (Ego)",
  "Gender: Female (Ego)",
  "Ethnicity: Miskito (Ego)",
  "Melanin Index (Ego)",
  "Body Mass Index (Ego)",
  
  "Age Similarity",
  "Same Gender",
  "Same Ethnicity",
  "Melanin Index Similarity",
  "Body Mass Index Similarity",

  "Out-degree Activity",
  "In-degree Popularity",
  "Out-degree Popularity",
  "Transitive Triplets",
  "Transitive Reciprocated Triplets",
  "Three Cycles",
  "Dense Triads",
  "Transitive Triplets Jumping HHs",
  "Shared Popularity"
)


## Create a list of data frames that contain the results from each fitted SIENA model object
siena.coefs <- lapply(X = rev(intersection.Nicaragua.sienaFits), ## Reverse the order of the fitted SIENA model objects in intersection.Nicaragua.sienaFits to left join with the output from full Model 8
                      ## Remove the first entry in the vector of effect names/estimates theta/std. errors/p-values as they all relate to the rate parameter which is fixed prior to estimation
                      FUN = function(x){ cbind.data.frame(effect = x$effects$effectName[-1], ## Name of each effect
                                                          beta_hat = x$theta[-1], ## Parameter estimates
                                                          se_beta = x$se[-1],  ## Standard error of each parameter estimate
                                                          p_value = 2*pnorm( abs(  x$theta[-1]/x$se[-1] ), lower.tail = FALSE), ## Two-sided p-value associated with each parameter estimate
                                                          
                                                          ####################### CALCULATE THE GLOBAL RELATIVE IMPORTANCE I_k(x) OF EACH EFFECT AT EACH OBSERVATION "WAVE" #######################
                                                          #### For this analysis, we will of course ignore shares of influence for the "second wave" as it is the first wave but altered with the
                                                          #### addition of the random tie as discussed above. Practically speaking, this tie will not make much difference to the results of sienaRI().
                                                          
                                                          #### Measure relative importance on the fly and then retrieve the global/average expected RI of each effect in the fitted model "x"
                                                          RI = sienaRI(multidata, x)$expectedRI[[1]] ## No need for [-1] as  I_k(x) is only returned for effects in the evaluation function, not the rate function
                                                          
                                                          , 
                                                          stringsAsFactors = FALSE) 
                      }
)

## reduce() comes from the library("purr"); left_join() comes from the library("dplyr")
siena.coefs <- reduce(.x = siena.coefs, .f = left_join, by = "effect") ## Left join; https://stackoverflow.com/questions/8091303/simultaneously-merge-multiple-data-frames-in-a-list


## Round all results to the thousandths place for tabular presentation; The first column contains the effect names (character class), hence [,-1]
siena.coefs[,-1] <- apply(siena.coefs[,-1], MARGIN = 2, FUN = function(x){sprintf("%.3f", x)})


## Basic column names
colnames(siena.coefs) <- c("effect",
                           "beta_hat_M8", "se_beta_M8", "p_value_M8", "RI_M8",
                           "beta_hat_M7", "se_beta_M7", "p_value_M7", "RI_M7",
                           "beta_hat_M6", "se_beta_M6", "p_value_M6", "RI_M6",
                           "beta_hat_M5", "se_beta_M5", "p_value_M5", "RI_M5",
                           "beta_hat_M4", "se_beta_M4", "p_value_M4", "RI_M4",
                           "beta_hat_M3", "se_beta_M3", "p_value_M3", "RI_M3",
                           "beta_hat_M2", "se_beta_M2", "p_value_M2", "RI_M2",
                           "beta_hat_M1", "se_beta_M1", "p_value_M1", "RI_M1"
)


## Reorder the columns of the data frame siena.coefs for presentation
siena.coefs <- siena.coefs[c("effect", 
                             "beta_hat_M1", "se_beta_M1", "p_value_M1", "RI_M1",
                             "beta_hat_M2", "se_beta_M2", "p_value_M2", "RI_M2",
                             "beta_hat_M3", "se_beta_M3", "p_value_M3", "RI_M3",
                             "beta_hat_M4", "se_beta_M4", "p_value_M4", "RI_M4",
                             "beta_hat_M5", "se_beta_M5", "p_value_M5", "RI_M5",
                             "beta_hat_M6", "se_beta_M6", "p_value_M6", "RI_M6",
                             "beta_hat_M7", "se_beta_M7", "p_value_M7", "RI_M7",
                             "beta_hat_M8", "se_beta_M8", "p_value_M8", "RI_M8"
)]


rownames(siena.coefs) <- all.pretty.effects.of.interest ## Make the effect names in the first column the official row names
siena.coefs <- siena.coefs[reorder.all.pretty.effects.of.interest, ] ## Indexing by row names, reorder the rows of siena.coefs for presentation
siena.coefs$effect <- NULL ## Remove the first column
siena.coefs[siena.coefs == "NA"] <- "" ## The cells associated with results for effects only in the second/sixth, third/seventh, and forth/eigth model specifications are NA for the other model specifications. Replace with nothing for pretty tabular presentation. 


print(siena.coefs) ## See how it all looks.


#### Use Microsoft Word's convert text to table option (tab delimited)
write.table(siena.coefs[, c("beta_hat_M1", "se_beta_M1", "p_value_M1",
                            "beta_hat_M2", "se_beta_M2", "p_value_M2",
                            "beta_hat_M3", "se_beta_M3", "p_value_M3",
                            "beta_hat_M4", "se_beta_M4", "p_value_M4")],
            file = "T2_PT1_ModelEstimates.txt", sep = "\t", quote = FALSE, row.names = TRUE) ## Main Models

write.table(siena.coefs[, c("beta_hat_M5", "se_beta_M5", "p_value_M5",
                            "beta_hat_M6", "se_beta_M6", "p_value_M6",
                            "beta_hat_M7", "se_beta_M7", "p_value_M7",
                            "beta_hat_M8", "se_beta_M8", "p_value_M8")],
            file = "T3_PT1_ModelEstimates.txt", sep = "\t", quote = FALSE, row.names = TRUE) ## Robustness Check

write.table(siena.coefs[, c("RI_M1", "RI_M5",
                            "RI_M2", "RI_M6",
                            "RI_M3", "RI_M7",
                            "RI_M4", "RI_M8")],
            file = "T4_RI_Effects.txt", sep = "\t", quote = FALSE, row.names = TRUE) ## Relative Importance of Effects in Main Models and Robustness Check Models




################################# Qualitatively Compare Results from the Models Using Lambda = 36 (i.e., Models 1, 2, 3, and 4) and the Models Using Lambda = 108 (i.e., Models 5, 6, 7, and 8) ################################# 
## Conventional Model
siena.coefs[, c("beta_hat_M1", "beta_hat_M5", 
                "se_beta_M1", "se_beta_M5",
                "p_value_M1", "p_value_M5",
                "RI_M1", "RI_M5")] 

## Extended Model
siena.coefs[, c("beta_hat_M2", "beta_hat_M6",
                "se_beta_M2", "se_beta_M6",
                "p_value_M2", "p_value_M6",
                "RI_M2", "RI_M6")]

## Networked Aid Model (Limited)
siena.coefs[, c("beta_hat_M3", "beta_hat_M7", 
                "se_beta_M3", "se_beta_M7",
                "p_value_M3", "p_value_M7",
                "RI_M3", "RI_M7")] 

## Networked Aid Model (Comprehensive)
siena.coefs[, c("beta_hat_M4", "beta_hat_M8", 
                "se_beta_M4", "se_beta_M8",
                "p_value_M4", "p_value_M8",
                "RI_M4", "RI_M8")]

## RI For All Models
print(siena.coefs[, c("RI_M1", "RI_M5", 
                      "RI_M2", "RI_M6",
                      "RI_M3", "RI_M7",
                      "RI_M4", "RI_M8")])



################################# TABLE 2 (PART 2) AND TABLE 3 (PART 2): SAOM GOODNESS OF FIT ################################# 
## RUN: ?sienaGOF
GeodesicDistribution <- function (i, data, sims, period, groupName,
                                  varName, levls = c(1:5, Inf), cumulative = FALSE, ...) {
  x <- networkExtraction(i, data, sims, period, groupName, varName)
  require(network)
  require(sna)
  # a <- geodist(symmetrize(x))$gdist ## http://faculty.ucr.edu/~hanneman/nettext/C7_Connection.html#geodesic
  a <- geodist(x)$gdist ## These are the geodesic distances for directed paths
  if (cumulative)
  {
    gdi <- sapply(levls, function(i){ sum(a <= i) })
  }
  else
  {
    gdi <- sapply(levls, function(i){ sum(a == i) })
  }
  names(gdi) <- as.character(levls)
  gdi
}


CliqueCensus <- function (i, obsData, sims, period, groupName, varName, levls = 1:5){
  require(sna)
  x <- networkExtraction(i, obsData, sims, period, groupName, varName)
  cc0 <- clique.census(x, mode = "graph", tabulate.by.vertex = FALSE,
                            enumerate = FALSE)[[1]]
  cc <- 0*levls
  names(cc) <- as.character(levls)
  levels.used <- as.numeric(intersect(names(cc0), names(cc)))
  cc[levels.used] <- cc0[levels.used]
  cc
}


maxInDegree <- max(degree(tangible.support.matrix.intersection, gmode = "digraph", cmode = "indegree"))
maxOutDegree <- max(degree(tangible.support.matrix.intersection, gmode = "digraph", cmode = "outdegree"))
maxGeodist <- max(geodist(tangible.support.matrix.intersection, inf.replace = -99)$gdist) ## Replace infinite geodesics with -99 to easily retrieve max
maxClique <- max(as.numeric(names(clique.census(tangible.support.matrix.intersection,  mode = "graph", tabulate.by.vertex = FALSE, enumerate = FALSE)$clique.count)))

intersection.Nicaragua.sienaGOFs.indegree <- list()
intersection.Nicaragua.sienaGOFs.outdegree <- list()
intersection.Nicaragua.sienaGOFs.geodist <- list()
intersection.Nicaragua.sienaGOFs.triadcensus <- list()
intersection.Nicaragua.sienaGOFs.cliquecensus <- list()
intersection.Nicaragua.sienaGOFs.consanguineous.ties <- list()

for(model in c("Model_1", "Model_2", "Model_3", "Model_4",
               "Model_5", "Model_6", "Model_7", "Model_8") ) {
  
  print(model)
  
  
  focal.fit <- intersection.Nicaragua.sienaFits[[model]]
  
  
  temp.gof.indegree <- sienaGOF(focal.fit, IndegreeDistribution,
                                varName = "support_net", cumulative = FALSE, levls = 0:maxInDegree, verbose = TRUE)
  print(temp.gof.indegree)
  
  
  temp.gof.outdegree <- sienaGOF(focal.fit, OutdegreeDistribution,
                                 varName = "support_net", cumulative = FALSE, levls = 0:maxOutDegree, verbose = TRUE)
  print(temp.gof.outdegree)
  
  
  temp.gof.geodist <- sienaGOF(focal.fit, GeodesicDistribution,
                               varName = "support_net", cumulative = FALSE, levls = c(1:maxGeodist, Inf), verbose = TRUE)
  print(temp.gof.geodist)
  
  
  temp.gof.triadcensus <- sienaGOF(focal.fit, TriadCensus,
                                   varName = "support_net", levls = 1:16, verbose = TRUE)
  print(temp.gof.triadcensus)
  
  
  temp.gof.cliquecensus <- sienaGOF(focal.fit, CliqueCensus,
                                    varName = "support_net", levls = 1:maxClique, verbose = TRUE)
  print(temp.gof.cliquecensus)
  
  
  temp.gof.consanguineous.ties <- sienaGOF(focal.fit, dyadicCov,
                                           varName = "support_net", dc = degree.of.relatedness, verbose = TRUE)
  print(temp.gof.consanguineous.ties)
  
  
  intersection.Nicaragua.sienaGOFs.indegree[[model]] <- temp.gof.indegree 
  intersection.Nicaragua.sienaGOFs.outdegree[[model]] <- temp.gof.outdegree
  intersection.Nicaragua.sienaGOFs.geodist[[model]] <- temp.gof.geodist
  intersection.Nicaragua.sienaGOFs.triadcensus[[model]] <- temp.gof.triadcensus
  intersection.Nicaragua.sienaGOFs.cliquecensus[[model]] <- temp.gof.cliquecensus
  intersection.Nicaragua.sienaGOFs.consanguineous.ties[[model]] <- temp.gof.consanguineous.ties
  
  
}


rm(model, focal.fit, 
   GeodesicDistribution, CliqueCensus,
   temp.gof.indegree, temp.gof.outdegree, temp.gof.geodist, temp.gof.triadcensus, temp.gof.cliquecensus, temp.gof.consanguineous.ties,
   maxInDegree, maxOutDegree, maxGeodist, maxClique)  



## Plot the distributions to visually compare fit (Main Models).
plot(intersection.Nicaragua.sienaGOFs.indegree[[1]], center = FALSE, scale = FALSE, violin = FALSE)
plot(intersection.Nicaragua.sienaGOFs.indegree[[2]], center = FALSE, scale = FALSE, violin = FALSE)
plot(intersection.Nicaragua.sienaGOFs.indegree[[3]], center = FALSE, scale = FALSE, violin = FALSE)
plot(intersection.Nicaragua.sienaGOFs.indegree[[4]], center = FALSE, scale = FALSE, violin = FALSE)

plot(intersection.Nicaragua.sienaGOFs.outdegree[[1]], center = FALSE, scale = FALSE, violin = FALSE)
plot(intersection.Nicaragua.sienaGOFs.outdegree[[2]], center = FALSE, scale = FALSE, violin = FALSE)
plot(intersection.Nicaragua.sienaGOFs.outdegree[[3]], center = FALSE, scale = FALSE, violin = FALSE)
plot(intersection.Nicaragua.sienaGOFs.outdegree[[4]], center = FALSE, scale = FALSE, violin = FALSE)

plot(intersection.Nicaragua.sienaGOFs.geodist[[1]], center = FALSE, scale = FALSE, violin = FALSE)
plot(intersection.Nicaragua.sienaGOFs.geodist[[2]], center = FALSE, scale = FALSE, violin = FALSE)
plot(intersection.Nicaragua.sienaGOFs.geodist[[3]], center = FALSE, scale = FALSE, violin = FALSE)
plot(intersection.Nicaragua.sienaGOFs.geodist[[4]], center = FALSE, scale = FALSE, violin = FALSE)

plot(intersection.Nicaragua.sienaGOFs.triadcensus[[1]], center = TRUE, scale = TRUE, violin = FALSE)
plot(intersection.Nicaragua.sienaGOFs.triadcensus[[2]], center = TRUE, scale = TRUE, violin = FALSE)
plot(intersection.Nicaragua.sienaGOFs.triadcensus[[3]], center = TRUE, scale = TRUE, violin = FALSE)
plot(intersection.Nicaragua.sienaGOFs.triadcensus[[4]], center = TRUE, scale = TRUE, violin = FALSE)

plot(intersection.Nicaragua.sienaGOFs.cliquecensus[[1]], center = FALSE, scale = FALSE, violin = FALSE)
plot(intersection.Nicaragua.sienaGOFs.cliquecensus[[2]], center = FALSE, scale = FALSE, violin = FALSE)
plot(intersection.Nicaragua.sienaGOFs.cliquecensus[[3]], center = FALSE, scale = FALSE, violin = FALSE)
plot(intersection.Nicaragua.sienaGOFs.cliquecensus[[4]], center = FALSE, scale = FALSE, violin = FALSE)

plot(intersection.Nicaragua.sienaGOFs.consanguineous.ties[[1]], center = TRUE, scale = TRUE, violin = FALSE)
plot(intersection.Nicaragua.sienaGOFs.consanguineous.ties[[2]], center = TRUE, scale = TRUE, violin = FALSE)
plot(intersection.Nicaragua.sienaGOFs.consanguineous.ties[[3]], center = TRUE, scale = TRUE, violin = FALSE)
plot(intersection.Nicaragua.sienaGOFs.consanguineous.ties[[4]], center = TRUE, scale = TRUE, violin = FALSE)



## Plot the distributions to visually compare fit (Robustness Check Models).
plot(intersection.Nicaragua.sienaGOFs.indegree[[5]], center = FALSE, scale = FALSE, violin = FALSE)
plot(intersection.Nicaragua.sienaGOFs.indegree[[6]], center = FALSE, scale = FALSE, violin = FALSE)
plot(intersection.Nicaragua.sienaGOFs.indegree[[7]], center = FALSE, scale = FALSE, violin = FALSE)
plot(intersection.Nicaragua.sienaGOFs.indegree[[8]], center = FALSE, scale = FALSE, violin = FALSE)

plot(intersection.Nicaragua.sienaGOFs.outdegree[[5]], center = FALSE, scale = FALSE, violin = FALSE)
plot(intersection.Nicaragua.sienaGOFs.outdegree[[6]], center = FALSE, scale = FALSE, violin = FALSE)
plot(intersection.Nicaragua.sienaGOFs.outdegree[[7]], center = FALSE, scale = FALSE, violin = FALSE)
plot(intersection.Nicaragua.sienaGOFs.outdegree[[8]], center = FALSE, scale = FALSE, violin = FALSE)

plot(intersection.Nicaragua.sienaGOFs.geodist[[5]], center = FALSE, scale = FALSE, violin = FALSE)
plot(intersection.Nicaragua.sienaGOFs.geodist[[6]], center = FALSE, scale = FALSE, violin = FALSE)
plot(intersection.Nicaragua.sienaGOFs.geodist[[7]], center = FALSE, scale = FALSE, violin = FALSE)
plot(intersection.Nicaragua.sienaGOFs.geodist[[8]], center = FALSE, scale = FALSE, violin = FALSE)

plot(intersection.Nicaragua.sienaGOFs.triadcensus[[5]], center = TRUE, scale = TRUE, violin = FALSE)
plot(intersection.Nicaragua.sienaGOFs.triadcensus[[6]], center = TRUE, scale = TRUE, violin = FALSE)
plot(intersection.Nicaragua.sienaGOFs.triadcensus[[7]], center = TRUE, scale = TRUE, violin = FALSE)
plot(intersection.Nicaragua.sienaGOFs.triadcensus[[8]], center = TRUE, scale = TRUE, violin = FALSE)

plot(intersection.Nicaragua.sienaGOFs.cliquecensus[[5]], center = FALSE, scale = FALSE, violin = FALSE)
plot(intersection.Nicaragua.sienaGOFs.cliquecensus[[6]], center = FALSE, scale = FALSE, violin = FALSE)
plot(intersection.Nicaragua.sienaGOFs.cliquecensus[[7]], center = FALSE, scale = FALSE, violin = FALSE)
plot(intersection.Nicaragua.sienaGOFs.cliquecensus[[8]], center = FALSE, scale = FALSE, violin = FALSE)

plot(intersection.Nicaragua.sienaGOFs.consanguineous.ties[[5]], center = TRUE, scale = TRUE, violin = FALSE)
plot(intersection.Nicaragua.sienaGOFs.consanguineous.ties[[6]], center = TRUE, scale = TRUE, violin = FALSE)
plot(intersection.Nicaragua.sienaGOFs.consanguineous.ties[[7]], center = TRUE, scale = TRUE, violin = FALSE)
plot(intersection.Nicaragua.sienaGOFs.consanguineous.ties[[8]], center = TRUE, scale = TRUE, violin = FALSE)




## Extract information from the sienaGOF objects and combine in one data frame
siena.GOFs <- rbind.data.frame(
  do.call(cbind, lapply(X = intersection.Nicaragua.sienaGOFs.indegree,
                        FUN = function(x){ cbind.data.frame(MHD = x$Joint$ObservedTestStat, p_value = x$Joint$p, stringsAsFactors = FALSE)  }
  ) ),
  do.call(cbind, lapply(X = intersection.Nicaragua.sienaGOFs.outdegree,
                        FUN = function(x){ cbind.data.frame(MHD = x$Joint$ObservedTestStat, p_value = x$Joint$p, stringsAsFactors = FALSE)  }
  ) ),
  do.call(cbind, lapply(X = intersection.Nicaragua.sienaGOFs.geodist,
                        FUN = function(x){ cbind.data.frame(MHD = x$Joint$ObservedTestStat, p_value = x$Joint$p, stringsAsFactors = FALSE)  }
  ) ),
  do.call(cbind, lapply(X = intersection.Nicaragua.sienaGOFs.triadcensus,
                        FUN = function(x){ cbind.data.frame(MHD = x$Joint$ObservedTestStat, p_value = x$Joint$p, stringsAsFactors = FALSE)  }
  ) ),
  do.call(cbind, lapply(X = intersection.Nicaragua.sienaGOFs.cliquecensus,
                        FUN = function(x){ cbind.data.frame(MHD = x$Joint$ObservedTestStat, p_value = x$Joint$p, stringsAsFactors = FALSE)  }
  ) ),
  do.call(cbind, lapply(X = intersection.Nicaragua.sienaGOFs.consanguineous.ties,
                        FUN = function(x){ cbind.data.frame(MHD = x$Joint$ObservedTestStat, p_value = x$Joint$p, stringsAsFactors = FALSE)  }
  ) )
  , stringsAsFactors = FALSE)


## Round all results to the thousandths place for tabular presentation
siena.GOFs <- apply(siena.GOFs , MARGIN = 2, FUN = function(x){sprintf("%.3f", x)})


## Basic column names
rownames(siena.GOFs) <- c("In-degree Distribution", "Out-degree Distribution", "Distribution of Geodesic Distances", "Triad Census", "Clique Census", "Consanguineous Ties")


print(siena.GOFs[, 1:8])
print(siena.GOFs[, 9:16])


#### Use Microsoft Word's convert text to table option (tab delimited)
write.table(siena.GOFs[,c("Model_1.MHD", "Model_1.p_value", 
                          "Model_2.MHD", "Model_2.p_value", 
                          "Model_3.MHD", "Model_3.p_value",
                          "Model_4.MHD", "Model_4.p_value")],
            file = "T2_PT2_ModelGOFs.txt", sep = "\t", quote = FALSE, row.names = TRUE) ## Main Models

write.table(siena.GOFs[,c("Model_5.MHD", "Model_5.p_value", 
                          "Model_6.MHD", "Model_6.p_value", 
                          "Model_7.MHD", "Model_7.p_value",
                          "Model_8.MHD", "Model_8.p_value")],
            file = "T3_PT2_ModelGOFs.txt", sep = "\t", quote = FALSE, row.names = TRUE) ## Robustness Check





################################# DESCRIPTIVE STATISTICS FOR TABLE 1 AND THE MANUSCRIPT TEXT ################################# 
## Monadic Covariates
table(attributes$sex) ## 1 == Female; Male == 0
table(attributes$ethnicity) ## 1 == Miskito; Mayangna == 0


stargazer( data.frame( HH.wealth.log = log( attributes$wealth_hh ), HH.wealth = attributes$wealth_hh, HH.size = attributes$HH.size, age = attributes$age, melanin_index = attributes$melanin_index,  BMI = attributes$bmi), 
           summary = T, summary.logical = T, digits = 2,
           summary.stat = c("n", "mean", "sd", "median", "min", "max"),
           type = "text")


## Dyadic Covariates 
stat.desc(geographic.distance[upper.tri(geographic.distance, diag = F)]) ## Geographic (Inter-household) Distance; symmetric dyadic covariate

stat.desc(degree.of.relatedness[upper.tri(degree.of.relatedness, diag = F)]) ## Consanguineal Relatedness; symmetric dyadic covariate
stat.desc(affinal.degree.of.relatedness[upper.tri(affinal.degree.of.relatedness, diag = F)]) ## Affinal Relatedness; symmetric dyadic covariate

table(fictive.relatedness[upper.tri(fictive.relatedness, diag = F)]) ## Godparental Relation; symmetric dyadic covariate
table(infidelity.relatedness[upper.tri(infidelity.relatedness, diag = F)]) ## Infidelity Relation; symmetric dyadic covariate

range(relative.wealth_hh.rank) ## asymmetric dyadic covariate


## How many unique dyads are comprised of close kin? (Presented in the Discussion of the paper)
degree.of.relatedness.at.least.0.125 <- degree.of.relatedness
degree.of.relatedness.at.least.0.125[degree.of.relatedness.at.least.0.125 < 0.125] <- 0 
degree.of.relatedness.at.least.0.125[degree.of.relatedness.at.least.0.125 >= 0.125] <- 1 

affinal.degree.of.relatedness.0.125 <- affinal.degree.of.relatedness
affinal.degree.of.relatedness.0.125[affinal.degree.of.relatedness.0.125 < 0.125] <- 0 
affinal.degree.of.relatedness.0.125[affinal.degree.of.relatedness.0.125 >= 0.125] <- 1 

close.kin <- degree.of.relatedness.at.least.0.125 + affinal.degree.of.relatedness.0.125
table(close.kin[upper.tri(close.kin, diag = F)]) ## How many close kin dyads?
table(close.kin[upper.tri(close.kin, diag = F)])/sum(table(close.kin[upper.tri(close.kin, diag = F)])) ## Proportion of dyads that are close kin?


## How many asymmetric ties occur between close kin?
table(tangible.support.matrix.intersection*close.kin) 


## How many dyads and triads amongst close kin? (Referenced in the Discussion of the paper)
dyad.census(tangible.support.matrix.intersection)
dyad.census(tangible.support.matrix.intersection*close.kin)

triad.census(tangible.support.matrix.intersection)
triad.census(tangible.support.matrix.intersection*close.kin)


## Snijders' Degree of Certainty
## Snijders, Tom A.B. ‘Explained Variation in Dynamic Network Models’. Mathématiques et Sciences Humaines, no. 168 (1 December 2004). https://doi.org/10.4000/msh.2938.

sienaRI(multidata, fit1.ans)
sienaRI(multidata, fit5.ans)

sienaRI(multidata, fit2.ans)
sienaRI(multidata, fit6.ans)

sienaRI(multidata, fit3.ans)
sienaRI(multidata, fit7.ans)

sienaRI(multidata, fit4.ans)
sienaRI(multidata, fit8.ans)





