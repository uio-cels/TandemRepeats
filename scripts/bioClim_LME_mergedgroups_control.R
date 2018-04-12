args <- commandArgs()

#define which genomic region to analyse
region <- args[6]

library(parallel)
library(purrr)
library(readr)
library(dplyr)

#read in the data frame containing the merged STR and bioclimatic data
allGenesAllAccBioClim <- read_csv(
file = "data/AllGenesAllAccBioClim.csv"
) %>% 
ungroup() %>%
mutate(GROUP = if_else(GROUP %in% c("admixed",
								  "south_sweden",
								  "central_europe",
								  "western_europe",
								  "spain",
								  "italy_balkan_caucasus"),
					 "merged",
					 GROUP))

print("read file: yes")

#make a filtering function to remove the bioclimatic variables which have a correlation 
#coefficient > 0.8 and where each population group contains less than three observations 
#and randomising the STR lengths 
linMixFilter <- function(x)
  x %>% 
  select(starts_with("BIO_"), SUM_VAR, GENE, CHR_START, UNIT, GROUP, WHERE) %>% 
  select(-BIO_TEMP_ANNUAL_RANGE, -BIO_MIN_TEMP_COLDEST_MONTH, -BIO_MEAN_TEMP_COLDEST_QUART, 
         -BIO_PREC_WETTEST_MONTH, -BIO_PREC_WETTEST_QUART, -BIO_PREC_DRIEST_MONTH) %>%
  group_by(GROUP) %>%
  filter(n() > 2) %>%
  ungroup() %>%
  mutate(SUM_VAR = sample(SUM_VAR))


print("define function: yes")

#make a list of all observations for each STR
startList <- split(allGenesAllAccBioClim, allGenesAllAccBioClim$CHR_START)

#find which STRs contain less than three different length variants
too_few_var <- which(map(startList, ~n_distinct(.$SUM_VAR)) < 3)

#allocate clusters
options(mc.cores = 4)

makeForkCluster(nnodes = getOption("mc.cores"))

startList <- mclapply(startList, linMixFilter)

print("apply function:yes")

#make the function for creating the linear mixed-effect model, returning the models for 
#each STR as a list
linMixEffNlme <- function(startpos){
  require(nlme)
  for (i in 1:(ncol(startpos)-6)) {
    lmeformu <- as.formula(paste('SUM_VAR ~ ', paste(colnames(startpos[i]))))
    u <- lme(lmeformu, random = ~1|GROUP, data = startpos) 
    v <- summary(u) #create object of the summary    
    testList[[i]] <- v 
  }
  return(testList)
}


#make the function for collecting the model results as a data frame 
lmeDataFrame <- function(resuList){
  bioClimVar <- list()
  modelResu <- list()
  corcoeff  <- list()
  randeffVar <- list()
  randeffStdDev <- list()
  for (i in seq(length(resuList))) {
    if (!is.null(resuList[[i]])) {
      #lapply(resuList[[i]], with, print(resuList[[i]][[j]]$coefficients))
      for (j in seq(length(resuList[[i]]))) {
        bioClimVar[[j]] <- resuList[[i]][[j]]$terms %>%
          as.list() %>%
          .[[3]] %>%
          as.character()
        
        modelResu[[j]] <- resuList[[i]][[j]]$tTable %>% as.data.frame() %>% .[2,]
        
        corcoeff[[j]] <- resuList[[i]][[j]]$corFixed %>% as.data.frame() %>% .[2,1]
        
        randeffVar[[j]] <- nlme::VarCorr(resuList[[i]][[j]]) %>% .[1,1]
        
        randeffStdDev[[j]] <- nlme::VarCorr(resuList[[i]][[j]]) %>% .[1,2]
        
        
      }
      
      modelResu <- bind_rows(lapply(modelResu, as.data.frame, stringsAsFactors=FALSE))
      bioClimVar <- bind_rows(lapply(bioClimVar, as.data.frame, stringsAsFactors=FALSE))
      corcoeff <- bind_rows(lapply(corcoeff, as.data.frame, stringsAsFactors = FALSE))
      randeffVar <- bind_rows(lapply(randeffVar, as.data.frame, stringAsFactors = FALSE))
      randeffStdDev <- bind_rows(lapply(randeffStdDev, as.data.frame, stringsAsFactors = FALSE))
      geneName <- rep(unique(resuList[[i]][[1]]$data$GENE), 19) %>% as.character()%>% as.data.frame(stringsAsFactors=FALSE)
      startLocus <- rep(unique(resuList[[i]][[1]]$data$CHR_START), 19) %>% as.character() %>% as.data.frame()
      unitLen <- rep(unique(resuList[[i]][[1]]$data$UNIT), 19) %>% as.character()%>% as.data.frame(stringsAsFactors=FALSE)
      e <- bind_cols(bioClimVar, modelResu, corcoeff, geneName, startLocus, unitLen, randeffVar, randeffStdDev) 
      names(e) <- c("VARIABLE", "Estimate", "Std. error", "DF", "t_value", "p_value", "Correlation", "Gene", "Start", "Unit", "Random_effect_Var", "Random_effect_stdDev")
      lmeResDF <<- bind_rows(lmeResDF, e) %>% dplyr::select(8,9,everything())
      
      bioClimVar <- list()
      modelResu <- list()
      corcoeff <- list()
      randeffVar <- list()
      randeffStdDev <- list()
    }
  }
}


testList <- list()




#apply the linear mixed-effect model function on the STRs in the defined region from calling 
#the script from the terminal
resuList <- startList[-too_few_var] %>% 
.[sapply(.,  function(x) all(x$WHERE==region))] %>%
.[sapply(., nrow) > 0] %>% 
mclapply(.,function(x) tryCatch(linMixEffNlme(x), error=function(e) NULL))

print("LME function applied")

lmeResDF <- data_frame()

#gather the model summaries in a data frame
resuList %>% lmeDataFrame()

#write to file
filename <- paste("data/lmeResDF_mergedgroups_control_", region, ".csv", sep = "")
write_csv(lmeResDF, filename)
print("analysed:yes")


