args <- commandArgs()

#define the range of data to be analysed 
from <- args[6]
to <- args[7]
print(args)

library(purrr)
library(readr)
library(dplyr)
library(nlme)
library(inline)

#read in the accession info
allGenesAllAccnoNAFilt <- read_csv("data/allGenesAllAccnoNAFilt.csv") %>%
	mutate(GROUP = if_else(GROUP %in% c("admixed",
								  "south_sweden",
								  "central_europe",
								  "western_europe",
								  "spain",
								  "italy_balkan_caucasus"),
					 "merged",
					 GROUP))

#read in the phenotype experiments
pheno <- read_csv("data/MergedPhenotypes.csv",
 col_names = FALSE) %>% 
  arrange(X1) %>% 
  rename(EXPERIMENT = X1,
         ECOTYPE = X2,
         NAME = X3,
         CS_CODE = X4,
         LONGITUDE = X5,
         LATITUDE = X6,
         COUNTRY = X7,
         PHENO_SCORE = X8)

#read list of experiments with a non-binary PHENO_SCORE distribution
experiments <- read_csv("data/pheno_experiments_non_binary.csv")

#define the functions needed for the analyses
linMixFilterPhenoControl <- function(x)
  x %>% 
  select(SUM_VAR, PHENO_SCORE, GENE, CHR_START, UNIT, GROUP, WHERE, LEN_REF) %>% 
  group_by(GROUP) %>%
  filter(n() > 2) %>%
  ungroup()

linMixEffNlmePheno <- function(startpos){
  u <- lme(PHENO_SCORE ~ SUM_VAR, random = ~1|GROUP, data = startpos) #make the model
  v <- summary(u) #create object of the summary
  
  return(v)
}

lmeDataFramePheno <- function(resuList){
  modelResu <- data_frame()
  corcoeff  <- data_frame()
  randeffVar <- data_frame()
  randeffStdDev <- data_frame()
  geneName <- data_frame()
  startLocus <- data_frame()
  unitLen <- data_frame()
  geneLoc <- data_frame()
  lenRef <- data_frame()
  collect_row  <- data_frame()
  for (i in seq(length(resuList))) {
    if (!is.null(resuList[[i]])) {
      #lapply(resuList[[i]], with, print(resuList[[i]][[j]]$coefficients))
      
      modelResu <- resuList[[i]]$tTable %>% as_data_frame() %>% .[2,]
      
      corcoeff <- resuList[[i]]$corFixed  %>% .[2,1]%>% as_data_frame()
      
      randeffVar <- nlme::VarCorr(resuList[[i]]) %>% .[1,1] %>% as_data_frame()
      
      randeffStdDev <- nlme::VarCorr(resuList[[i]])  %>% .[1,2] %>% as_data_frame()
      
      
    
    
    geneName <- resuList[[i]]$data$GENE %>% unique() %>% as.character()%>% as.data.frame(stringsAsFactors=FALSE)
    startLocus <- resuList[[i]]$data$CHR_START %>% unique() %>% as.character() %>% as.data.frame(stringsAsFactors = FALSE)
    unitLen <- resuList[[i]]$data$UNIT %>% unique() %>% as.character()%>% as.data.frame(stringsAsFactors=FALSE)
    geneLoc <- resuList[[i]]$data$WHERE %>% unique() %>% as.character() %>% as_data_frame()
    lenRef <- resuList[[i]]$data$LEN_REF %>% unique() %>% as.character() %>% as_data_frame()
    collect_row <- bind_cols(modelResu, corcoeff, geneName, startLocus, geneLoc, unitLen, lenRef,randeffVar, randeffStdDev) 
    
    names(collect_row) <- c("Estimate", "Std. error", "DF", "t_value", "p_value", "Correlation", "Gene", "Start",
                  "Where", "Unit", "LenRef", "Random_effect_Var", "Random_effect_stdDev")
    lmeResDF <<- bind_rows(lmeResDF, collect_row)
    }
    modelResu <- data_frame()
    corcoeff <- data_frame()
    randeffVar <- data_frame()
    randeffStdDev <- data_frame()
    geneName <- data_frame()
    startLocus <- data_frame()
    unitLen <- data_frame()
    geneLoc <- data_frame()
    lenRef <- data_frame()
    collect_row  <- data_frame()
  }
}


#read through the number, as defined by commandArgs, of the entries in experiment list 
for (line in experiments$value[from:to]) {
  #make the dataframe for each experiment
  
  paste(line, "start") %>% print()
  
  temp <- pheno %>% 
    filter(EXPERIMENT ==line) 
  
  #merge with repeat info DF
  temp2 <- inner_join(allGenesAllAccnoNAFilt, temp %>% select(ECOTYPE, EXPERIMENT, PHENO_SCORE), by = "ECOTYPE")
  
  startList <- list()
  
  #make each repeat a DF in a list
  startList <- split(temp2, temp2$CHR_START)
	
  
  #apply the filtering function to each DF in the list
  startList <- lapply(startList, linMixFilterPhenoControl)
  
  testList <- list()
  
  resuList <- list()
  
  too_few_var <- which(map(startList, ~n_distinct(.$SUM_VAR)) < 3)
  
  #apply the lme modelling function to each DF
  startList_filt <- startList[-too_few_var] %>% 
    .[sapply(., nrow) > 0]
  

  
  resuList <- lapply(startList_filt,function(x) tryCatch(linMixEffNlmePheno(x), error=function(e) NULL))
    
  print(!is.null(resuList[[1]]))
  
  #gather the results in a DF
  lmeResDF <- data_frame()

  resuList %>% lmeDataFramePheno()
  
  filename <- paste("data/phenotypes/observed/lmeResDF_", line, ".csv", sep = "")
  #write the results to a separate DF
  write_csv(lmeResDF, filename)
  
  paste(line, "end") %>% print()
    
  gc()
  
}


  
  