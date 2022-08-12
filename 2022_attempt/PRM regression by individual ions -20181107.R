#####################################################################
## 2018-10-29: (1) Add s.e. when reporting lm models  
##             (2) maxit = 50 in rlm
#####################################################################
require(dplyr)

## Setup of data Directory
setwd("c:\\Results\\AUDIT")
current.Dir <- getwd()
new.Dir <- ".\\Regression_inputs"
dir.create(file.path(current.Dir, new.Dir), showWarnings = FALSE)
dir.create(file.path(current.Dir, "Inputs"), showWarnings = FALSE)

Sample.Types <- c("Neat", "TD", "PCS") # Sample types to differentiate
Curve.Types <- c("Fwd", "Rev") # Sample types to differentiate

## Step 1: conversion of PRM outputs for Regression inputs
sum.intensity <- function(df, Curve) {
  df$fragment_ions <- paste0("+",df$Precursor.Charge,df$Fragment.Ion)
  require(plyr)
  require(reshape2)
  dfw <- ddply(df, .(Peptide, File_Name, Isotope.Label.Type, Replicate, SampleGroup, Theor_conc, fragment_ions), summarise,
               Sum=sum(as.numeric(Area)))
  # dfw$Peptide <- paste(dfw$Peptide, dfw$fragment_ions, sep=".")
  
  dfw <- reshape(dfw, idvar=c("Peptide", "fragment_ions", "SampleGroup", "Replicate"), timevar="Isotope.Label.Type", direction="wide")
  dfw <-  dfw[,!grepl("fragment_ions.Light|Theor_conc.Light|File_Name.Light", names(dfw), ignore.case=TRUE)]
  colnames(dfw) <- gsub("\\.","_",names(dfw))
  colnames(dfw) <- gsub("Theor_conc_heavy","Theor_conc",names(dfw), ignore.case=TRUE)
  colnames(dfw) <- gsub("fragment_ions_heavy","fragment_ions",names(dfw), ignore.case=TRUE)
  colnames(dfw) <- gsub("File_Name_heavy","File_Name",names(dfw), ignore.case=TRUE)
  colnames(dfw) <- gsub("Sum_Heavy","Sum_heavy",names(dfw), ignore.case=TRUE)
  colnames(dfw) <- gsub("Sum_Light","Sum_light",names(dfw), ignore.case=TRUE)

  dfw <- dfw[,!names(dfw)=="SampleGroup"]
  if(Curve=="Rev") {
    dfw$"PAR_(H/L)" <- dfw$Sum_heavy/dfw$Sum_light
    dfw$Meas_conc <- dfw$"PAR_(H/L)"*df$IS.Spike[1]
  } else {
    dfw$"PAR_(L/H)" <- dfw$Sum_light/dfw$Sum_heavy
    dfw$Meas_conc <- dfw$"PAR_(L/H)"*df$IS.Spike[1]
  } 

  ColNames <- names(dfw)
  dfw <- do.call(data.frame,lapply(dfw, function(x) replace(x, is.infinite(x),NA)))  # avoid infinite values
  colnames(dfw) <- ColNames
  return(dfw)
}

filelist = list.files(path = ".\\Inputs", pattern="*\\.csv$")  # Find all csv files
output.prefix <- data.frame(do.call('rbind', strsplit(as.character(filelist),'.csv',fixed=TRUE))) # file name prefix for outputs
for (i in 1:length(filelist)) {
  assign("df", read.csv(file.path(".\\Inputs", filelist[i]), check.names=TRUE, header = T, comment.char="#"))  # Import data files
  for (j in 1:length(Sample.Types)) { if (grepl(Sample.Types[j], filelist[i], ignore.case=TRUE)) Sample <- Sample.Types[j]} # Determine sample type
  for (j in 1:length(Curve.Types)) { if (grepl(Curve.Types[j], filelist[i], ignore.case=TRUE)) Curve <- Curve.Types[j]} # Determine curve type

  # ifelse(grepl("Rev", filelist[i]), IS.Reversed <- TRUE, IS.Reversed <- FALSE)
  df$Peptide <- paste(df$"Protein.Name",df$"Peptide.Sequence",sep=".")
  colnames(df)[names(df)=="File.Name"] <- "File_Name"
  colnames(df)[names(df)=="Concentration"] <- "Theor_conc"
  df$Replicate <- as.factor(df$Replicate)
  df$Theor_conc <- as.factor(df$Theor_conc)
  df <- df[complete.cases(df),]

  dfw <- sum.intensity(df, Curve)
  # dfw <- dfw[dfw$Meas_conc!=0,]
  write.table(dfw, file.path(current.Dir, new.Dir, paste0(output.prefix[i,1], " regression inputs", ".txt")), sep="\t", col.names = TRUE, row.names=FALSE)	
}
## End of conversion



## Step 2: Data fitting by Regression 
setwd("c:\\Results\\AUDIT\\Regression_inputs")

####################################################################################
#  LoB estimation by multiple methods
#  Jingqin (Rosy) Luo            07/26/2013
#####################################################################################

LoB <- function(x)
{
  ##method1:Linett paper
  ##x: a vecotr of values
  x <- na.exclude(sort(x))
  Nb <- length(x)
  obs.idx <- Nb*0.95+0.5
  if(obs.idx%in%(1:Nb)) quan <- x[obs.idx]
  else if (obs.idx>1 & obs.idx<Nb)      
  {
    #quan <- 0.5*(x[floor(obs.idx)]+x[ceiling(obs.idx)])
    quan <- x[floor(obs.idx)]+(obs.idx-floor(obs.idx))*(x[ceiling(obs.idx)]-x[floor(obs.idx)])##linear interpolation
  }
  else if(obs.idx>Nb) quan <- x[Nb]
  else if(obs.idx<1) quan <- x[1]
  else
  {
    stop("sth wrong!")
  }
  return(quan)
}

LoB2 <- function(x)
{
  ##method2: Linett paper, remove the "0.5" in obs.idx
  ##x: a vecotr of values
  
  x <- na.exclude(sort(x))
  Nb <- length(x)
  obs.idx <- Nb*0.95##remove 0.5
  if(obs.idx%in%(1:Nb)) quan <- x[obs.idx]
  else if (obs.idx>1 & obs.idx<Nb)      
  {
    #quan <- 0.5*(x[floor(obs.idx)]+x[ceiling(obs.idx)])
    quan <- x[floor(obs.idx)]+(obs.idx-floor(obs.idx))*(x[ceiling(obs.idx)]-x[floor(obs.idx)])##linear interpolation
  }
  else if(obs.idx>Nb) quan <- x[Nb]
  else if(obs.idx<1) quan <- x[1]
  else {
    browser()
    stop("sth wrong!")
  }
  return(quan)
}

LoB3 <- function(x)
{
  ##method3: complement with negative values
  ##x: a vecotr of values
  
  x <- na.exclude(sort(x))
  Nb <- length(x)
  neg.x <- -1*x
  quan <- as.numeric(quantile(c(neg.x,x),prob=0.95))
  return(quan)
}

LoB4 <- function(x)
{
  ##method4: Armbruster paper, mean+1.645*sd of blank obs
  ##x: a vecotr of values
  quan <- as.numeric(mean(x,na.rm=T)+1.645*sd(x,na.rm=T))
  return(quan)
}

LoB.Est.method <- function(x)
{
  ###wrap for all the four methods
  est.quan1 <- LoB(x)
  est.quan2 <- LoB2(x)
  est.quan3 <- LoB3(x)
  est.quan4 <- LoB4(x)
  return(c(method1=est.quan1,method2=est.quan2,method3=est.quan3,method4=est.quan4))
  
}

file_name = list.files(path = ".", pattern="*regression inputs\\.txt$")  # Find all input text files
summary_data <- data.frame()
individual_data <- NULL

for (current_file in file_name) {
  # current_file <- file_name[1]
  cptac <- read.table(current_file, sep="\t", header=TRUE)
  cptac <- cptac[!is.na(cptac$Peptide),]
  
  ## filter out peptidew with less than two injections
  # temp <- cptac[cptac$Peptide =="G6PI_HUMAN.VWYVSNIDGTHIAK",]
  temp <- ddply(cptac, .(Peptide, Theor_conc), summarise, freq=length(File_Name))
  temp <- temp[temp$freq >= 3,]
  cptac <- cptac[cptac$Peptide %in% temp$Peptide,]
  rm(temp)

  peptide_name <- unique(cptac[ ,'Peptide'])
  ion_name <- cptac[ , 'fragment_ions']
  unique_ion <- unique(ion_name)
  

  for (i in peptide_name) {
    # i <- peptide_name[1] 
    label_ion <- i
    cptac1 <- cptac[cptac$Peptide==i,]  # data of peptide i
	cptac1 <- cptac1[!is.na(cptac1$Meas_conc),]
    Fragment_ions <- unique(cptac1$fragment_ions)
    Fragment_ions <- paste(Fragment_ions, collapse = '')
    cptac1 <- cptac1[order(cptac1$Theor_conc, cptac1$Replicate),]
	cptac1_display <- cptac1[order(cptac1$Theor_conc, cptac1$Replicate), c('Peptide', 'Replicate', 'Theor_conc', 'Meas_conc', 'fragment_ions')]
	
	## LOB calculation
	blanks <- cptac1[with(cptac1, {Theor_conc==0}), 'Meas_conc']
	LOB <- LoB3(blanks)
	if(is.na(LOB)) LOB <- 0
	
	## repalce zero with NA
	temp <- cptac1[with(cptac1, {Theor_conc > 0}), "Meas_conc", drop = FALSE]
	temp[temp$Meas_conc==0, ] <- NA
	cptac1[with(cptac1, {Theor_conc > 0}), "Meas_conc"] <- temp
	rm(temp)
	
	## degree of freedom
	deg_f <- ddply(cptac1[with(cptac1, {Theor_conc > 0}), ], .(Theor_conc, fragment_ions), summarise, freq = length(Theor_conc), SD = sd(Meas_conc, na.rm=T))
	deg_f$f <- deg_f$freq - 1
	deg_f <- deg_f[deg_f$f >= 2, ]
	
	## Minimum concentation for LOD calculation
	temp <- ddply(deg_f, .(Theor_conc), summarise, N_fragments = length(fragment_ions), N_na = sum(is.na(SD))) 
	temp <- temp[with(temp, {N_na <= N_fragments/2}),]
	min_conc <- min(temp$Theor_conc)
	rm(temp)
	
	## LOD calculation
	z_score <- 1.645 # one-sided z-score at 95% percentile
	
	## method (1): LOD for each transition -> mean LOD
	cptac1_lod <- deg_f[with(deg_f, {Theor_conc == min_conc & !is.na(SD)}), ]
	cptac1_lod$LOD <- LOB + z_score *1/(1 - 1/(4*cptac1_lod$f))*cptac1_lod$SD 
	LOD_1 <- mean(cptac1_lod$LOD, na.rm=T)
	
	individual_data <- rbind(individual_data, cbind.data.frame(File_Name=current_file, Peptide=i, cptac1_lod))
	
	## method (2) mean standard deviation --> mean LOD
	mean_sd <- mean(cptac1_lod$SD)
	tot_f <- sum(cptac1_lod$f) # total degree of freedom 
	LOD_2 <- LOB + z_score *1/(1 - 1/(4*tot_f))*mean_sd

	## use method (1)
	LOD <- LOD_1

    # regression all points above LOD
    above_LOD <- cptac1[cptac1$Meas_conc>=LOD, ]
    above_LOD_sorted <- above_LOD[order(above_LOD$Meas_conc), c('Peptide', 'Replicate', 'Theor_conc', 'Meas_conc', 'fragment_ions')]
    
    # no NA included in regression calculation
    above_LOD <- na.omit(above_LOD)
    above_LOD_sorted <- na.omit(above_LOD_sorted)
    above_LOD_sorted_count <- nrow(above_LOD_sorted)
    
	run.scripts <- FALSE
	if (run.scripts) {
		# bottom_6_samples <- cptac1[cptac1$Theor_conc<=0.25,'Meas_conc']
		bottom_6_samples <- cptac1[1:ifelse(nrow(cptac1)>=6,6,nrow(cptac1)),'Meas_conc']
		bottom_6_samples_sorted <- sort(bottom_6_samples, decreasing=TRUE)
		Theor_concentration_list_of_bottom_6 <- cptac1[cptac1$Theor_conc<=0.25,'Theor_conc']
		
		Theor_concentration_list_of_bottom_6 <- unique(Theor_concentration_list_of_bottom_6)
		Theor_concentration_list_of_bottom_6 <- sort(Theor_concentration_list_of_bottom_6)

		Theor_concentration_list_above_bottom_6 <- cptac1[cptac1$Theor_conc>0.25,'Theor_conc']
		Theor_concentration_list_above_bottom_6 <- unique(Theor_concentration_list_above_bottom_6)
		Theor_concentration_list_above_bottom_6 <- sort(Theor_concentration_list_above_bottom_6)

		# calculate LOB using Rosy's methods for bottom 6 Theor_conc sets (Theor_conc<=0.25)
		LOB_results <- LoB.Est.method(bottom_6_samples)
		LOB_method1 <- LoB(bottom_6_samples)
		LOB_method2 <- LoB2(bottom_6_samples)
		LOB_method3 <- LoB3(bottom_6_samples)
		LOB_method4 <- LoB4(bottom_6_samples)
		
		LOB_for_LOD_calc <- LOB_method3
		# bottom_6_sample_count <- length( bottom_6_samples_sorted)

		
		##  LOD calculation
		min_counter <- 0
		LOD_values <- data.frame()
		LOD_Theor <- c()
		
		for (x in Theor_concentration_list_above_bottom_6){  # go through different concentrations
		  current_set <- cptac1[cptac1$Theor_conc==x, c('Theor_conc', 'Meas_conc')]
		  if (nrow(current_set) != 3){ print("*** Error: current Theor_conc set does not contain 3 replicates *** ", quote=FALSE) }
		  is_na <- any(is.na(current_set$Meas_conc))
		  
		  if(is_na){
			print("*** Theor_conc set contains NA ... skipping *** ", quote=FALSE)
			print (current_set$Theor_conc, quote=FALSE)
		  }
		  else{
			if (min(current_set$Meas_conc) > LOB_for_LOD_calc & min_counter <3) {
			  min_counter <- min_counter + 1
			  LOD_values <- rbind(LOD_values, current_set)
			  LOD_Theor <- c(LOD_Theor, x)
			}
		  } #end else
		} # end for loop
		##  End of LOD calculation

		SDs_top <- c() # standard deviations

		for (y in LOD_Theor){
		  current_SDs <- cptac1[cptac1$Theor_conc==y, c('Meas_conc')]
		  # print (current_SDs, quote=FALSE)
		  SDs_top <- c(SDs_top, 2*sd(current_SDs)*sd(current_SDs))
		}
		
		# bottom of SDs2 equation = number of samples = number sample sets X degree of freedom for each sample = 3 X (3 - 1) = 6
		SDs_bottom <- 6
		SDs <- sqrt(sum(SDs_top)/SDs_bottom)
		z1_b <- 1.645
		f <- 6
		middleLOD_term <- (z1_b/(1 - (1/(4*f))))
		# calculate (1 - 1/(4 X f))
		# (1 - 1/(4 X 6))
		# 4 * 6 = 24    0.042
		LOD_old <- LOB_method3 + (1.645/0.958) * SDs
		LOD <- LOB_method3 + middleLOD_term * SDs
		# ******* end LOD calculation *******
		
		# regression all points above LOD
		above_LOD <- cptac1[cptac1$Meas_conc>=LOD & cptac1$Theor_conc>0.25, ]
		above_LOD_sorted <- above_LOD[order(above_LOD$Meas_conc), c('Peptide', 'Replicate', 'Theor_conc', 'Meas_conc', 'fragment_ions')]
		
		# no NA included in regression calculation
		above_LOD <- na.omit(above_LOD)
		above_LOD_sorted <- na.omit(above_LOD_sorted)
		above_LOD_sorted_count <- nrow(above_LOD_sorted)
		
	
	
	}
	

    # ***** Regression *****
    if(length(unique(above_LOD$Theor_conc)) > 2) {
      model = lm(above_LOD$Meas_conc ~ above_LOD$Theor_conc)
      model_coef <- coef(model)
      model_intercept <- as.numeric(model_coef[1])
      model_slope <- as.numeric(model_coef[2])
      model_Rsquared <- summary(model)$r.squared
	  model_Sigma <- summary(model)$sigma
    } else {
      model_intercept <- NA
      model_slope <-NA
      model_Rsquared <- NA
	  model_Sigma <- NA
    }
    
    # fits a model using MM-estimation
    library(MASS)
    if(length(unique(above_LOD$Theor_conc)) > 2) {
      model_robust1 = rlm(above_LOD$Meas_conc ~ above_LOD$Theor_conc, method  = "MM", maxit = 50)
      model_robust1_coef <- coef(model_robust1)
      model_robust1_intercept <- as.numeric(model_robust1_coef[1])
      model_robust1_slope <- as.numeric(model_robust1_coef[2])
	  
	  model_robust1_summary <- summary(model_robust1)
	  model_robust1_Rsquared <- model_robust1_summary$r.squared
	  model_robust1_intercept_SE <- model_robust1_summary$coefficients[, 2][1]
	  model_robust1_slope_SE <- model_robust1_summary$coefficients[, 2][2]
	  model_robust1_Sigma <- model_robust1_summary$sigma
    } else {
      model_robust1_coef <- NA
      model_robust1_intercept <-NA
	  model_robust1_intercept_SE <- NA
      model_robust1_slope <- NA
	  model_robust1_slope_SE <- NA
      model_robust1_Rsquared <- NA
	  model_robust1_Sigma <- NA
    }

    # fits a model using MS- and S-estimation
    library(robust)
    if(length(unique(above_LOD$Theor_conc)) > 2) {
      model_robust2 = lmRob(above_LOD$Meas_conc ~ above_LOD$Theor_conc)
      model_robust2_coef <- coef(model_robust2)
      model_robust2_intercept <- as.numeric(model_robust2_coef[1])
      model_robust2_slope <- as.numeric(model_robust2_coef[2])
      model_robust2_Rsquared <- summary(model_robust2)$r.squared
	  model_robust2_Sigma <- summary(model_robust2)$sigma
    } else {
      model_robust2_coef <- NA
      model_robust2_intercept <-NA
      model_robust2_slope <- NA
      model_robust2_Rsquared <- NA
	  model_robust2_Sigma <- NA
    }
    
    # fits a model using LQS, includes least median squares
    library(robust)
    if(length(unique(above_LOD$Theor_conc)) > 2) {
      model_robust3 = lqs(above_LOD$Meas_conc ~ above_LOD$Theor_conc)
      model_robust3_coef <- coef(model_robust3)
      model_robust3_intercept <- as.numeric(model_robust3_coef[1])
      model_robust3_slope <- as.numeric(model_robust3_coef[2])
    } else {
      model_robust3_coef <- NA
      model_robust3_intercept <-NA
      model_robust3_slope <- NA
    }
    
    
    if(length(unique(above_LOD$Theor_conc)) > 2) {
      model_robust4 = lqs(above_LOD$Meas_conc ~ above_LOD$Theor_conc, method="lms")
      model_robust4_coef <- coef(model_robust4)
      model_robust4_intercept <- as.numeric(model_robust4_coef[1])
      model_robust4_slope <- as.numeric(model_robust4_coef[2])
      model_robust4_coef <- coef(model_robust4)
    } else {
      model_robust4_coef <- NA
      model_robust4_intercept <-NA
      model_robust4_slope <- NA
    }
    # ***** end Regression *****
    
	for (j in 1:length(Sample.Types)) { if (grepl(Sample.Types[j], current_file, ignore.case=TRUE)) Sample <- Sample.Types[j]} # Determine sample type
	for (j in 1:length(Curve.Types)) { if (grepl(Curve.Types[j], current_file, ignore.case=TRUE)) Curve <- Curve.Types[j]} # Determine curve type
	
	Pep_display <- gsub(".*\\|", "", i) %>% gsub("(.*\\.[A-Z]{3}).*", "\\1", .)
	
    summary_data <- rbind(summary_data, data.frame(Peptide = i, Pep_display = Pep_display, Fragment_ions = Fragment_ions, 
                                                   # LoB_method1 = LOB_method1, LoB_method2 = LOB_method2, LoB_method3 = LOB_method3, LoB_method4 = LOB_method4, 
												   LOB = LOB, LOD = LOD,
                                                   lm_intercept = model_intercept, lm_slope = model_slope, lm_Rsquared = model_Rsquared, lm_Sigma = model_Sigma, 
                                                   rlm_intercept = model_robust1_intercept, rlm_intercept_SE = model_robust1_intercept_SE, rlm_slope = model_robust1_slope, rlm_slope_SE = model_robust1_slope_SE, 
												     rlm_Rsquared = model_robust1_Rsquared, rlm_Sigma = model_robust1_Sigma,                 
												   lmRob_intercept = model_robust2_intercept, lmRob_slope = model_robust2_slope, lmRob_Rsquared = model_robust2_Rsquared, lmRob_Sigma = model_robust2_Sigma, 
                                                   lqs_intercept = model_robust3_intercept, lqs_slope = model_robust3_slope,
                                                   lsq_lms.only_intercept = model_robust4_intercept, lsq_lms.only_slope = model_robust4_slope,
                                                   Curve = Curve, Sample=Sample, Filename=current_file))
	
	

  } # end fragment ions
  
  
} # end file for loop

write.table(summary_data, "Summary_data.txt", sep="\t", col.names = TRUE, row.names=FALSE)	
write.table(individual_data, "individual_data.txt", sep="\t", col.names = TRUE, row.names=FALSE)	


## Step 3: Data plotting
library(plotly)
library(ggplot2)
library(gcookbook)
library(plyr)
library(devtools)
library(dplyr)
library(ggExtra)
library(MASS)

my_theme <- theme_bw() +
  theme(
    axis.text.x  = element_text(angle=0, vjust=0.5, size=12),   # rotate x-axis label and define label size
    # axis.ticks.x  = element_blank(), # x-axis ticks
    axis.text.y  = element_text(angle=0, vjust=0.5, size=12),   # rotate y-axis label and define label size
    axis.title.x = element_text(colour="black", size=16),  # x-axis title size
    axis.title.y = element_text(colour="black", size=16),  # y-axis title size
    plot.title = element_text(face="bold", colour="black", size=20, hjust=0.5, vjust=0.5),   # Main title size
    
    strip.text = element_text(size=8), 
    
    panel.grid.major.x = element_blank(), 
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(), 
    
    legend.key = element_rect(colour = NA, fill = 'transparent'), 
    legend.background = element_rect(colour = NA,  fill = "transparent"),
    # legend.position = c(0.2, 0.3),  
    legend.position = "none",  # to remove all legends
    # legend.title = element_blank(),
    legend.title = element_text(colour="black", size=16),
    legend.text = element_text(colour="black", size=14),
    legend.text.align = 0, 
    legend.box = NULL
  )


regression_eqn_r2 <- function (summary_data) { # obtain the regression equation
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(summary_data$rlm_intercept, digits = 2), 
                        b = format(summary_data$rlm_slope, digits = 2), 
                        r2 = format(summary_data$rlm_Rsquared, digits = 2)))
  as.character(as.expression(eq));
}

regression_eqn_r <- function (summary_data) { # obtain the regression equation
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)~"="~r2, 
                   list(a = format(summary_data$rlm_intercept, digits = 2), 
                        b = format(summary_data$rlm_slope, digits = 2), 
                        r2 = format( sqrt(summary_data$rlm_Rsquared) , digits = 2)))
  as.character(as.expression(eq));
}

setwd("C:\\Results\\AUDIT\\Regression_inputs")
dir.create(file.path("C:\\Results\\AUDIT", "Plots"), showWarnings = FALSE)
file_name = list.files(path = ".", pattern="*regression inputs\\.txt$")  # Find all txt files
summary_data <- read.table("Summary_data.txt", sep="\t", header=TRUE)

for (current_file in file_name) {
  # current_file <- file_name[3]
  for (j in 1:length(Sample.Types)) { if (grepl(Sample.Types[j], current_file, ignore.case=TRUE)) Sample <- Sample.Types[j]} # Determine sample type
  for (j in 1:length(Curve.Types)) { if (grepl(Curve.Types[j], current_file, ignore.case=TRUE)) Curve <- Curve.Types[j]} # Determine curve type

  summary_data_sub <- summary_data[summary_data$Filename == current_file,]
  df <- read.table(current_file, sep="\t", header=TRUE)
  
  if (Curve == "Fwd") {
	df <- df[!is.na(df$PAR_.L.H.),]
  } else {
	df <- df[!is.na(df$PAR_.H.L.),]
  }

  df <- df[!is.na(df$Meas_conc),]

  output.prefix <- data.frame(do.call('rbind', strsplit(as.character(current_file),'.txt',fixed=TRUE))) # File name prefix for outputs
  eq <- ddply(summary_data_sub,.(Peptide),regression_eqn_r)
  
  p <- ggplot() + 
    geom_point(df, mapping= aes(x= Theor_conc, y= Meas_conc), shape=21, fill="#377eb8", size=2, alpha=.5) + 
    geom_abline(data=summary_data_sub, aes(intercept=rlm_intercept, slope=rlm_slope), color="red", size=.75) + 
    # geom_hline(data=summary_data_sub, aes(yintercept=LoB_method3), linetype = "dashed", colour="#984ea3", size=.5) + 
	geom_hline(data=summary_data_sub, aes(yintercept=LOB), linetype = "dashed", colour="#984ea3", size=.5) + 
    geom_hline(data=summary_data_sub, aes(yintercept=LOD), linetype = "solid", colour="#4daf4a", size=.5) + 
    geom_text(data=eq,aes(x = 100, y = 175,label=V1), parse = TRUE, inherit.aes=FALSE, size=3) + 
    scale_x_continuous(limits = c(0, 200), breaks=(seq(0,200,by=100))) + 
    scale_y_continuous(limits = c(0, 200), breaks=(seq(0,200,by=100))) + 
    labs(title=paste(Curve, Sample, sep=": "), x=expression("Theoretical concentration (fmol/mg)"), y=expression("Measured concentation (fmol/mg)")) + 		  
    facet_wrap(~Peptide) + 
    my_theme
  
  filename <- paste(output.prefix[1,1], "Metabolic rev")
  ggsave(file.path("C:\\Results\\AUDIT\\Plots", paste(filename, "png", sep=".")), p, width = 20, height = 16, units="in", dpi=600)
  # ggsave(file.path("C:\\Results\\AUDIT\\Plots", paste(filename, "pdf", sep=".")), p, width = 20, height = 16, units="in")



  ## Plots at log10 scales
  df.abline <- summary_data_sub[, c(which(names(summary_data_sub)=="Peptide"), which(names(summary_data_sub)=="rlm_intercept"), which(names(summary_data_sub)=="rlm_slope"))]
  x.min <- 1E-2; x.max <- 200 # range of x values

  x.to.y <- function(x.min, x.max, rlm_intercept, rlm_slope) { # generate series of x and y values
    x <- seq(from=x.min, to=x.max, by=.1)
    y <- sapply(x, function(x, rlm_intercept, rlm_slope) y <- rlm_intercept + x*rlm_slope, rlm_intercept, rlm_slope)
    return(data.frame(x,y))
  }
  
  regession_curves <- NULL # data points for regression curves
  for (i in 1:nrow(df.abline)) {
    xy <- x.to.y(x.min, x.max, df.abline$rlm_intercept[i], df.abline$rlm_slope[i])
    regession_curves <- rbind(regession_curves, data.frame(Peptide=rep(df.abline$Peptide[i], nrow(xy)), Theor_conc=xy$x, Meas_conc=xy$y))
  }

  p.log10 <- ggplot() + 
    geom_point(df, mapping= aes(x= Theor_conc, y= Meas_conc), shape=21, fill="#377eb8", size=2, alpha=.5) + 
    geom_line(regession_curves, mapping= aes(x= Theor_conc, y= Meas_conc), colour="red", size=1) + 
    # geom_abline(data=summary_data_sub, aes(intercept=rlm_intercept, slope=rlm_slope), color="red", size=.75) + 
    # geom_hline(data=summary_data_sub, aes(yintercept=LoB_method3), linetype = "dashed", colour="#984ea3", size=.5) + 
	geom_hline(data=summary_data_sub, aes(yintercept=LOB), linetype = "dashed", colour="#984ea3", size=.5) + 
    geom_hline(data=summary_data_sub, aes(yintercept=LOD), linetype = "solid", colour="#4daf4a", size=.5) + 
    geom_text(data=eq,aes(x = 1, y = 50,label=V1), parse = TRUE, inherit.aes=FALSE, size=3) + 
    # scale_x_continuous(limits = c(0, 200), breaks=(seq(0,200,by=100))) + 
    # scale_y_continuous(limits = c(0, 200), breaks=(seq(0,200,by=100))) + 
    scale_x_log10(limits=c(x.min, x.max), breaks=c(0.01, 0.1, 1, 10, 100), labels=c(0.01, 0.1, 1, 10, 100)) + 
    scale_y_log10(limits=c(x.min, x.max), breaks=c(0.01, 0.1, 1, 10, 100), labels=c(0.01, 0.1, 1, 10, 100)) + 
    labs(title=paste(Curve, Sample, sep=": "), x=expression("Theoretical concentration (fmol/mg)"), y=expression("Measured concentation (fmol/mg)")) + 		  
    facet_wrap(~Peptide) + 
    my_theme
  # p.log10

  filename <- paste(output.prefix[1,1], "Metabolic rev", "log10")
  ggsave(file.path("C:\\Results\\AUDIT\\Plots", paste(filename, "png", sep=".")), p.log10, width = 20, height = 16, units="in", dpi=600)
  # ggsave(file.path("C:\\Results\\AUDIT\\Plots", paste(filename, "pdf", sep=".")), p.log10, width = 20, height = 16, units="in")

}  

