# npv - bootstrap confidence intervals - Table 3 and Table S3
# intervention arm: normals and take double
# intervention and contorl arm: BMD
library(dplyr)

npv.function <- function(data, ascus.strategy, nilm.strategy, reps=1000, endpoint='CIN3+'){
  total <- dim(data)[1]
  
  if (endpoint =='CIN3+'){
    cases.1 <- 'first.round.cin3plus'
    cases.all <- 'cin3plus.cases'
  }else{
    cases.1 <- 'first.round.cin2plus'
    cases.all <- 'cin2plus.cases'
  }
  
  num.not.referred <- rep(NA, 1000)
  num.no.cin.1st.round <- rep(NA, 1000)
  num.no.cin.ever <- rep(NA, 1000)
  
  for (i in 1:reps){
    boot.indxs <- sample(1:total, total, replace=T)
    bootstrap.data <- data[boot.indxs,]
    
    ascus.results <- ppv.cases.function(bootstrap.data[bootstrap.data$baseline.cyt=="ASC-US/LSIL", ], ascus.strategy)
    nilm.results <- ppv.cases.function(bootstrap.data[bootstrap.data$baseline.cyt=="NILM", ], nilm.strategy)
    not.referred <- bootstrap.data$id[!bootstrap.data$id %in% c(ascus.results$first.round.referred, nilm.results$first.round.referred)]
    
    num.not.referred[i] <- length(not.referred)
  
    num.no.cin.1st.round[i] <- length(bootstrap.data$id[bootstrap.data$id %in% not.referred & bootstrap.data[[cases.1]] ==0])
    num.no.cin.ever[i] <- length(bootstrap.data$id[bootstrap.data$id %in% not.referred & bootstrap.data[[cases.all]] ==0])
  }

  npv.first.round <- num.no.cin.1st.round/num.not.referred
  npv.first.round.est <- mean(npv.first.round)
  npv.first.round.var <- var(npv.first.round)
  npv.first.round.ci <- npv.first.round.est + c(-1, 1)*1.96*sqrt(npv.first.round.var)
  
  npv.cumulative <- num.no.cin.ever/num.not.referred
  npv.cumulative.est <- mean(npv.cumulative)
  npv.cumulative.var <- var(npv.cumulative)
  npv.cumulative.ci <- npv.cumulative.est + c(-1, 1)*1.96*sqrt(npv.cumulative.var)
  
  results <- list(npv.first.round = npv.first.round.est, npv.first.round.lower = npv.first.round.ci[1], npv.first.round.upper = npv.first.round.ci[2],
                  npv.cumulative = npv.cumulative.est, npv.cumulative.lower = npv.cumulative.ci[1], npv.cumulative.upper = npv.cumulative.ci[2])
    
  return(results)
}


ascus.strategy <- c(rep("all", 5), rep("extended", 4), rep("partial", 3), rep("repeat", 2))
nilm.strategy <- c("all", "extended", "partial", "repeat", "extended.repeat", "extended", "partial", "repeat", "extended.repeat",
                   "partial", "repeat", "extended.repeat", "repeat", "extended.repeat")

npv.data <- ref.data.cens[ref.data.cens$baseline.cyt == "ASC-US/LSIL" | (ref.data.cens$grp =='i' & ref.data.cens$baseline.cyt != "HSIL" ),]

npv.double.data <- data.frame(rbind(npv.data, npv.data[npv.data$baseline.cyt=='NILM',]))

results.npv <- list()
for (i in 1:14){
  ascus <- ascus.strategy[i]
  nilm <- nilm.strategy[i]
  results.npv[[i]] <- npv.function(npv.double.data, ascus, nilm, "CIN3+")
}

{
  # Table 3:  NPV for endpoint CIN3+ of 14 different strategies
  npv.results <- round(data.frame(matrix(unlist(results.npv), nrow=length(results.npv), byrow = T))*100, 2)
  colnames(npv.results) <- names(results.npv[[1]])
  npv.results$ascus.strategy <- factor(ascus.strategy)
  npv.results$nilm.strategy <- factor(nilm.strategy, levels=c("all", "extended", "partial", "repeat", "extended.repeat"))
  npv.results$npv.first.round.upper[npv.results$npv.first.round.upper > 100] <- 100.00
  npv.results %>%
    mutate(NPV.1st.round = paste0(sprintf("%.1f", npv.first.round), " (", sprintf("%.1f", npv.first.round.lower), " - ", sprintf("%.1f", npv.first.round.upper), ")"),
           NPV.two.rounds = paste0(sprintf("%.1f", npv.cumulative), " (", sprintf("%.1f", npv.cumulative.lower), " - ", sprintf("%.1f", npv.cumulative.upper), ")")) %>%
    select(ascus.strategy,   nilm.strategy, NPV.1st.round, NPV.two.rounds) %>%
   write.csv(., "~/Desktop/PhD/Projects/Colposcopy Referral Rates/github/Table3_NPVresults.csv")

}

#---------
results.npvCIN2 <- list()
for (i in 1:14){
  ascus <- ascus.strategy[i]
  nilm <- nilm.strategy[i]
  results.npvCIN2[[i]] <- npv.function(npv.double.data, ascus, nilm, endpoint= "CIN2+")
}

{
  # Table S3:  NPV for endpoint CIN3+ of 14 different strategies 
  npv.resultsCIN2 <- round(data.frame(matrix(unlist(results.npvCIN2), nrow=length(results.npvCIN2), byrow = T))*100, 2)
  colnames(npv.resultsCIN2) <- names(results.npvCIN2[[1]])
  npv.resultsCIN2$ascus.strategy <- factor(ascus.strategy)
  npv.resultsCIN2$nilm.strategy <- factor(nilm.strategy, levels=c("all", "extended", "partial", "repeat", "extended.repeat"))
  npv.resultsCIN2$npv.first.round.upper[npv.resultsCIN2$npv.first.round.upper > 100] <- 100.00
  npv.resultsCIN2 %>% 
    mutate(NPV.1st.round = paste0(sprintf("%.1f", npv.first.round), " (", sprintf("%.1f", npv.first.round.lower), " - ", sprintf("%.1f", npv.first.round.upper), ")"),
           NPV.two.rounds = paste0(sprintf("%.1f", npv.cumulative), " (", sprintf("%.1f", npv.cumulative.lower), " - ", sprintf("%.1f", npv.cumulative.upper), ")")) %>% 
    select(ascus.strategy,   nilm.strategy, NPV.1st.round, NPV.two.rounds) %>% 
    write.csv(., "~/Desktop/PhD/Projects/Colposcopy Referral Rates/github/TableS3_NPVresults_CIN2plus.csv")
  
}
