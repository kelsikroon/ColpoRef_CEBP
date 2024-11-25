# ppv - Table 2
library(DescTools)
library(stringr)

ppv.cases.function <- function(data, strategy){
  baseline.cyt <- data$baseline.cyt
  first.repeat.cyt <- data$first.repeat.cyt
  hpv.1618.pos <- data$hpv.1618pos
  hpv.genotype.pos <- data$hpv.genotype.pos
  second.round.hpv <- data$second.round.hpv
  first.round.cin3plus <- data$first.round.cin3plus
  ue.cin2plus.first.round <- ifelse((data$ue.detected==1 & data$hist.time/365.25 <= 4) | (data$first.round.cin2plus), T, F)
  
  if (strategy == "all"){
    baseline.referred <- data$id 
    repeat.test.referred <- c()
    second.round.hpv.pos <- c()
    
  }else if (strategy == "partial"){
    baseline.referred <- data$id[which(hpv.1618.pos)]
    repeat.test.referred <- data$id[which(!hpv.1618.pos & first.repeat.cyt=='ASC-US/LSIL+')]
    second.round.hpv.pos <- data$id[which(!hpv.1618.pos & first.repeat.cyt=='NILM' & !ue.cin2plus.first.round & second.round.hpv =='+')]
    second.round.test <- data$id[which(!hpv.1618.pos & first.repeat.cyt=='NILM' & !ue.cin2plus.first.round & second.round.hpv !='')]

  }else if (strategy =="extended"){
    baseline.referred <- data$id[which(hpv.genotype.pos)]
    repeat.test.referred <- data$id[which(! hpv.genotype.pos & first.repeat.cyt=='ASC-US/LSIL+')]
    second.round.hpv.pos <- data$id[which(!hpv.genotype.pos & first.repeat.cyt=='NILM' & !ue.cin2plus.first.round &  second.round.hpv =='+')]
    second.round.test <- data$id[which(!hpv.genotype.pos & first.repeat.cyt=='NILM' & !ue.cin2plus.first.round & second.round.hpv !='')]

  }else if (strategy == "extended.repeat"){
    baseline.referred <- c()
    repeat.test.referred <- data$id[which(hpv.genotype.pos & first.repeat.cyt =='ASC-US/LSIL+')]
    second.round.hpv.pos <-  data$id[which(!ue.cin2plus.first.round & second.round.hpv =='+' & (hpv.genotype.pos & first.repeat.cyt =='NILM') | (! hpv.genotype.pos))]
    second.round.test <-  data$id[which(!ue.cin2plus.first.round & second.round.hpv !='' & ((hpv.genotype.pos & first.repeat.cyt =='NILM') | (! hpv.genotype.pos )))]
    
  }else if (strategy =='repeat'){
    baseline.referred <- c()
    repeat.test.referred <- data$id[which(first.repeat.cyt =='ASC-US/LSIL+')]
    second.round.hpv.pos <-  data$id[which((first.repeat.cyt =='NILM' & !ue.cin2plus.first.round & second.round.hpv =='+') )]
    second.round.test <-  data$id[which((first.repeat.cyt =='NILM' & !ue.cin2plus.first.round & second.round.hpv !='') )]
  }

  cin3plus.baseline <- sum(data[data$id %in% baseline.referred,]$first.round.cin3plus)
  cin3plus.repeat <-  sum(data[data$id %in% repeat.test.referred,]$first.round.cin3plus) 
  
  
  cin2plus.baseline <- sum(data[data$id %in% baseline.referred,]$first.round.cin2plus)
  cin2plus.repeat <-  sum(data[data$id %in% repeat.test.referred,]$first.round.cin2plus) 
  
  
  first.round.referred <- c(baseline.referred, repeat.test.referred)
  anytime.referred <- c(first.round.referred, second.round.hpv)
  results <- list(baseline.referred, repeat.test.referred,first.round.referred, cin3plus.baseline, 
                  cin3plus.repeat, anytime.referred, cin2plus.baseline, cin2plus.repeat)
  names(results) <- c("baseline.referred", "repeat.test.referred", "first.round.referred", "cin3plus.baseline", 
                      "cin3plus.repeat", "anytime.referred", "cin2plus.baseline", "cin2plus.repeat")
  return(results)
}


ascus.strategy <- c(rep("all", 5), rep("extended", 4), rep("partial", 3), rep("repeat", 2))
nilm.strategy <- c("all", "extended", "partial", "repeat", "extended.repeat", "extended", "partial", "repeat", "extended.repeat", 
                   "partial", "repeat", "extended.repeat", "repeat", "extended.repeat")

# ppv first round:
ppv.1st.round.function <- function(data, ascus.strategy, nilm.strategy, endpoint='CIN3+'){
  ascus.results <- ppv.cases.function(data[data$baseline.cyt=="ASC-US/LSIL", ], ascus.strategy)
  nilm.results <- ppv.cases.function(data[data$baseline.cyt=="NILM", ], nilm.strategy)
  
  if (endpoint=='CIN3+'){
    ppv.baseline.bmd <- BinomCI(ascus.results$cin3plus.baseline, length(ascus.results$baseline.referred), method='wilson')
    ppv.baseline.norm <- BinomCI(nilm.results$cin3plus.baseline, length(nilm.results$baseline.referred), method='wilson')
    ppv.repeat.bmd <- BinomCI(ascus.results$cin3plus.repeat, length(ascus.results$repeat.test.referred), method='wilson')
    ppv.repeat.norm <- BinomCI(nilm.results$cin3plus.repeat, length(nilm.results$repeat.test.referred), method='wilson')
  }else{
    ppv.baseline.bmd <- BinomCI(ascus.results$cin2plus.baseline, length(ascus.results$baseline.referred), method='wilson')
    ppv.baseline.norm <- BinomCI(nilm.results$cin2plus.baseline, length(nilm.results$baseline.referred), method='wilson')
    ppv.repeat.bmd <- BinomCI(ascus.results$cin2plus.repeat, length(ascus.results$repeat.test.referred), method='wilson')
    ppv.repeat.norm <- BinomCI(nilm.results$cin2plus.repeat, length(nilm.results$repeat.test.referred), method='wilson')
  }

  
  results <- list(ppv.baseline.bmd[1], ppv.baseline.bmd[2], ppv.baseline.bmd[3], 
                  ppv.baseline.norm[1], ppv.baseline.norm[2], ppv.baseline.norm[3], 
                  ppv.repeat.bmd[1], ppv.repeat.bmd[2], ppv.repeat.bmd[3],
                  ppv.repeat.norm[1], ppv.repeat.norm[2], ppv.repeat.norm[3])
  names(results) <- c("ppv.baseline.bmd", "ppv.baseline.bmd.low", "ppv.baseline.bmd.upp", 
                      "ppv.baseline.norm", "ppv.baseline.norm.low", "ppv.baseline.norm.upp",
                      "ppv.repeat.bmd", "ppv.repeat.bmd.low", "ppv.repeat.bmd.upp",
                      "ppv.repeat.norm", "ppv.repeat.norm.low", "ppv.repeat.norm.upp")
  return(results)
}

ppv.data1 <- ref.data.cens[ref.data.cens$baseline.cyt == "ASC-US/LSIL" | (ref.data.cens$grp =='i' & ref.data.cens$baseline.cyt != "HSIL" ),]

ppv.data2 <- ref.data.cens[ref.data.cens$baseline.cyt == "ASC-US/LSIL" | 
                            (ref.data.cens$baseline.cyt =="NILM" & ref.data.cens$hpv.genotype.pos & ref.data.cens$grp =='i') |
                            (ref.data.cens$baseline.cyt =="NILM" & ! ref.data.cens$hpv.genotype.pos & ref.data.cens$grp =='c'),]


results.ppv <- list()
for (i in 1:14){
  ascus <- ascus.strategy[i]
  nilm <- nilm.strategy[i]
  if (str_detect(nilm, "extended.repeat")){
    dat <- ppv.data2
  }else {
    dat <- ppv.data1
  }
  
  results.ppv[[i]] <- ppv.1st.round.function(dat, ascus, nilm, "CIN3+")
}

# ppv second round: remove UE/CIN2+ first round cases from the data or had missing test results 
second.round.data <- ppv.data1[!ppv.data1$ue.cin2plus.first.round & ppv.data1$second.round.hpv!='', ]

ppv.second.rond.norm <- round(BinomCI(sum(second.round.data$second.round.cin3plus[second.round.data$baseline.cyt=='NILM']), 
                                      sum(second.round.data$second.round.hpv[second.round.data$baseline.cyt=='NILM'] =='+'), method='wilson')*100, 2)

ppv.second.round.bmd <- round(BinomCI(sum(second.round.data$second.round.cin3plus[second.round.data$baseline.cyt=='ASC-US/LSIL']), 
                                      sum(second.round.data$second.round.hpv[second.round.data$baseline.cyt=='ASC-US/LSIL'] =='+'), method='wilson')*100, 2)



# Table 2: ppv of different strategies (uses control arm for 7-types negative normal cytology)
ppv.results <- round(data.frame(matrix(unlist(results.ppv), nrow=length(results.ppv), byrow = T))*100, 2)
colnames(ppv.results) <- names(results.ppv[[1]])
ppv.results$ascus.strategy <- factor(ascus.strategy)
ppv.results$nilm.strategy <- factor(nilm.strategy, levels=c("all", "extended", "partial", "repeat", "extended.repeat"))

ppv.results %>%   as.data.frame() %>% mutate_if(is.numeric, round, 1) %>% 
  mutate(ppv.baseline.bmd = paste0(ppv.baseline.bmd, " (", ppv.baseline.bmd.low, " - ", ppv.baseline.bmd.upp, ")"), 
         ppv.repeat.bmd = paste0(ppv.repeat.bmd, " (", ppv.repeat.bmd.low, " - ", ppv.repeat.bmd.upp, ")"),
         ppv.baseline.norm = paste0(ppv.baseline.norm, " (", ppv.baseline.norm.low, " - ", sprintf("%.1f", ppv.baseline.norm.upp), ")"), 
         ppv.repeat.norm = paste0(sprintf("%.1f", ppv.repeat.norm), " (", ppv.repeat.norm.low, " - ", ppv.repeat.norm.upp, ")"),
         ppv.2nd.bmd = paste0(ppv.second.round.bmd[1], " (", ppv.second.round.bmd[2], " - ", ppv.second.round.bmd[3], ")"),
         ppv.2nd.norm = paste0(ppv.second.rond.norm[1], " (", ppv.second.rond.norm[2], " - ", ppv.second.rond.norm[3], ")")) %>% 
  select(ascus.strategy,   nilm.strategy, ppv.baseline.bmd, ppv.baseline.norm, ppv.repeat.bmd, ppv.repeat.norm, ppv.2nd.bmd, ppv.2nd.norm) %>% 
  write.csv(., "~/Desktop/PhD/Projects/Colposcopy Referral Rates/github/results/Table2_PPVresults.csv")



# Table S2: PPV of 14 strategies for endpoint CIN2+

results.ppvCIN2 <- list()
for (i in 1:14){
  ascus <- ascus.strategy[i]
  nilm <- nilm.strategy[i]
  if (str_detect(nilm, "extended.repeat")){
    dat <- ppv.data2
  }else {
    dat <- ppv.data1
  }
  
  results.ppvCIN2[[i]] <- ppv.1st.round.function(dat, ascus, nilm, "CIN2+")
}

#------------ CIN2+
# ppv second round: remove UE/CIN2+ first round cases from the data or had missing test results 
ppv.second.rond.normCIN2 <- round(BinomCI(sum(second.round.data$second.round.cin2plus[second.round.data$baseline.cyt=='NILM']), 
                                      sum(second.round.data$second.round.hpv[second.round.data$baseline.cyt=='NILM'] =='+'), method='wilson')*100, 2)

ppv.second.round.bmdCIN2 <- round(BinomCI(sum(second.round.data$second.round.cin2plus[second.round.data$baseline.cyt=='ASC-US/LSIL']), 
                                      sum(second.round.data$second.round.hpv[second.round.data$baseline.cyt=='ASC-US/LSIL'] =='+'), method='wilson')*100, 2)



# Table 3: ppv of different strategies (uses control arm for 7-types negative normal cytology)
ppv.resultsCIN2 <- round(data.frame(matrix(unlist(results.ppvCIN2), nrow=length(results.ppvCIN2), byrow = T))*100, 2)
colnames(ppv.resultsCIN2) <- names(results.ppvCIN2[[1]])
ppv.resultsCIN2$ascus.strategy <- factor(ascus.strategy)
ppv.resultsCIN2$nilm.strategy <- factor(nilm.strategy, levels=c("all", "extended", "partial", "repeat", "extended.repeat"))

ppv.resultsCIN2 %>%   as.data.frame() %>% mutate_if(is.numeric, round, 1) %>% 
  mutate(ppv.baseline.bmd = paste0(ppv.baseline.bmd, " (", ppv.baseline.bmd.low, " - ", ppv.baseline.bmd.upp, ")"), 
         ppv.repeat.bmd = paste0(ppv.repeat.bmd, " (", ppv.repeat.bmd.low, " - ", ppv.repeat.bmd.upp, ")"),
         ppv.baseline.norm = paste0(ppv.baseline.norm, " (", ppv.baseline.norm.low, " - ", sprintf("%.1f", ppv.baseline.norm.upp), ")"), 
         ppv.repeat.norm = paste0(sprintf("%.1f", ppv.repeat.norm), " (", ppv.repeat.norm.low, " - ", ppv.repeat.norm.upp, ")"),
         ppv.2nd.bmd = paste0(ppv.second.rond.normCIN2[1], " (", ppv.second.rond.normCIN2[2], " - ", ppv.second.rond.normCIN2[3], ")"),
         ppv.2nd.norm = paste0(ppv.second.round.bmdCIN2[1], " (", ppv.second.round.bmdCIN2[2], " - ", ppv.second.round.bmdCIN2[3], ")")) %>% 
  select(ascus.strategy,   nilm.strategy, ppv.baseline.bmd, ppv.baseline.norm, ppv.repeat.bmd, ppv.repeat.norm, ppv.2nd.bmd, ppv.2nd.norm) %>% 
  write.csv(., "~/Desktop/PhD/Projects/Colposcopy Referral Rates/github/results//TableS2_PPVresults_CIN2plus.csv")

