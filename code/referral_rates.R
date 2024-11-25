# referral rates calculation - Figure 2 and Table S1
 
#------------------------------------------
# Function to calculate the posterior probability using sampling from beta prior 
#------------------------------------------
post.prob <- function(props){
  x <- props[[1]][1]
  n <- props[[1]][2]
  
  sample <- rbeta(100000, x + 0.5, n-x+0.5)
  if (length(props)>1){
    for (i in 2:length(props)){
      x <- props[[i]][1]
      n <- props[[i]][2]
      sample <- sample*rbeta(100000, x + 0.5, n-x+0.5)
    }
  }
  
  return(sample)
  #return(c(mean =mean(sample), quantile(sample, 0.025), quantile(sample, 0.975)))
}

#------------------------------------------
# Function to calculate referral rate by multiplying proportions at each test for each type of cytology & genotype result 
#------------------------------------------
# possible strategies are:
# (1) all, (2) partial, (3) extended, (4) partial.repeat, (5) extended.repeat, (6) repeat 
strategy.function <- function(data, total, strategy){
  baseline.cyt <- data$baseline.cyt
  first.repeat.cyt <- data$first.repeat.cyt
  hpv.1618.pos <- data$hpv.1618pos
  hpv.genotype.pos <- data$hpv.genotype.pos
  second.round.hpv <- data$second.round.hpv
  subset.total <- dim(data)[1]
  first.round.cin3plus <- data$first.round.cin3plus
  
  if (strategy == "all"){
    first.round.referred <- data$id
    second.round.hpv.pos <- c()
    
    catchup.2nd.round <- 0
    
    baseline.samples <- post.prob(list(c(subset.total, total)))
    repeat.samples <- 0
    second.round.samples <- 0

    
  }else if (strategy == "partial"){
    first.round.referred <- data$id[which(hpv.1618.pos | (! hpv.1618.pos & first.repeat.cyt=='ASC-US/LSIL+'))]
    #second.round.hpv.pos <- data$id[which(!hpv.1618.pos & first.repeat.cyt=='NILM' & second.round.hpv =='+')]
    
    catchup.2nd.round <- sum(! data$id %in% first.round.referred & first.round.cin3plus==1)
    prob.2nd.round.test <- sum(!hpv.1618.pos & first.repeat.cyt=='NILM' & second.round.hpv !='')/sum(!hpv.1618.pos & first.repeat.cyt=='NILM')
    
    baseline.samples <- post.prob(list(c(subset.total, total), c(sum(hpv.1618.pos), subset.total)))
    
    repeat.samples <- post.prob(list(c(subset.total, total), c(sum(! hpv.1618.pos), subset.total), 
                                     c(sum(! hpv.1618.pos & first.repeat.cyt=='ASC-US/LSIL+'), sum(! hpv.1618.pos & first.repeat.cyt!='Missing')))) 
    # second round referral rate = those who are HPV+ at second round + proportion of missed CIN3+ cases
    second.round.samples <- post.prob(list(c(subset.total, total), c(sum(! hpv.1618.pos), subset.total),
                                           c(sum(!hpv.1618.pos & first.repeat.cyt=='NILM'), sum(!hpv.1618.pos & first.repeat.cyt!='Missing')),
                                           c(sum(!hpv.1618.pos & first.repeat.cyt=='NILM' & second.round.hpv =='+') + catchup.2nd.round*prob.2nd.round.test, 
                                             sum(!hpv.1618.pos & first.repeat.cyt=='NILM' & second.round.hpv!='') + catchup.2nd.round*prob.2nd.round.test))) 
    
  }else if (strategy =="extended"){
    first.round.referred <- data$id[which(hpv.genotype.pos | (! hpv.genotype.pos & first.repeat.cyt=='ASC-US/LSIL+'))]
    second.round.hpv.pos <- data$id[which(!hpv.genotype.pos & first.repeat.cyt=='NILM' & second.round.hpv =='+')]
    
    catchup.2nd.round <- sum(! data$id %in% first.round.referred & first.round.cin3plus==1)
    prob.2nd.round.test <- sum(!hpv.genotype.pos & first.repeat.cyt=='NILM' &  second.round.hpv !='')/sum(!hpv.genotype.pos & first.repeat.cyt=='NILM')
    
    baseline.samples <- post.prob(list(c(subset.total, total), c(sum(hpv.genotype.pos), subset.total)))
    
    repeat.samples <- post.prob(list(c(subset.total, total), c(sum(! hpv.genotype.pos), subset.total), 
                                     c(sum(! hpv.genotype.pos & first.repeat.cyt=='ASC-US/LSIL+'), sum(! hpv.genotype.pos & first.repeat.cyt!='Missing')) ))
    
    # second round referral rate = those who are HPV+ at second round + proportion of missed CIN3+ cases
    second.round.samples <- post.prob(list(c(subset.total, total), c(sum(!hpv.genotype.pos), subset.total),
                                           c(sum(!hpv.genotype.pos & first.repeat.cyt=='NILM'), sum(!hpv.genotype.pos & first.repeat.cyt!='Missing')),
                                           c(sum(!hpv.genotype.pos & first.repeat.cyt=='NILM' & second.round.hpv =='+') + catchup.2nd.round*prob.2nd.round.test,
                                             sum(!hpv.genotype.pos & first.repeat.cyt=='NILM' & second.round.hpv!='') + catchup.2nd.round*prob.2nd.round.test))) #+ 
    
  }else if (strategy == "extended.repeat"){
    first.round.referred <- data$id[which(hpv.genotype.pos & first.repeat.cyt =='ASC-US/LSIL+')]
    second.round.hpv.pos <-  data$id[which((hpv.genotype.pos & first.repeat.cyt =='NILM' & second.round.hpv =='+') | (! hpv.genotype.pos & second.round.hpv =='+'))]
    
    catchup.2nd.round1 <- sum(hpv.genotype.pos & first.repeat.cyt =='NILM' & first.round.cin3plus==1)
    prob.2nd.round.test1 <- sum(hpv.genotype.pos & first.repeat.cyt =='NILM' &  second.round.hpv !='')/sum(hpv.genotype.pos & first.repeat.cyt =='NILM')
    
    catchup.2nd.round2 <- sum(! hpv.genotype.pos & first.round.cin3plus==1)
    prob.2nd.round.test2 <- sum(! hpv.genotype.pos &  second.round.hpv !='')/sum(! hpv.genotype.pos)
    
    catchup.2nd.round <- catchup.2nd.round1 + catchup.2nd.round2
    
    baseline.samples <- 0 # repeat referral strategy so baseline rate =0
    repeat.samples <- post.prob(list(c(subset.total, total), c(sum(hpv.genotype.pos), subset.total), 
                                     c(sum(hpv.genotype.pos & first.repeat.cyt =='ASC-US/LSIL+'), sum(hpv.genotype.pos & first.repeat.cyt!='Missing'))))
    
    # second round referral rate = those who are HPV+ at second round + proportion of missed CIN3+ cases
    # referral in second round after extended positive:
    # (1) HPVh7+ --> Repeat cyt = Normal --> second round HPV+
    # (2) HPVh7- --> second round HPV+ 
    second.round.samples <- post.prob(list(c(subset.total, total), c(sum(hpv.genotype.pos), subset.total),
                                           c(sum(hpv.genotype.pos & first.repeat.cyt =='NILM'), sum(hpv.genotype.pos & first.repeat.cyt!='Missing')),
                                           c(sum(hpv.genotype.pos & first.repeat.cyt =='NILM' & second.round.hpv =='+') + catchup.2nd.round1*prob.2nd.round.test1, 
                                             sum(hpv.genotype.pos & first.repeat.cyt =='NILM' & second.round.hpv !='')+ catchup.2nd.round1*prob.2nd.round.test1))) +
      post.prob(list(c(subset.total, total), c(sum(! hpv.genotype.pos), subset.total), 
                     c(sum(! hpv.genotype.pos & second.round.hpv =='+') + catchup.2nd.round2*prob.2nd.round.test2, 
                       sum(! hpv.genotype.pos & second.round.hpv !='') + catchup.2nd.round2*prob.2nd.round.test2))) 
    
  }else if (strategy =='repeat'){
    first.round.referred <- data$id[which(first.repeat.cyt =='ASC-US/LSIL+')]
    second.round.hpv.pos <-  data$id[which((first.repeat.cyt =='NILM' & second.round.hpv =='+') )]

    catchup.2nd.round <- sum(! data$id %in% first.round.referred & first.round.cin3plus==1)
    prob.2nd.round.test <- sum(first.repeat.cyt=='NILM' &  second.round.hpv !='')/sum(first.repeat.cyt=='NILM')
    
    baseline.samples <- 0 # repeat referral strategy so baseline rate =0
    repeat.samples <-  post.prob(list(c(subset.total, total), c(sum(first.repeat.cyt =='ASC-US/LSIL+'), sum(first.repeat.cyt != 'Missing'))))
    
    # second round referral rate = those who are HPV+ at second round + proportion of missed CIN3+ cases
    second.round.samples <- post.prob(list(c(subset.total, total), c(sum(first.repeat.cyt=='NILM'), sum(first.repeat.cyt!='Missing')),
                                           c(sum(first.repeat.cyt=='NILM' & second.round.hpv=='+')  + catchup.2nd.round*prob.2nd.round.test,
                                             sum(first.repeat.cyt=='NILM' & second.round.hpv!='')  + catchup.2nd.round*prob.2nd.round.test))) 
  }


  results <- list(baseline.samples, repeat.samples, second.round.samples, catchup.2nd.round)
  names(results) <- c("baseline.samples", "repeat.samples", "second.round.samples", "catchup.2nd.round")
  
  return(results)
}


#------------------------------------------
# Function to combine referral rates of baseline cytology results  
#------------------------------------------
referral.function <- function(data, ascus.strategy, nilm.strategy, diff.calc=T){
  total <- dim(data)[1]
  ascus.results <- strategy.function(data[data$baseline.cyt=="ASC-US/LSIL", ], total, ascus.strategy)
  nilm.results <- strategy.function(data[data$baseline.cyt=="NILM", ], total, nilm.strategy)
  
  baseline.samples <-  ascus.results$baseline.samples + nilm.results$baseline.samples 
  baseline.rate <- c(mean(baseline.samples), quantile(baseline.samples, 0.025), quantile(baseline.samples, 0.975))
  
  repeat.samples <-  ascus.results$repeat.samples + nilm.results$repeat.samples
  repeat.rate <- c(mean(repeat.samples), quantile(repeat.samples, 0.025), quantile(repeat.samples, 0.975))
  
  second.round.samples <-  ascus.results$second.round.samples + nilm.results$second.round.samples
  second.round.rate <- c(mean(second.round.samples), quantile(second.round.samples, 0.025), quantile(second.round.samples, 0.975))
  
  first.round.rate <- c(mean(baseline.samples + repeat.samples), 
                        quantile(baseline.samples + repeat.samples, 0.025), 
                        quantile(baseline.samples + repeat.samples, 0.975))
  
  overall.hpv.rate <- c(mean(baseline.samples + repeat.samples + second.round.samples), 
                        quantile(baseline.samples + repeat.samples + second.round.samples, 0.025),
                        quantile(baseline.samples + repeat.samples + second.round.samples, 0.975))
  
  catchup.2nd.round <- ascus.results$catchup.2nd.round + nilm.results$catchup.2nd.round
  
  ### difference between reference at 3 time points + 95% CI 
  if (diff.calc){
    baseline.diff <- c(mean(baseline.samples - (ascus.ref$baseline.samples + nilm.ref$baseline.samples)),
                       quantile(baseline.samples - (ascus.ref$baseline.samples + nilm.ref$baseline.samples), 0.025),
                       quantile(baseline.samples - (ascus.ref$baseline.samples + nilm.ref$baseline.samples), 0.975))
    
    first.round.diff <- c(mean(baseline.samples+repeat.samples - 
                                 (ascus.ref$baseline.samples + nilm.ref$baseline.samples + ascus.ref$repeat.samples + nilm.ref$repeat.samples)),
                          quantile(baseline.samples+repeat.samples   - 
                                     (ascus.ref$baseline.samples + nilm.ref$baseline.samples + ascus.ref$repeat.samples + nilm.ref$repeat.samples), 0.025),
                          quantile(baseline.samples+repeat.samples  - 
                                     (ascus.ref$baseline.samples + nilm.ref$baseline.samples + ascus.ref$repeat.samples + nilm.ref$repeat.samples), 0.975))
    
    two.round.diff <-  c(mean(baseline.samples+repeat.samples + second.round.samples- 
                                (ascus.ref$baseline.samples + nilm.ref$baseline.samples + ascus.ref$repeat.samples + nilm.ref$repeat.samples +
                                   ascus.ref$second.round.samples + nilm.ref$second.round.samples)),
                         quantile(baseline.samples+repeat.samples + second.round.samples - 
                                    (ascus.ref$baseline.samples + nilm.ref$baseline.samples + ascus.ref$repeat.samples + nilm.ref$repeat.samples +
                                       ascus.ref$second.round.samples + nilm.ref$second.round.samples), 0.025),
                         quantile(baseline.samples+repeat.samples + second.round.samples - 
                                    (ascus.ref$baseline.samples + nilm.ref$baseline.samples + ascus.ref$repeat.samples + nilm.ref$repeat.samples +
                                       ascus.ref$second.round.samples + nilm.ref$second.round.samples), 0.975))
    
    results <- list(baseline.rate[1], baseline.rate[2], baseline.rate[3], repeat.rate[1], repeat.rate[2], repeat.rate[3],
                    second.round.rate[1], second.round.rate[2], second.round.rate[3], overall.hpv.rate[1], overall.hpv.rate[2], overall.hpv.rate[3], 
                    first.round.rate[1],  first.round.rate[2],  first.round.rate[3], baseline.diff[1], baseline.diff[2], baseline.diff[3],
                    first.round.diff[1], first.round.diff[2], first.round.diff[3], two.round.diff[1], two.round.diff[2], two.round.diff[3], catchup.2nd.round)
    names(results) <- c( "baseline.rate","baseline.rate.lower", "baseline.rate.upper","repeat.rate", "repeat.rate.lower", "repeat.rate.upper",
                         "second.round.rate", "second.round.rate.lower", "second.round.rate.upper", "overall.hpv.rate", "overall.hpv.rate.lower", "overall.hpv.rate.upper",
                         "first.round.rate", "first.round.rate.lower", "first.round.rate.upper", "baseline.diff", "baseline.diff.lower", "baseline.diff.upper",
                         "first.round.diff", "first.round.diff.lower", "first.round.diff.upper", "two.round.diff", "two.round.diff.lower", "two.round.diff.upper", "catchup.2nd.round")
  }
  else{
    results <- list(baseline.rate[1], baseline.rate[2], baseline.rate[3], repeat.rate[1], repeat.rate[2], repeat.rate[3],
                    second.round.rate[1], second.round.rate[2], second.round.rate[3], overall.hpv.rate[1], overall.hpv.rate[2], overall.hpv.rate[3], 
                    first.round.rate[1],  first.round.rate[2],  first.round.rate[3], catchup.2nd.round)
    names(results) <- c( "baseline.rate","baseline.rate.lower", "baseline.rate.upper", "repeat.rate", "repeat.rate.lower", "repeat.rate.upper",
                         "second.round.rate", "second.round.rate.lower", "second.round.rate.upper", "overall.hpv.rate", "overall.hpv.rate.lower", "overall.hpv.rate.upper",
                         "first.round.rate", "first.round.rate.lower", "first.round.rate.upper", "catchup.2nd.round")
  }
  return(results)
}


diff.conservative <- function(data1, ascus.strategy1, nilm.strategy1, data2, ascus.strategy2, nilm.strategy2){
  total1 <- dim(data1)[1]
  ascus.results1 <- strategy.function(data1[data1$baseline.cyt=="ASC-US/LSIL", ], total1, ascus.strategy1)
  nilm.results1 <- strategy.function(data1[data1$baseline.cyt=="NILM", ], total1, nilm.strategy1)
  
  baseline.samples1 <-  ascus.results1$baseline.samples + nilm.results1$baseline.samples 
  repeat.samples1 <-  ascus.results1$repeat.samples + nilm.results1$repeat.samples
  second.round.samples1 <-  ascus.results1$second.round.samples + nilm.results1$second.round.samples
  
  total2 <- dim(data2)[1]
  ascus.results2 <- strategy.function(data2[data2$baseline.cyt=="ASC-US/LSIL", ], total2, ascus.strategy2)
  nilm.results2 <- strategy.function(data2[data2$baseline.cyt=="NILM", ], total2, nilm.strategy2)
  
  baseline.samples2 <-  ascus.results2$baseline.samples + nilm.results2$baseline.samples 
  repeat.samples2 <-  ascus.results2$repeat.samples + nilm.results2$repeat.samples
  second.round.samples2 <-  ascus.results2$second.round.samples + nilm.results2$second.round.samples
  
  two.round.diff <-  c(mean(baseline.samples1+repeat.samples1 + second.round.samples1- 
                              (baseline.samples2 + repeat.samples2 + second.round.samples2)),
                       quantile(baseline.samples1+repeat.samples1 + second.round.samples1 - 
                                  (baseline.samples2 + repeat.samples2 + second.round.samples2), 0.025),
                       quantile(baseline.samples1+repeat.samples1 + second.round.samples1 - 
                                  (baseline.samples2 + repeat.samples2 + second.round.samples2), 0.975))
  results <- list(two.round.diff[1], two.round.diff[2], two.round.diff[3])
  names(results) <- c("two.round.diff", "two.round.diff.lower", "two.round.diff.upper")
  return(results)
}

#------------------------------------------
# Split data into 2 parts - one for extended genotyping strategies and one for other strategies 
#------------------------------------------
# - use intervention for BMD and Normal (use control for extended genotypes normal cytology)

ref.data1 <- ref.data.cens[ref.data.cens$grp =='i' & ref.data.cens$baseline.cyt != "HSIL" ,]

ref.data2 <- ref.data.cens[(ref.data.cens$baseline.cyt == "ASC-US/LSIL" & ref.data.cens$grp=='i') | 
                             (ref.data.cens$baseline.cyt =="NILM" & ref.data.cens$hpv.genotype.pos & ref.data.cens$grp =='i') |
                             (ref.data.cens$baseline.cyt =="NILM" & ! ref.data.cens$hpv.genotype.pos & ref.data.cens$grp =='c'),]

#old.ref.data1 <- ppv.data1
#old.ref.data2 <- ppv.data2
#------------------------------------------
# apply functions to the data 
#------------------------------------------
ascus.strategy <- c(rep("all", 5), rep("extended", 4), rep("partial", 3), rep("repeat", 2))
nilm.strategy <- c("all", "extended", "partial", "repeat", "extended.repeat", "extended", "partial", "repeat", "extended.repeat", 
                   "partial", "repeat", "extended.repeat", "repeat", "extended.repeat")
ascus.ref <- strategy.function(ref.data1[ref.data1$baseline.cyt=='ASC-US/LSIL',], dim(ref.data1)[1] , "all")
nilm.ref <- strategy.function(ref.data1[ref.data1$baseline.cyt=='NILM',], dim(ref.data1)[1] , "repeat")

results.referral <- list()
for (i in 1:14){
  ascus <- ascus.strategy[i]
  nilm <- nilm.strategy[i]
  if (str_detect(nilm, "extended.repeat")){
    dat <- ref.data2
    }else {
    dat <- ref.data1
  }
  # use control arm for 7-types negative normals for referral rates and PPV calculations 
  results.referral[[i]] <- referral.function(dat, ascus, nilm)
}


#diff.conservative(ref.data2, "partial", "extended.repeat", ref.data2, "repeat", "extended.repeat")
#unlist(diff.conservative(ref.data1, "partial", "partial", ref.data2, "partial", "extended.repeat"))*100

#------------------------------------------
# referral rates - organise the data
#------------------------------------------
{
  referral.results <- round(data.frame(matrix(unlist(results.referral), nrow=length(results.referral), byrow=T)),3) 
  colnames(referral.results) <- names(results.referral[[1]])
  
  referral.results$ascus.strategy <- factor(ascus.strategy)
  referral.results$nilm.strategy <- factor(nilm.strategy, levels=c("all", "extended", "partial", "repeat", "extended.repeat"))
  
  referral.results <- data.frame(rbind(referral.results[1:4,], referral.results[4:length(ascus.strategy),] ))
  referral.results$grp <- NA
  referral.results$grp <- ifelse(referral.results$first.round.rate > referral.results$first.round.rate[4], "More aggressive strategies", "More conservative strategies") 
  referral.results$grp[4] <- "More aggressive strategies"
  referral.results$N.data <- ifelse(str_detect(referral.results$nilm.strategy, "extended"), dim(ref.data2)[1], dim(ref.data1)[1])
  results.reference <- referral.results[4,]
  #referral.results <- referral.results[-5,]
}
# Table 2: referral rates 
referral.results[-5,][order(referral.results$first.round.rate[-5]),] %>%  mutate_if(., is.numeric, function(col) round(100*col, 2)) %>% 
  mutate(baseline.rate = paste0(baseline.rate, " (", baseline.rate.lower, " - ", baseline.rate.upper, ")"),
         repeat.rate = paste0(repeat.rate, " (", repeat.rate.lower, " - ", repeat.rate.upper, ")"),
         second.round.rate = paste0(second.round.rate, " (", second.round.rate.lower, " - ", second.round.rate.upper, ")"),
         overall.hpv.rate = paste0(overall.hpv.rate, " (", overall.hpv.rate.lower, " - ", overall.hpv.rate.upper, ")"),
         first.round.rate = paste0(first.round.rate, " (", first.round.rate.lower, " - ", first.round.rate.upper, ")"), 
         baseline.diff = paste0(baseline.diff, " (", baseline.diff.lower, " - ", baseline.diff.upper, ")"),
         first.round.diff = paste0(first.round.diff, " (", first.round.diff.lower, " - ", first.round.diff.upper, ")"),
         second.round.diff = paste0(two.round.diff, " (", two.round.diff.lower, " - ", two.round.diff.upper, ")")) %>% 
  select(c("baseline.rate", "repeat.rate", "second.round.rate", "overall.hpv.rate", "first.round.rate", "baseline.diff", "first.round.diff", "second.round.diff")) %>% 
  cbind(referral.results[-5,][order(referral.results$first.round.rate[-5]),c("ascus.strategy", "nilm.strategy")], .) %>% 
  write.csv(., "Results/refRates_Oct2023.csv")

#------------------------------------------
# calculate cumulative referral 
#------------------------------------------
{
  cumulative.ref <- data.frame(rate = unlist((data.frame(t(referral.results[,c(1, 13, 10)])))))
  cumulative.ref$lower <- data.frame(rate = unlist((data.frame(t(referral.results[,c(2, 14, 11)])))))
  cumulative.ref$upper <- data.frame(rate = unlist((data.frame(t(referral.results[,c(3, 15, 12)])))))
  
  cumulative.ref$ascus <- factor(rep(referral.results$ascus.strategy, each=3))
  cumulative.ref$nilm <-  factor(rep(referral.results$nilm.strategy, each=3))
  cumulative.ref$grp <-  rep(referral.results$grp, each=3)
  cumulative.ref$strategy <- rep(paste(referral.results$ascus.strategy, referral.results$nilm.strategy), each=3)
  cumulative.ref$time <- rep(c(0.5, 1, 1.5), 15)
  cumulative.ref$code <- rep(c("A1", "A2", "A3", "R", "R", "C1", "A4", "A5", "C2", "C3", "A6", "C4", "C5", "C6", "C7"), each=3)
}

ggplot_refrate <-function(){
  ggplot(cumulative.ref, aes(x=time, y=rate*100, col=ascus, group=strategy, label=code)) + 
  geom_point(aes(shape=nilm,  col=ascus), size=2) +
  geom_line(size=1, aes(col=ascus)) + 
  #geom_ribbon(aes(x=time, ymin=min*100,ymax=max*100, fill=grp), alpha=0.1) + 
  theme_classic() +  
  scale_shape_manual(values=c(9, 15, 17, 16, 8), labels=c("Refer all", "Refer 7 types", "Refer 16/18", "Repeat all", "Repeat 7 types")) + 
  scale_x_continuous(breaks=c(0.5, 1, 1.5), labels=c("Baseline", "First \nrepeat", "5 years")) +
  scale_y_continuous(breaks=seq(0, 100, 10), limits=c(-5, 110), expand=c(0, 0)) +
  scale_color_manual(values = c("#F8766D", "#7CAE00", "#00BFC4","#C77CFF"), labels=c("Refer all", "Refer 7 types", "Refer 16/18", "Repeat all")) +
  theme(legend.position = c(0.8, 0.85), legend.direction='vertical', legend.box='horizontal', 
        text = element_text(size=15), panel.spacing = unit(1, "lines"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), legend.title=element_text(size=10), 
        legend.text=element_text(size=9),
        strip.background = element_blank(), strip.text.x = element_blank()) + 
  guides(col=guide_legend("BMD cytology" ,override.aes = list(shape = NA)), shape=guide_legend("Normal cytology")) + 
  ylab("Cumulative referral rate (%)") + xlab("Screening time point")  +
  facet_wrap("grp") +
  geom_text(data = cumulative.ref[cumulative.ref$time=='1' & cumulative.ref$grp=='More conservative strategies',][1,], 
            aes(x=time, label=grp), col='#505050', y= 50) + 
  geom_text(data = cumulative.ref[cumulative.ref$time=='1' & cumulative.ref$grp=='More aggressive strategies',][1,], 
            aes(x=time, label=grp), col='#505050', y= 105) + 
  geom_text(data=cumulative.ref[cumulative.ref$time=='1.5' & cumulative.ref$grp=='More aggressive strategies',],
            aes(x=time+0.02, y=rate*100, label=code),  col='#505050', hjust = 0, vjust = 0.5, show.legend = F)+
  geom_text(data=cumulative.ref[cumulative.ref$time=='0.5' & cumulative.ref$grp=='More conservative strategies',],
            aes(x=time-0.05, y=rate*100 + rep(c(1.5, -2), 4), label=code), col='#505050', show.legend = F)#,  vjust = 0.5)
}


ggsave(
  "~/Desktop/PhD/Projects/Colposcopy Referral Rates/github/plots/Figure_2.png",
  ggplot_refrate(),
  width = 12,
  height = 6.5,
  dpi = 1000
)
