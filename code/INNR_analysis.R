# INNR analysis - Figure 3
# 
library(ggthemes)
library(ggplot2)
library(ggrepel)
library(ggtext)
# library(dampack)
innr.function <- function(data, ascus.strategy, nilm.strategy){
  total <- dim(data)[1]
  
  ascus.results <- ppv.cases.function(data[data$baseline.cyt=="ASC-US/LSIL", ], ascus.strategy)
  nilm.results <- ppv.cases.function(data[data$baseline.cyt=="NILM", ], nilm.strategy)
  referred <- data$id[data$id %in% c(ascus.results$first.round.referred, nilm.results$first.round.referred)]
    
  #detected.cin2plus <- sum(data[data$id %in% referred, ]$first.round.cin2plus)
  detected.cin3plus <- sum(data[data$id %in% referred, ]$first.round.cin3plus)
  return(list(length(referred), detected.cin3plus))
}


innr.data.temp <- ref.data.cens[ref.data.cens$baseline.cyt == "ASC-US/LSIL" | (ref.data.cens$grp =='i' & ref.data.cens$baseline.cyt != "HSIL" ),]
innr.double.data <- data.frame(rbind(innr.data.temp, innr.data.temp[innr.data.temp$baseline.cyt=='NILM',]))

results.innr <- list()
for (i in 1:14){
  ascus <- ascus.strategy[i]
  nilm <- nilm.strategy[i]
  results.innr[[i]] <- innr.function(innr.double.data, ascus, nilm)
}

#------------------------------------------
# INNR plot 
#------------------------------------------
# Figure 2: effiency frontier plot 
{
innr.temp <- round(data.frame(matrix(unlist(results.innr), nrow=length(results.innr), byrow = T)), 3)
colnames(innr.temp) <- c("referred", "detected")
innr.temp$ascus.strategy <- factor(ascus.strategy)
innr.temp$nilm.strategy <- factor(nilm.strategy, levels=c("all", "extended", "partial", "repeat", "extended.repeat"))

innr.temp
strategies <- paste(c(rep("refer all", 5), rep("refer extended", 4), rep("refer 16/18", 3), rep("all repeat", 2)), 
                    c("refer all", "refer extended", "refer 16/18", "all repeat", "repeat extended",
                      "refer extended", "refer 16/18", "all repeat", "repeat extended", 
                      "refer 16/18", "all repeat", "repeat extended", "all repeat", "repeat extended"), sep=' - ')
code <- c("A1", "A2", "A3", "R", "C1", "A4", "A5", "C2", "C3", "A6", "C4", "C5", "C6", "C7", "nothing")

innr.data <- data.frame(strategy = c(strategies, "nothing"),
                        referred = c(innr.temp$referred, 0),
                        detected = c(innr.temp$detected,0), 
                        code=code)

innr.data$strategy.long <- c("A1\nBMD: Refer all \nNormal: Refer all", "A2\nBMD: Refer all \nNormal: Refer 7 types", "A3\nBMD: Refer all \nNormal: Refer 16/18",
                             "R\nBMD: Refer all \nNormal: Repeat all", "C1\nBMD: Refer all \nNormal: Repeat 7 types", "A4\nBMD: Refer 7 types \nNormal: Refer 7 types",
                             "A5\nBMD: Refer 7 types \nNormal: Refer 16/18", "C2\nBMD: Refer 7 types \nNormal: Repeat all", "C3\nBMD: Refer 7 types \nNormal: Repeat 7 types", 
                             "A6\nBMD: Refer 16/18 \nNormal: Refer 16/18", "C4\nBMD: Refer 16/18 \nNormal: Repeat all", "C5\nBMD: Refer 16/18 \nNormal: Repeat 7 types", 
                             "C6\nBMD: Repeat all \nNormal: Repeat all", "C7\nBMD: Repeat all \nNormal: Repeat 7 types", "No screening")

innr.data <- innr.data[rev(order(innr.data$referred)),]

innr.results <- calculate_icers(cost=innr.data$referred, 
                                effect=innr.data$detected,
                                strategies=innr.data$strategy)
innr.results <- innr.results[rev(order(innr.results$Cost)),]
innr.data$Status <- innr.results$Status

status_expand <- c("D" = "Dominated", "ED" = "Extended Dominated", "ND" = "Efficient Frontier", "ref" = "ref")
plot_lines <- c("Dominated" = "blank", "Extended Dominated" = "blank", "Efficient Frontier" = "solid")
col_lines <- c("Dominated" = "#ababab", "Extended Dominated" = "#ababab", "Efficient Frontier" = "darkred")

innr.results$Status <- factor(status_expand[innr.results$Status], ordered = FALSE,
                              levels = c("Dominated", "Extended Dominated", "Efficient Frontier"))
innr.data$status <- innr.results$Status
innr.data$innr <- innr.results$ICER
innr.data$innr.label <- paste0("mPPV: " , round(100/innr.data$innr, 1), "%")
innr.plot.dat <- innr.data[-15,]
innr.plot.dat$strategy.long.nnr <- sapply(strsplit(paste0(innr.plot.dat$strategy.long, "\n", innr.plot.dat$innr.label), '\n'), 
                                          function(x) paste0("atop(", "atop('", x[1], "','", x[2], "'), atop('", x[3], "',bold('", x[4], "')))"))
# innr.plot.dat$strategy.long.nnr <- ifelse(innr.plot.dat$Status!='ND', sapply(strsplit(paste0(innr.plot.dat$strategy.long, "\n", innr.plot.dat$innr.label), '\n'), 
#                                                                          function(x) str_c(x[1], "\n", x[2], "\n", x[3]) ),  innr.plot.dat$strategy.long.nnr)
#                                                                            #paste0("atop(", "atop('", x[1], "','", x[2], "'), atop('", x[3], "'))")) ,) 

innr.plot.dat$strategy.long.nnr <- ifelse(innr.plot.dat$Status!='ND', innr.plot.dat$strategy.long, innr.plot.dat$strategy.long.nnr)

}
innr.plot.dat
# Figure 2:  INNR (first round referred and detected)
par(lheight=0.2)  
ggplot_innr <- function() {
  ggplot(innr.plot.dat, aes(x = referred, y = detected, label=strategy.long.nnr, shape = status, linetype=status, color=status, group = status)) +
  geom_point(size = 2) +
  geom_line(lwd=1) +
  scale_colour_manual("", values=col_lines) +
  scale_linetype_manual("", values=plot_lines) +
  scale_shape_manual("", values=c(19,17, 15)) + 
  labs(y = paste0("Number of CIN3+ cases detected"), x = paste0("Number of colposcopy referrals")) + theme_classic() +
  theme(legend.position = c(0.87, 0.13), text = element_text(size=14), axis.text=element_text(colour="black")) +
  geom_text_repel(data=subset(innr.plot.dat, Status=='ND'), parse = T,
                  arrow = arrow(length = unit(0.015, "npc")), point.padding=unit(0.8,'lines'), 
                  box.padding=1, col='black',segment.color = 'black', size=5, show.legend = FALSE, nudge_x=-130, nudge_y =8, segment.size=0.2) +
  geom_text(data=subset(innr.plot.dat, Status!='ND'),  hjust=0, 
                  aes(label=strategy.long, x=referred - c(5, 5, 5, 5, -15, -16, 5, 5), 
                      y=detected -c(1.5, 1.5, 1.5, 1.5, -0.4, -0.4, 1.5, 1.5)), #point.padding=unit(0.2,'lines'), 
                  col='#ababab', size=2, show.legend = FALSE) + 
  scale_x_continuous(limits=c(120, 1900), breaks=seq(200, 1900, 400)) +
  scale_y_continuous(limits=c(100, 135), breaks=seq(100, 135, 5)) 
}


ggsave(
  "~/Desktop/PhD/Projects/Colposcopy Referral Rates/github/plots/Figure_3.png",
  ggplot_innr(),
  width = 12,
  height = 6.5,
  dpi = 1000
)

