library(CubSFS)
library(parallel)
library(tidyverse)


data <- read_table2("popsize_1mil.txt")

ggplot(data) +
  geom_ribbon(aes(x=time_E,ymin=LCL_E,ymax=UCL_E),alpha=0.2,color="black",size=0.1) +
  theme_bw(base_size = 4) + scale_y_log10(expand=c(0,0)) + ylab("Effective population size") + xlab ("Time (years)") +
  theme(axis.title=element_text(face="bold"), legend.title = element_text(face="bold")) + theme(aspect.ratio = 1) + 
  theme(legend.key.size = unit(0.5,"line")) +
#geom_line(aes(x=time,y=Nt),size=0.5) +
  geom_line(aes(x=time_E,y=median_E),size=0.5,linetype="dashed") + xlab ("Time (years)") +
  theme(axis.title=element_text(face="bold")) +
  geom_ribbon(aes(x=time_W,ymin=LCL_W,ymax=UCL_W,fill='#b5b0b0'),alpha=0.2,color="black",size=0.1) +
  theme_bw(base_size = 4) + scale_y_log10(expand=c(0,0)) + ylab("Effective population size") + xlab ("Time (years)") +
  theme(axis.title=element_text(face="bold"), legend.title = element_text(face="bold")) + theme(aspect.ratio = 1) + 
  theme(legend.key.size = unit(0.5,"line")) +
#geom_line(aes(x=time,y=Nt),size=0.5) +
  geom_line(aes(x=time_W,y=median_W),size=0.5,linetype="dashed",color='Red') + xlab ("Time (years)") +
  theme(axis.title=element_text(face="bold"))



ggsave("Both_pop_1mil.pdf",width=8.7,height=17.8,units="cm")





data <- read_table2("popsize_200.txt")

ggplot(data) +
  geom_ribbon(aes(x=time_E,ymin=LCL_E,ymax=UCL_E),alpha=0.2,color="black",size=0.1) +
  theme_bw(base_size = 4) + scale_y_log10(expand=c(0,0)) + ylab("Effective population size") + xlab ("Time (years)") +
  theme(axis.title=element_text(face="bold"), legend.title = element_text(face="bold")) + theme(aspect.ratio = 1) + 
  theme(legend.key.size = unit(0.5,"line")) +
#geom_line(aes(x=time,y=Nt),size=0.5) +
  geom_line(aes(x=time_E,y=median_E),size=0.5,linetype="dashed") + xlab ("Time (years)") +
  theme(axis.title=element_text(face="bold")) +
  geom_ribbon(aes(x=time_W,ymin=LCL_W,ymax=UCL_W,fill='#b5b0b0'),alpha=0.2,color="black",size=0.1) +
  theme_bw(base_size = 4) + scale_y_log10(expand=c(0,0)) + ylab("Effective population size") + xlab ("Time (years)") +
  theme(axis.title=element_text(face="bold"), legend.title = element_text(face="bold")) + theme(aspect.ratio = 1) + 
  theme(legend.key.size = unit(0.5,"line")) +
#geom_line(aes(x=time,y=Nt),size=0.5) +
  geom_line(aes(x=time_W,y=median_W),size=0.5,linetype="dashed",color='Red') + xlab ("Time (years)") +
  theme(axis.title=element_text(face="bold"))



ggsave("Both_pop_200.pdf",width=8.7,height=17.8,units="cm")



