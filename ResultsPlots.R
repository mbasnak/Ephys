# Code to make plots with the whole-cell data

library(readxl)

pulsePattern <- c(seq(-180,200, 20) , seq(250, 550, 50))
currents <- rep(pulsePattern,28)


#Load the data###########################################

# voltages <- read_excel("Z:/MICROSCOPE/Melanie/ephys/Mel/experiment/everyCell/Voltages.xls", 
#                        col_names = FALSE)

voltages <- read_excel("D:/Doctorado/rotations/Bernardo/Ephys/everyCell/Voltages.xls", 
                       col_names = FALSE)

# status <- read_excel("Z:/MICROSCOPE/Melanie/ephys/Mel/experiment/everyCell/groups.xls", 
#                        col_names = FALSE)

status <- read_excel("D:/Doctorado/rotations/Bernardo/Ephys/everyCell/groups.xls", 
                     col_names = FALSE)

regroupedAPnum <- read_excel("D:/Doctorado/rotations/Bernardo/Ephys/everyCell/regroupedAPnum.xls", 
                     col_names = FALSE)

regroupedInstFR <- read_excel("D:/Doctorado/rotations/Bernardo/Ephys/everyCell/regroupedInstFR.xls", 
                             col_names = FALSE)

regroupedtotFR <- read_excel("D:/Doctorado/rotations/Bernardo/Ephys/everyCell/regroupedtotFR.xls", 
                             col_names = FALSE)


cellName = sort(rep(1:28,27))

Data = cbind(currents,voltages,status,cellName,regroupedAPnum,regroupedInstFR,regroupedtotFR)

colnames(Data) <- c('currents','voltages','status','cellName','APnum','InstFR','totFR')


library(ggplot2)

#IV curves######################################################################33

ggplot(Data, aes(voltages, currents, group=cellName, colour = status)) + 
  geom_point() + geom_line()

library(dplyr)     

summaryStatsVolt = Data %>% group_by(status, currents) %>% summarise(n = n(), 
                                                                 mean_volt = mean(voltages),
                                                                 sem = sd(voltages, na.rm=T)/sqrt(n),
                                                                 ymin = mean_volt-sem,
                                                                 ymax = mean_volt+sem)

#mean IV curves
ggplot(summaryStatsVolt, aes(currents, mean_volt,colour=status)) + geom_line(aes(y=mean_volt), lwd=2)+
  geom_ribbon(aes(ymin=ymin, ymax=ymax,fill=status), alpha=0.4)


#############IF curves##############################

library(cowplot)

plotAPnum <- ggplot(Data, aes(currents, APnum,  group=cellName, colour = status)) + 
  geom_point() + geom_line() + theme(legend.position="none")

plotInstFR <- ggplot(Data, aes(currents, InstFR,  group=cellName, colour = status)) + 
  geom_point() + geom_line() + theme(legend.position="none")

plottotFR <- ggplot(Data, aes(currents, totFR,  group=cellName, colour = status)) + 
  geom_point() + geom_line() + theme(legend.position="none")

plot_grid(plotAPnum, plotInstFR, plottotFR, labels = c("A", "B", "C"),ncol = 3, align = "v")


summaryStatsAPnum = Data %>% group_by(status, currents) %>% summarise(n = n(), 
                                                                     mean_APnum = mean(APnum),
                                                                     sem = sd(APnum, na.rm=T)/sqrt(n),
                                                                     ymin = mean_APnum-sem,
                                                                     ymax = mean_APnum+sem)

summaryStatsInstFR = Data %>% group_by(status, currents) %>% summarise(n = n(), 
                                                                      mean_InstFR = mean(InstFR),
                                                                      sem = sd(InstFR, na.rm=T)/sqrt(n),
                                                                      ymin = mean_InstFR-sem,
                                                                      ymax = mean_InstFR+sem)

summaryStatstotFR = Data %>% group_by(status, currents) %>% summarise(n = n(), 
                                                                      mean_totFR = mean(totFR),
                                                                      sem = sd(totFR, na.rm=T)/sqrt(n),
                                                                      ymin = mean_totFR-sem,
                                                                      ymax = mean_totFR+sem)

plot.meanAPnum <- ggplot(summaryStatsAPnum, aes(currents, mean_APnum,colour=status)) + geom_line(aes(y=mean_APnum), lwd=2)+
  geom_ribbon(aes(ymin=ymin, ymax=ymax,fill=status), alpha=0.4) + theme(legend.position='none')

plot.meanInstFR <- ggplot(summaryStatsInstFR, aes(currents, mean_InstFR,colour=status)) + geom_line(aes(y=mean_InstFR), lwd=2)+
  geom_ribbon(aes(ymin=ymin, ymax=ymax,fill=status), alpha=0.4) + theme(legend.position='none')

plot.meantotFR <- ggplot(summaryStatstotFR, aes(currents, mean_totFR,colour=status)) + geom_line(aes(y=mean_totFR), lwd=2)+
  geom_ribbon(aes(ymin=ymin, ymax=ymax,fill=status), alpha=0.4) + theme(legend.position='none')


plot_grid(plot.meanAPnum, plot.meanInstFR, plot.meantotFR, labels = c("A", "B", "C"),ncol = 3, align = "v")


####################### resting properties ########################

Vrest <- read_excel("D:/Doctorado/rotations/Bernardo/Ephys/everyCell/Vrest.xls", 
                       col_names = FALSE)
Ihold <- read_excel("D:/Doctorado/rotations/Bernardo/Ephys/everyCell/Ihold.xls", 
                       col_names = FALSE)
cellType <- read_excel("D:/Doctorado/rotations/Bernardo/Ephys/everyCell/cellType.xls", 
                       col_names = FALSE)
Vrest <- t(Vrest)
Ihold <- t(Ihold)

restingProp <- cbind(Vrest,Ihold,cellType)
colnames(restingProp) <- c('Vrest','Ihold','cellType')

plot.Vrest <- ggplot(restingProp, aes(cellType,Vrest,fill=cellType)) + geom_boxplot() + geom_jitter() + theme(legend.position='none')

plot.Ihold <- ggplot(restingProp, aes(cellType,Ihold,fill=cellType)) + geom_boxplot() + geom_jitter() + theme(legend.position='none')

plot_grid(plot.Vrest, plot.Ihold, labels = c("A", "B"),ncol = 2, align = "v")


##################### AP properties ############################

APproperties <- read_excel("D:/Doctorado/rotations/Bernardo/Ephys/everyCell/APproperties.xls", 
                       col_names = FALSE)
APproperties <- t(APproperties)

APprop <- cbind(APproperties, cellType)

colnames(APprop) <- c('maxInstFR','APthreshold','latency','APamplitude','APhalfwidth','maxtotFR','APthrough','cellType')

plot.maxInstFR <- ggplot(APprop, aes(cellType,maxInstFR,fill=cellType)) + geom_boxplot() + geom_jitter() + theme(legend.position='none')
plot.APthreshold <- ggplot(APprop, aes(cellType,APthreshold,fill=cellType)) + geom_boxplot() + geom_jitter() + theme(legend.position='none')
plot.latency <- ggplot(APprop, aes(cellType,latency,fill=cellType)) + geom_boxplot() + geom_jitter() + theme(legend.position='none')
plot.APamplitude <- ggplot(APprop, aes(cellType,APamplitude,fill=cellType)) + geom_boxplot() + geom_jitter() + theme(legend.position='none')
plot.APhalfwidth <- ggplot(APprop, aes(cellType,APhalfwidth,fill=cellType)) + geom_boxplot() + geom_jitter() + theme(legend.position='none')
plot.maxtotFR <- ggplot(APprop, aes(cellType,maxtotFR,fill=cellType)) + geom_boxplot() + geom_jitter() + theme(legend.position='none')

plot_grid(plot.maxInstFR, plot.APthreshold, plot.latency, plot.APamplitude, plot.APhalfwidth, plot.maxtotFR, labels = c("A", "B", "C", "D", "E", "F"),ncol = 3)


######################## PCA #################################

Resistance <- read_excel("D:/Doctorado/rotations/Bernardo/Ephys/everyCell/Resistance.xls", 
                           col_names = FALSE)
Resistance <- as.data.frame(t(Resistance))
colnames(Resistance) <- c('Rin','Cm')

PCAvariables = cbind(APprop$maxtotFR,APprop$APamplitude,APprop$APthreshold,APprop$latency,restingProp$Vrest,Resistance$Rin,Resistance$Cm);
PCAvariables <- as.data.frame(PCAvariables)
colnames(PCAvariables) <- c('maxtotFR','APamplitude','APthreshold','latency','Vrest','Rin','Cm')

dominance.pca <- prcomp(PCAvariables,
                 center = TRUE,
                 scale. = TRUE)
print(dominance.pca)
plot(dominance.pca, type = "l")

library(devtools)
install_github("ggbiplot", "vqv")

library(ggbiplot)
g <- ggbiplot(dominance.pca, obs.scale = 1, var.scale = 1, 
              groups = APprop$cellType, ellipse = TRUE )
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
print(g)
