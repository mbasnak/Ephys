# Code to make plots with the whole-cell data

library(readxl)

# Generate a vector with the current pattern
pulsePattern <- c(seq(-180,200, 20) , seq(250, 550, 50))
currents <- rep(pulsePattern,28)

#Load the data###########################################

# voltages <- read_excel("Z:/MICROSCOPE/Melanie/ephys/Mel/experiment/everyCell/Voltages.xls", 
#                        col_names = FALSE)

voltages <- read_excel("D:/Doctorado/rotations/Bernardo/Ephys/everyCell/Voltages.xls", 
                       col_names = FALSE)

# status <- read_excel("Z:/MICROSCOPE/Melanie/ephys/Mel/experiment/everyCell/groups.xls", 
#                        col_names = FALSE)

Status <- read_excel("D:/Doctorado/rotations/Bernardo/Ephys/everyCell/groups.xls", 
                     col_names = FALSE)

regroupedAPnum <- read_excel("D:/Doctorado/rotations/Bernardo/Ephys/everyCell/regroupedAPnum.xls", 
                     col_names = FALSE)

regroupedInstFR <- read_excel("D:/Doctorado/rotations/Bernardo/Ephys/everyCell/regroupedInstFR.xls", 
                             col_names = FALSE)

regroupedtotFR <- read_excel("D:/Doctorado/rotations/Bernardo/Ephys/everyCell/regroupedtotFR.xls", 
                             col_names = FALSE)


cellName = sort(rep(1:28,27))

# Save all of the data into a data.frame
Data = cbind(currents,voltages,Status,cellName,regroupedAPnum,regroupedInstFR,regroupedtotFR)
colnames(Data) <- c('currents','voltages','Status','cellName','APnum','InstFR','totFR')


library(ggplot2)

#IV curves######################################################################

### individual curves, colored by dominance status ###

ggplot(Data, aes(voltages, currents, group=cellName, colour = Status)) +
  geom_point() + geom_line() + 
  scale_fill_brewer(palette='Set1') +
  labs(x = "Volatge (mV)", y = "Current (pA)",title='Individual IV curves')+ 
  theme(plot.title = element_text(hjust = 0.5))+ theme(legend.title = element_text(face="bold"))

ggplot(Data, aes(currents, voltages, group=cellName, colour = Status)) +
  geom_point() + geom_line() + 
  scale_fill_brewer(palette='Set1') +
  labs(y = "Volatge (mV)", x = "Current (pA)",title='Individual IV curves')+ 
  theme(plot.title = element_text(hjust = 0.5)) + theme(legend.title = element_text(face="bold"))+
  geom_vline(xintercept=0)



library(dplyr)     

# Save the mean and sem of the voltage
summaryStatsVolt = Data %>% group_by(Status, currents) %>% summarise(n = n(), 
                                                                 mean_volt = mean(voltages),
                                                                 sem = sd(voltages, na.rm=T)/sqrt(n),
                                                                 ymin = mean_volt-sem,
                                                                 ymax = mean_volt+sem)

# mean IV curves
ggplot(summaryStatsVolt, aes(currents, mean_volt,colour=Status)) + geom_line(aes(y=mean_volt), lwd=2)+
  geom_ribbon(aes(ymin=ymin, ymax=ymax,fill=Status), alpha=0.4) +
  scale_fill_brewer(palette='Set1') +
  labs(y = "Volatge (mV)", x = "Current (pA)",title='Mean IV curves')+ 
  theme(plot.title = element_text(hjust = 0.5)) + theme(legend.title = element_text(face="bold"))+
  geom_vline(xintercept=0)


############# IF curves ##############################

library(cowplot)

plotAPnum <- ggplot(Data, aes(currents, APnum,  group=cellName, colour = Status)) + 
  geom_point() + geom_line() + theme(legend.position="none") +
  scale_fill_brewer(palette='Set1') + labs(x= 'Current (pA)',y='AP number')

plotInstFR <- ggplot(Data, aes(currents, InstFR,  group=cellName, colour = Status)) + 
  geom_point() + geom_line() + theme(legend.position="none") +
  scale_fill_brewer(palette='Set1') + labs(x= 'Current (pA)',y='Instantaneous firing rate (Hz)')

plottotFR <- ggplot(Data, aes(currents, totFR,  group=cellName, colour = Status)) + 
  geom_point() + geom_line() + theme(legend.position="none") +
  scale_fill_brewer(palette='Set1') + labs(x= 'Current (pA)',y='Total firing rate (Hz)')

plot_grid(plotAPnum, plotInstFR, plottotFR, labels = c("A", "B", "C"),ncol = 3, align = "v")

# Save the summary stats of the different firing rates
summaryStatsAPnum = Data %>% group_by(Status, currents) %>% summarise(n = n(), 
                                                                     mean_APnum = mean(APnum),
                                                                     sem = sd(APnum, na.rm=T)/sqrt(n),
                                                                     ymin = mean_APnum-sem,
                                                                     ymax = mean_APnum+sem)

summaryStatsInstFR = Data %>% group_by(Status, currents) %>% summarise(n = n(), 
                                                                      mean_InstFR = mean(InstFR),
                                                                      sem = sd(InstFR, na.rm=T)/sqrt(n),
                                                                      ymin = mean_InstFR-sem,
                                                                      ymax = mean_InstFR+sem)

summaryStatstotFR = Data %>% group_by(Status, currents) %>% summarise(n = n(), 
                                                                      mean_totFR = mean(totFR),
                                                                      sem = sd(totFR, na.rm=T)/sqrt(n),
                                                                      ymin = mean_totFR-sem,
                                                                      ymax = mean_totFR+sem)

plot.meanAPnum <- ggplot(summaryStatsAPnum, aes(currents, mean_APnum,colour=Status)) + geom_line(aes(y=mean_APnum), lwd=2)+
  geom_ribbon(aes(ymin=ymin, ymax=ymax,fill=Status), alpha=0.4) +
  theme(legend.position='none') + scale_fill_brewer(palette='Set1') +
  labs(x='Current (pA)', y='AP number')

plot.meanInstFR <- ggplot(summaryStatsInstFR, aes(currents, mean_InstFR,colour=Status)) + geom_line(aes(y=mean_InstFR), lwd=2)+
  geom_ribbon(aes(ymin=ymin, ymax=ymax,fill=Status), alpha=0.4) +
  theme(legend.position='none') + scale_fill_brewer(palette='Set1') +
  labs(x='Current (pA)', y='Instantaneous firing rate (Hz)')

plot.meantotFR <- ggplot(summaryStatstotFR, aes(currents, mean_totFR,colour=Status)) + geom_line(aes(y=mean_totFR), lwd=2)+
  geom_ribbon(aes(ymin=ymin, ymax=ymax,fill=Status), alpha=0.4) +
  theme(legend.position='none') + scale_fill_brewer(palette='Set1') +
  labs(x='Current (pA)', y='Total firing rate (Hz)')


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

plot.Vrest <- ggplot(restingProp, aes(cellType,Vrest,fill=cellType)) + 
  geom_boxplot() + scale_fill_brewer(palette='Set1') +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.85, fill='gray40') +
  theme(legend.position='none') +
  labs(x='',y='Resting membrane potential (mV)')

plot.Ihold <- ggplot(restingProp, aes(cellType,Ihold,fill=cellType)) + 
  geom_boxplot() + scale_fill_brewer(palette='Set1') +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.85, fill='gray40') +
  theme(legend.position='none') +
  labs(x='',y='Holding current (pA)')

plot_grid(plot.Vrest, plot.Ihold, labels = c("A", "B"),ncol = 2, align = "v")


##################### AP properties ############################

APproperties <- read_excel("D:/Doctorado/rotations/Bernardo/Ephys/everyCell/APproperties.xls", 
                       col_names = FALSE)
APproperties <- t(APproperties)

APprop <- cbind(APproperties, cellType)

colnames(APprop) <- c('maxInstFR','APthreshold','latency','APamplitude','APhalfwidth','maxtotFR','APthrough','sag','cellType')

plot.maxInstFR <- ggplot(APprop, aes(cellType,maxInstFR,fill=cellType)) +
  geom_boxplot() + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.85, fill='gray40') +
  theme(legend.position='none') + scale_fill_brewer(palette='Set1') +
  labs(x='',y='Max. instantaneous FR (Hz)')

plot.APthreshold <- ggplot(APprop, aes(cellType,APthreshold,fill=cellType)) +
  geom_boxplot() + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.85, fill='gray40') +
  theme(legend.position='none') + scale_fill_brewer(palette='Set1') +
  labs(x='',y='AP threshold (mV)')

plot.latency <- ggplot(APprop, aes(cellType,latency,fill=cellType)) +
  geom_boxplot() + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.85, fill='gray40') +
  theme(legend.position='none') + scale_fill_brewer(palette='Set1') +
  labs(x='',y='Latency to first AP (ms)')

plot.APamplitude <- ggplot(APprop, aes(cellType,APamplitude,fill=cellType)) +
  geom_boxplot() + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.85, fill='gray40') +
  theme(legend.position='none') + scale_fill_brewer(palette='Set1') +
  labs(x='',y='AP amplitude (mV)')

plot.APhalfwidth <- ggplot(APprop, aes(cellType,APhalfwidth,fill=cellType)) +
  geom_boxplot() + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.85, fill='gray40') +
  theme(legend.position='none') + scale_fill_brewer(palette='Set1') +
  labs(x='',y='AP halfwidth (ms)')

plot.maxtotFR <- ggplot(APprop, aes(cellType,maxtotFR,fill=cellType)) +
  geom_boxplot() + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.85, fill='gray40') +
  theme(legend.position='none') + scale_fill_brewer(palette='Set1') +
  labs(x='',y='Max. total FR (Hz)')

plot.APthrough <- ggplot(APprop, aes(cellType,APthrough,fill=cellType)) +
  geom_boxplot() + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.85, fill='gray40') +
  theme(legend.position='none') + scale_fill_brewer(palette='Set1') +
  labs(x='',y='AP through (mV)')

plot.sag <- ggplot(APprop, aes(cellType,sag,fill=cellType)) +
  geom_boxplot() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.85, fill='gray40') +
  theme(legend.position='none') + scale_fill_brewer(palette='Set1') +
  labs(x='',y='Sag ratio for first current')

plot_grid(plot.maxInstFR, plot.APthreshold, plot.latency, plot.APamplitude, 
          plot.APhalfwidth, plot.maxtotFR, plot.APthrough, plot.sag,
          labels = c("A", "B", "C", "D", "E", "F","G","H"),ncol = 4)


######################## PCA #################################

Resistance <- read_excel("D:/Doctorado/rotations/Bernardo/Ephys/everyCell/Resistance.xls", 
                           col_names = FALSE)
Resistance <- as.data.frame(t(Resistance))
colnames(Resistance) <- c('Rin','Cm')

PCAvariables = cbind(APprop$maxtotFR,APprop$APamplitude,APprop$APthreshold,APprop$latency,APprop$sag,restingProp$Vrest,Resistance$Rin,Resistance$Cm);
PCAvariables <- as.data.frame(PCAvariables)
colnames(PCAvariables) <- c('maxtotFR','APamplitude','APthreshold','latency','sag','Vrest','Rin','Cm')

dominance.pca <- prcomp(PCAvariables,
                 center = TRUE,
                 scale. = TRUE)
print(dominance.pca)

library(devtools)
library(ggbiplot)

g <- ggbiplot(dominance.pca, obs.scale = 1, var.scale = 1, 
              groups = APprop$cellType, ellipse = TRUE, varname.size = 5)
g <- g  + scale_colour_manual(values = c("red", "blue")) 
print(g)

