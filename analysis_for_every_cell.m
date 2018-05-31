% Analysis for all the cells together

close all; clear all;

% choose directory of cell to be analyzed
CellPath = uigetdir('\\files.med.harvard.edu\Neurobio\MICROSCOPE\Melanie\ephys\Mel','Choose folder')

% set the current path to that directory
path = cd(CellPath);

%% For every folder inside a day, import the struct

% Get list of all subfolders.
allSubFolders = genpath(CellPath);
% Parse into a cell array.
remain = allSubFolders;
listOfFolderNames = {};
while true
	[singleSubFolder, remain] = strtok(remain, ';');
	if isempty(singleSubFolder)
		break;
	end
	listOfFolderNames = [listOfFolderNames singleSubFolder];
end
numberOfFolders = length(listOfFolderNames);

for i= 2:length(listOfFolderNames) 
    cellData{i-1} = load(strcat(listOfFolderNames{i}, '\cellProp.mat'));  
end

% Extract specifically the struct
for i = 1:size(cellData,2)
    Data{i} = cellData{1,i}.cellProp;
end

%% IV curves with every cell

pulseStart=5000; % when in the protocol is our pulse starting (in points)
pulseEnd=15000; % when is it ending

% Extract the current pulses, neuron's voltage response and firing rates
% from the data
for i = 1:length(Data)
    voltage(:,i) = Data{1,i}.voltage;
    current(:,i) = Data{1,i}.currents;
    APnum(:,i) = Data{1,i}.APnum;

end


for i = 1:length(Data)
    for j = 1:size(current,1)
            InstFR(j,i) = Data{1,i}(j).InstFR;
            totFR(j,i) = Data{1,i}(j).totFR;
    end
end

% Plot the IV curve
%figure, set(gcf,'units','points','position',[100,100,1000,600]); %if I run it in lab
figure, set(gcf,'units','points','position',[80,80,600,350]); %if I run it in my laptop
plot(voltage,current,'o')
hold on
plot(voltage,current)
hline = refline([0 0]);
hline.Color = 'k';
title('IV curve');
ylim([-200 600]); xlim([-150 50]);
xlabel('V (mV)');ylabel('Current (pA)');

saveas(gcf,'IVcurvesEveryCell.png');

%% IF curves three ways with every cell

% Plot the IF curve
%figure, set(gcf,'units','points','position',[100,100,1000,600]); %if I run it in lab
figure, set(gcf,'units','points','position',[80,80,600,350]); %if I run it in my laptop
subplot(1,3,1)
plot(current,APnum,'^')
hold on
plot(current,APnum)
title('If curve');
xlim([0 600]);
ylabel('AP number');xlabel('Current (pA)');

subplot(1,3,2)
plot(current,InstFR,'^')
hold on
plot(current,InstFR)
title('If curve');
xlim([0 600]);
ylabel('Instantaneous firing rate (spikes/s)');xlabel('Current (pA)');

subplot(1,3,3)
plot(current,totFR,'^')
hold on
plot(current,totFR)
title('If curve');
xlim([0 600]);
ylabel('Total firing rate (spikes/s)');xlabel('Current (pA)');

%% Plot with cells grouped by dominance

%If we saved the folders as Dcellx and Scellx, all of the Dominant
%subfolders are going to be saved first, and the subordinate second

% 1) Reshape the voltage and current data into vectors
regroupedV = reshape(voltage,[size(voltage,1)*size(voltage,2),1]);
regroupedC = reshape(current,[size(current,1)*size(current,2),1]);

% 2) Generate a char vector with information about the group
    % a) Check to see how many D and S subfolders we have
    names = dir;
    for i = 1:length(names)
        FolderNames{i} = names(i).name;
    end
    
    dominant = strfind(FolderNames,'DCell'); dominantNum = sum(cell2mat(dominant));
    subordinate = strfind(FolderNames,'SCell'); subordinateNum = sum(cell2mat(subordinate));
    
    % b) Make the vector
groups = cell(size(regroupedV,1),1);
groups(1:(size(voltage,1)*dominantNum),1) = {'Dominant'};
groups((size(voltage,1)*dominantNum)+1:end,1) = {'Subordinate'};

% Plot IV curves sorted by dominance
%figure, set(gcf,'units','points','position',[100,100,1000,600]); %if I run it in lab
figure, set(gcf,'units','points','position',[80,80,600,350]); %if I run it in my laptop
gscatter(regroupedV,regroupedC,groups)
hold on
plot(voltage,current)
title('IV curve');
ylim([-200 600]); xlim([-150 50]);
xlabel('V (mV)');ylabel('Current (pA)');

saveas(gcf,'IVcurveByStatus.png');

%% IF curves by status

regroupedAPnum = reshape(APnum,[size(APnum,1)*size(APnum,2),1]);
regroupedInstFR = reshape(InstFR,[size(InstFR,1)*size(InstFR,2),1]);
regroupedtotFR = reshape(totFR,[size(totFR,1)*size(totFR,2),1]);

figure, set(gcf,'units','points','position',[80,80,600,350]); %if I run it in my laptop
subplot(1,3,1)
gscatter(regroupedC,regroupedAPnum,groups)
hold on
plot(current,APnum)
title('If curve');
xlim([0 600]);
ylabel('AP num (mV)');xlabel('Current (pA)');

subplot(1,3,2)
gscatter(regroupedC,regroupedInstFR,groups)
hold on
plot(current,InstFR)
title('If curve');
xlim([0 600]);
ylabel('Instantaneous firing rate (spikes/s)');xlabel('Current (pA)');

subplot(1,3,3)
gscatter(regroupedC,regroupedtotFR,groups)
hold on
plot(current,totFR)
title('If curve');
xlim([0 600]);
ylabel('Total firing rate (spikes/s)');xlabel('Current (pA)');

%% Plot IV and IF mean curves by dominance status

meanCD = mean(current(:,1:dominantNum),2);
meanCS = mean(current(:,dominantNum+1:end),2);

meanVD = mean(voltage(:,1:dominantNum),2);
meanVS = mean(voltage(:,dominantNum+1:end),2);
stdVD = std(voltage(:,1:dominantNum),[],2);
stdVS = std(voltage(:,dominantNum+1:end),[],2);

regroupedMeanC = vertcat(meanCD,meanCS);
regroupedMeanV = vertcat(meanVD,meanVS);

regroupedStdV = vertcat(stdVD,stdVS);

status = cell(size(voltage,1)*2,1);
status(1:length(voltage),1) = {'Dominant'};
status(length(voltage)+1:end,1) = {'Subordinate'};
status = char(status);

%figure, set(gcf,'units','points','position',[100,100,1000,600]); %if I run it in lab
figure, set(gcf,'units','points','position',[80,80,600,350]); %if I run it in my laptop
gscatter(regroupedMeanV,regroupedMeanC,status)
title('IV curve');
ylim([-200 600]); xlim([-150 50]);
xlabel('V (mV)');ylabel('Current (pA)');

%% IF curves

meanFRD = mean(APnum(:,1:dominantNum),2);
meanFRS = mean(APnum(:,dominantNum+1:end),2);
stdFRD = std(APnum(:,1:dominantNum),[],2);
stdFRS = std(APnum(:,dominantNum+1:end),[],2);

regroupedMeanFR = vertcat(meanFRD,meanFRS);
regroupedStdFR = vertcat(stdFRD,stdFRS);

meanInstFRD = mean(InstFR(:,1:dominantNum),2);
meanInstFRS = mean(InstFR(:,dominantNum+1:end),2);
stdInstFRD = std(InstFR(:,1:dominantNum),[],2);
stdInstFRS = std(InstFR(:,dominantNum+1:end),[],2);

regroupedMeanInstFR = vertcat(meanInstFRD,meanInstFRS);
regroupedStdInstFR = vertcat(stdInstFRD,stdInstFRS);

meantotFRD = mean(totFR(:,1:dominantNum),2);
meantotFRS = mean(totFR(:,dominantNum+1:end),2);
stdtotFRD = std(totFR(:,1:dominantNum),[],2);
stdtotFRS = std(totFR(:,dominantNum+1:end),[],2);

regroupedMeantotFR = vertcat(meantotFRD,meantotFRS);
regroupedStdtotFR = vertcat(stdtotFRD,stdtotFRS);

%figure, set(gcf,'units','points','position',[100,100,1000,600]); %if I run it in lab
figure, set(gcf,'units','points','position',[80,80,600,350]); %if I run it in my laptop
subplot(1,3,1)
gscatter(regroupedMeanC,regroupedMeanFR,status)
xlim([0 600]);
ylabel('AP number');xlabel('Current (pA)');

subplot(1,3,2)
gscatter(regroupedMeanC,regroupedMeanInstFR,status)
xlim([0 600]);
ylabel('Instantaneous firing rate (spikes/s)');xlabel('Current (pA)');

subplot(1,3,3)
gscatter(regroupedMeanC,regroupedMeantotFR,status)
xlim([0 600]);
ylabel('Total firing rate (spikes/s)');xlabel('Current (pA)');

%% Compare the resting properties

% Extract the variables from the data
for i = 1:length(Data)
    Vrest(:,i) = Data{1,i}(1).Vrest;
    Ihold(:,i) = Data{1,i}(1).Ihold;
end

% Make boxplots
    % group by cellnumber
cellType = cell(size(Data,2),1);
cellType(1:dominantNum,1) = {'Dominant'};
cellType(dominantNum+1:end,1) = {'Subordinate'};

%figure, set(gcf,'units','points','position',[100,100,1000,600]); %if I run it in lab
figure, set(gcf,'units','points','position',[80,80,600,350]); %if I run it in my laptop
subplot(1,2,1)
boxplot(Vrest,cellType)
ylabel('Voltage (mV)'); title('Vrest values');
subplot(1,2,2)
boxplot(Ihold,cellType)
ylabel('Current (pA)'); title('Iholding current values');

saveas(gcf,'CellPropertiesComparison.png');

%% Compare AP properties

% Extract the variables from the data
for i = 1:length(Data)
    maxInstFR(:,i) = Data{1,i}(1).maxInstFR;
    APthreshold(:,i) = Data{1,i}.APthreshold;
    latency(:,i) = Data{1,i}.Latency;
    APamplitude(:,i) = Data{1,i}.APamplitude;
    APhalfwidth(:,i) = Data{1,i}.APhalfwidth;
    maxtotFR(:,i) = Data{1,i}(1).maxtotFR;
end


% Make boxplots

%figure, set(gcf,'units','points','position',[100,100,1000,600]); %if I run it in lab
figure, set(gcf,'units','points','position',[80,80,600,350]); %if I run it in my laptop
subplot(2,3,1)
boxplot(APthreshold,cellType)
ylabel('Voltage (mV)'); title('AP threshold');
subplot(2,3,2)
boxplot(latency,cellType)
ylabel('Latency (ms)'); title('Latency to first AP');
subplot(2,3,3)
boxplot(APamplitude,cellType)
ylabel('Voltage (mV)'); title('AP amplitude');
subplot(2,3,4)
boxplot(APhalfwidth,cellType)
ylabel('Time (ms)'); title('AP halfwidth');
ylim([-2 2]);
subplot(2,3,5)
boxplot(maxtotFR,cellType)
ylabel('Spikes/s'); title('max total firing rate');
subplot(2,3,6)
boxplot(InstFR,cellType)
ylabel('Spikes/s'); title('Max Instantaneous firing rate');

saveas(gcf,'APPropertiesComparison.png');



%% Clustering analysis with all of the parameters

% 1) Generate matrix with the different variables and their corresponding
% values for each cell

% Extract the variables from the data
for i = 1:length(Data)
    Cage(:,i) = Data{1,i}(1).Cage;
    mouseID(:,i) = Data{1,i}(1).mouseID;
    Date(:,i) = Data{1,i}(1).Date;
end

Parameters = horzcat(maxtotFR',APamplitude',APthreshold',Ihold',latency',Vrest',Cage',mouseID',Date');

% Make vector with variable names

Variables = char('maxtotFR','APamplitude','APthreshold','Ihold','latency','Vrest','Cage','mouseID','Date');

% 2) Run PCA analysis

[coeff,score,latent,tsquared,explained,mu] = pca(Parameters);

% 3) Plot the first two components

figure, biplot(coeff(:,1:2),'scores',score(:,1:2),'Varlabels',Variables);

figure, plot(score(:,1),score(:,2),'ro');
xlabel('1st Principal Component');
ylabel('2nd Principal Component');
title('PCA analysis');

figure()
pareto(explained)
xlabel('Principal Component');
ylabel('Variance Explained (%)');
title('Variance explained per component');
