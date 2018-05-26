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

%% IV and If curves with every cell

pulseStart=5000; % when in the protocol is our pulse starting (in points)
pulseEnd=15000; % when is it ending

% Extract the current pulses, neuron's voltage response and firing rates
% from the data
for i = 1:length(Data)
    voltage(:,i) = Data{1,i}.voltage;
    current(:,i) = Data{1,i}.currents;
    firingRate(:,i) = Data{1,i}.firingRate;
end

% Plot the IV curve
figure,  set(gcf,'units','points','position',[100,100,1000,600]);
subplot(1,2,1)
plot(current,voltage,'o')
hold on
plot(current,voltage)
title('IV curve');
xlim([-200 450]); ylim([-150 50]);
ylabel('V (mV)');xlabel('Current (pA)');

% Plot the IF curve
subplot(1,2,2);
plot(current,firingRate,'^')
hold on
plot(current,firingRate)
title('If curve');
xlim([0 450]);
ylabel('Firing rate (spikes/s)');xlabel('Current (pA)');

saveas(gcf,'IVandIFcurvesEveryCell.png');
saveas(gcf,'IVandIFcurvesEveryCell.svg');

%% Plot with cells grouped by dominance

%If we saved the folders as Dcellx and Scellx, all of the Dominant
%subfolders are going to be saved first, and the subordinate second

% 1) Reshape the voltage and current data into vectors
regroupedV = reshape(voltage,[size(voltage,1)*size(voltage,2),1]);
regroupedC = reshape(current,[size(current,1)*size(current,2),1]);
regroupedFR = reshape(firingRate,[size(firingRate,1)*size(firingRate,2),1]);

% 2) Generate a char vector with information about the group
    % a) Check to see how many D and S subfolders we have
    names = dir;
    for i = 1:length(names)
        FolderNames{i} = names(i).name;
    end
    
    dominant = strfind(FolderNames,'Dcell'); dominantNum = sum(cell2mat(dominant));
    subordinate = strfind(FolderNames,'Scell'); subordinateNum = sum(cell2mat(subordinate));
    
    % b) Make the vector
groups = cell(size(regroupedV,1),1);
groups(1:(size(voltage,1)*dominantNum),1) = {'Dominant'};
groups((size(voltage,1)*dominantNum)+1:end,1) = {'Subordinate'};

% Plot IV and IF curves sorted by dominance
figure,  set(gcf,'units','points','position',[100,100,1000,600]);
subplot(1,2,1)
gscatter(regroupedC,regroupedV,groups)
hold on
plot(current,voltage)
title('IV curve');
xlim([-200 450]); ylim([-150 50]);
ylabel('V (mV)');xlabel('Current (pA)');

subplot(1,2,2);
gscatter(regroupedC,regroupedFR,groups)
hold on
plot(current,firingRate)
title('If curve');
xlim([0 450]);
ylabel('Firing rate (spikes/s)');xlabel('Current (pA)');

saveas(gcf,'IVandIFcurvesByStatus.png');
saveas(gcf,'IVandIFcurvesByStatus.svg');

%% Compare the resting properties

% Extract the variables from the data
for i = 1:length(Data)
    Vrest(:,i) = str2num(Data{1,i}(1).Vrest);
    Ihold(:,i) = str2num(Data{1,i}(1).Ihold);
end

% Make boxplots
    % group by cellnumber
cellType = cell(size(Data,2),1);
cellType(1:dominantNum,1) = {'Dominant'};
cellType(dominantNum+1:end,1) = {'Subordinate'};

figure,  set(gcf,'units','points','position',[100,100,1000,600]);
subplot(1,2,1)
boxplot(Vrest,cellType)
ylabel('Voltage (mV)'); title('Vrest values');
subplot(1,2,2)
boxplot(Ihold,cellType)
ylabel('Current (pA)'); title('Iholding current values');

saveas(gcf,'CellPropertiesComparison.png');
saveas(gcf,'CellPropertiesComparison.svg');

%% Compare AP properties

% Extract the variables from the data
for i = 1:length(Data)
    rheobase(:,i) = Data{1,i}.rheobase;
    APthreshold(:,i) = Data{1,i}.APthreshold;
    latency(:,i) = Data{1,i}.Latency;
    APamplitude(:,i) = Data{1,i}.APamplitude;
    APhalfwidth(:,i) = Data{1,i}.APhalfwidth;
    maxfreq(:,i) = Data{1,i}.maxfreq;
end

% Make boxplots

figure,  set(gcf,'units','points','position',[100,100,1000,600]);
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
ylabel('Voltage (mV)'); title('AP halfwidth');
subplot(2,3,5)
boxplot(rheobase,cellType)
ylabel('Current (pA)'); title('Rheobase');
subplot(2,3,6)
boxplot(maxfreq,cellType)
ylabel('Current (pA)'); title('Max firing freq');

saveas(gcf,'APPropertiesComparison.png');
saveas(gcf,'APPropertiesComparison.svg');
