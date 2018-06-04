% ANALYZING THE CC RESULTS

clear all;close all;

% choose directory of cell to be analyzed
CellPath = uigetdir() 
% set the current path to that directory
path = cd(CellPath);

files = dir('*.mat'); %list the .mat files in the cwd 
vars = cell(size(files)); %initialize a cell the size of the number of files we have

for ii = 1:numel(files) %for each file
   fields{ii} = files(ii).name; %save the name of the files 
end

[OrderedFileNames,index] = sort_nat(fields); %order the files by name in ascending order
waveformFiles = strfind(OrderedFileNames,'AD'); %find the files corresponding to waveforms
OrderedFileNames = OrderedFileNames(1:sum(cell2mat(waveformFiles)));

for ii = 1:numel(files) %for each file
   vars{ii} = load(OrderedFileNames{ii}); %in each cell of the array 'vars' load the struct with the data from each file
   fieldsID{ii} = OrderedFileNames{ii}(1:end-4); %substract the .mat from the name
   fieldsID{ii} = regexprep(OrderedFileNames{ii},'.mat',''); %substract the .mat from the name to use later on for reading the struct names
   data{1,ii} = vars{ii,1}.(fieldsID{ii}).data(1:30000); %in a new cell array called 'data', store in each cell the info corresponding to the waveforms in each file. 
   data{2,ii} = OrderedFileNames{ii};
end

% for i=1:size(data,2)
%     data{1,i} = data{1,i}(1:30000);
% end

% % Read the ephys files from a chosen directory
% [data,CellPath] = readEphysFiles(false);

%% Define pulses and responses

[responses,correctedPulses] = DefinePulsesAndResponses(data);

%% plot the raw results, all together

%figure, set(gcf,'units','points','position',[100,100,1000,600]); %if I run it in lab
figure, set(gcf,'units','points','position',[80,80,600,350]); %if I run it in my laptop
subplot(2,1,1)
hold on
cellfun(@plot,correctedPulses)
title('Current pulses'), ylabel('Current (pA)');
ylim([-300, 900]);
subplot(2,1,2)
hold on
cellfun(@plot,responses)
title('Responses'), ylabel('Voltage (mV)'); xlabel('Time (ms)');

saveas(gcf,fullfile(CellPath,'AllTraces.png'));
%saveas(gcf,'AllTraces.svg');

%% Find peaks (APs) in the individual responses

% Find the APs with an automatic function, and if you want, plot the
% individual results
[peakLoc,peakMag] = myPeakFinder(correctedPulses,responses,true,0);

%% AP analysis

close all;

pulseStart = 5000; % when in the protocol is our pulse starting (in points)
pulseEnd = 15000; % when is it ending
threshold = 10; % define a threshold to look for AP

[APlocation,APpeak,APthreshold,latency,APthrough,APhalfwidth,APamplitude,sweepnumberwithfirstAP] = FindAPProperties(pulseStart,pulseEnd,responses,threshold,peakLoc,peakMag);

    %% IV curve

% Plot the IV curve using the "plotIVcurve" function
[pulseV,DifCurrents] = plotIVcurve(responses,pulseStart,pulseEnd,correctedPulses,0);

%% Hyperpolarization parameters
% Calculate the hypolarized steady state and the sag using the
% hyperpolParameters function

[hyperpolsteadystate,sag] = hyperpolParameters(DifCurrents,responses,pulseStart,pulseEnd);

%% first sweep with AP

for ii = 1:size(responses,2)
APnum(ii) = size(peakLoc{ii},2);
end
% AP number in sweep with first AP
spikenumberfirsttrain = APnum(sweepnumberwithfirstAP);
%rheobase
rheobase = DifCurrents(sweepnumberwithfirstAP);
%% Firing rate

[maxAPnum,InstFR,totFR,maxInstFR,maxtotFR] = firingRateAnalysis(APnum,peakLoc);

%% If curves with the three different frequencies

plotIFcurves(DifCurrents,APnum,InstFR,maxInstFR,totFR,maxtotFR,0);

%% Load session's properties and save outputs to struct

[cellProp] = restingProp(CellPath,APthrough,InstFR,totFR,sag,maxtotFR,maxInstFR,rheobase,correctedPulses,responses,peakLoc,peakMag,pulseV,DifCurrents,APnum,APthreshold,latency,APamplitude,APhalfwidth,maxAPnum);
