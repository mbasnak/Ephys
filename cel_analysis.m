% ANALYZING THE CC RESULTS

clear all;close all;

% choose directory of cell to be analyzed
CellPath = uigetdir() 

% set the current path to that directory
path = cd(CellPath);

%% extracting the data into one array

files = dir('*.mat'); %list the .mat files in the cwd 

vars = cell(size(files)); %initialize a cell the size of the number of files we have

for ii = 1:numel(files) %for each file
   fields{ii} = files(ii).name; %save the name of the files 
end

[OrderedFileNames,index] = sort_nat(fields); %order the files by name in ascending order

for ii = 1:numel(files) %for each file
   vars{ii} = load(OrderedFileNames{ii}); %in each cell of the array 'vars' load the struct with the data from each file
   fieldsID{ii} = OrderedFileNames{ii}(1:end-4); %substract the .mat from the name
   %fieldsID{ii} = regexprep(OrderedFileNames{ii},'.mat',''); %substract the .mat from the name to use later on for reading the struct names
   data{ii} = vars{ii,1}.(fieldsID{ii}).data; %in a new cell array called 'data', store in each cell the info corresponding to the waveforms in each file. 
end


%% Define the responses and the pulses

% Since the responses are always AD0 and the pulses are AD4, the responses
% are loaded first, so the first half of the cell array contains responses
% and the second half, pulses

halfSize = (size(data,2))/2;

for r = 1:halfSize
responses{r} = data{1,r}(1:30000);
pulses{r} = data{1,halfSize+r}(1:30000);
end

% Correct the pulses so that they are in pA
correctedPulses = cellfun(@(x) x*2000,pulses,'un',0);

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

saveas(gcf,'AllTraces.png');
%saveas(gcf,'AllTraces.svg');

%% Find peaks (APs) in the individual responses

% Find the APs with an automatic function
for ii = 1:size(responses,2)
 [peakLoc{ii},peakMag{ii}] = peakfinder(responses{ii},[],-10); % I added 10 as a threshold above which the data should be considered to look for peaks, cause otherwise it gives me weird peaks that are not there
end

% plot the data with the APs extracted and add or substract peaks that
% haven't been properly identified
% (agregar corrector de picos de Vero)
%     
% for ii =1:size(responses,2)
% %figure, set(gcf,'units','points','position',[100,100,1000,600]); %if I run it in lab
% figure, set(gcf,'units','points','position',[80,80,600,350]); %if I run it in my laptop
% subplot(2,1,1)
% plot(correctedPulses{ii},'k')
% title('Current pulses'), ylabel('Current (pA)');
% ylim([-400,1000]);
% 
% subplot(2,1,2)
% plot(responses{ii},'k')
% hold on
% plot(peakLoc{ii},peakMag{ii},'ro')
% title('Responses'), ylabel('Voltage (mV)'); xlabel('Time (s)');
% ylim([-120,60]);
% shg
% end

%% AP analysis

close all;

pulseStart = 5000; % when in the protocol is our pulse starting (in points)
pulseEnd = 15000; % when is it ending
threshold = 10; % define a threshold to look for AP

% Find trace with first AP
myData = cell2mat(responses);
myData2 = reshape(myData,[30000,length(pulses)]);
[row,col] = find(myData2(pulseStart:pulseEnd,:)>threshold); %look for every trace the points above threshold
sweepwithfirstAP = myData2(:,col(1)); % take the first column in myData2 (= first sweep) with points above threshold
sweepnumberwithfirstAP = col(1); % which waveform number is the first with an AP?
sweepnumberwithfirstAP_stimOnly = sweepwithfirstAP(pulseStart:pulseEnd);

% Plot it and add stuff to it as the analysis progresses
figure, plot(sweepnumberwithfirstAP_stimOnly)
xlabel('Time (points)');ylabel('Voltage (mV)'); xlim([0,2000]);
title('First trace with an AP (stimulation time only)');
hold on

% APthreshold
% Calculated as the maximum of the second derivative
firstDerFirstAP = diff(sweepnumberwithfirstAP_stimOnly); %first derivative
secDerFirstAP = diff(firstDerFirstAP); %second derivative
[row,col] = max(secDerFirstAP); %look for the max of the second derivative
APthreshold = sweepnumberwithfirstAP_stimOnly(col); %the threshold is the value that the Voltage takes in the position of the max of the second derivative
thresholdlocation = col + pulseStart - 1;
plot((thresholdlocation + pulseStart -1),APthreshold,'*')

% latency to first AP
latency = col/10; %the latency is the position of the threshold since the pulse started, divided by 10 to get it in ms.
plot([0,latency*10],[APthreshold,APthreshold]);

% AP amplitude
APamplitude = minus(max(sweepwithfirstAP(thresholdlocation:thresholdlocation+10)),APthreshold);
[row] = find(round(sweepwithfirstAP(pulseStart:pulseEnd)) == round(APamplitude+APthreshold));
APlocation = row(1)+(pulseStart-1);
plot([APlocation-pulseStart+1,APlocation-pulseStart+1],[-80,APamplitude]);

% AP through
%(minimum value of the membrane potential between the peak and the next AP)
if numel(peakLoc{1,sweepnumberwithfirstAP})>1 % if there are at least two peaks, take the minimum between them
    APthrough = min(sweepwithfirstAP(APlocation:peakLoc{1,sweepnumberwithfirstAP}(2)));
else % otherwise take the minimum in the next few ms
    APthrough = min(sweepwithfirstAP(APlocation:APlocation+40));
end
[APthroughLocation] = find(sweepwithfirstAP == APthrough);
APthroughLocation = APthroughLocation(1);

plot(APthroughLocation-pulseStart+1,APthrough,'bo')


% AP 1/2 width
% [row] = find(sweepwithfirstAP(thresholdlocation:(thresholdlocation+10))>(rdivide(APamplitude,2)+APthreshold));
% P2 = (row(1))+thresholdlocation-1;
% P1 = P2-1;
% [row] = find(sweepwithfirstAP(APlocation:APlocation+20)<(rdivide(APamplitude,2)+APthreshold));
% P4 = row(1)+APlocation-1;
% P3 = P4-1;
% yAPhalfwidth = APthreshold+rdivide(APamplitude,2);
% APhalfwidth = (rdivide((yAPhalfwidth-sweepwithfirstAP(P3)),(sweepwithfirstAP(P4)-sweepwithfirstAP(P3)))+P3)-(rdivide((yAPhalfwidth-sweepwithfirstAP(P1)),(sweepwithfirstAP(P2)-sweepwithfirstAP(P1)))+P1);
% xvectorforplot = [((rdivide((yAPhalfwidth-sweepwithfirstAP(P1)),(sweepwithfirstAP(P2)-sweepwithfirstAP(P1)))+P1)-thresholdlocation+11),((rdivide((yAPhalfwidth-sweepwithfirstAP(P3)),(sweepwithfirstAP(P4)-sweepwithfirstAP(P3)))+P3)-thresholdlocation+11)];
% yvectorforplot = [yAPhalfwidth,yAPhalfwidth];
% APhalfwidth = rdivide(APhalfwidth,10);

% The half-width of the AP is the width (in time units) of the trace at
% half height. Here we considered the spike amplitude from the threshold,
% so we will take the half height as the half amplitude (but we might want
% to recompute things).

% 1) Calculate the half-height (or in this case half-amplitude) of the
% action potential
halfAmplitude = APthreshold + APamplitude/2;

% 2) Look for the closest points at either side of the peak that match that height
[c halfAmpInLocation] = min(abs(sweepwithfirstAP(APlocation-20:APlocation)-halfAmplitude));
[d halfAmpEndLocation] = min(abs(sweepwithfirstAP(APlocation:APlocation+50)-halfAmplitude));

% 3) Define the width as the time interval between them
APhalfwidth = (halfAmpInLocation + halfAmpEndLocation)/10;
plot([APlocation-pulseStart+1-(20-halfAmpInLocation),APlocation-pulseStart+1+halfAmpEndLocation],[halfAmplitude,halfAmplitude])

%It seems as though I don't have enough data points as I would like.
% maybe I can do some sort of fit?


%AHP
   
 


    %% IV curve

for ii = 1:size(responses,2)
    pulseV(ii) = mode(round(responses{1,ii}(pulseStart:pulseEnd))); %calculate the voltage of the plateau, using the mode of the round value of V during the pulse
    pulseVmed(ii) = median(responses{1,ii}(pulseStart:pulseEnd)); %calculate it using the median instead
    currents(ii) = correctedPulses{1,ii}(10300);
    DifCurrents(ii) = currents(ii) - correctedPulses{1,ii}(300); %add an extra correction, given that they start shifted from zero sometimes. Make it read the difference instead.   
    APnum(ii) = size(peakLoc{ii},2);
end

%figure, set(gcf,'units','points','position',[100,100,1000,600]); %if I run it in lab
figure, set(gcf,'units','points','position',[80,80,600,350]); %if I run it in my laptop
plot(pulseV,DifCurrents,'ro')
hold on
plot(pulseV,DifCurrents,'r')
%plot(pulseVmed,DifCurrents,'b') they are really similar, so I'm using the mode
title('IV curve'); ylim([-250 600]); xlim([-120 0]);
xlabel('V (mV)');ylabel('Current (pA)');
hline = refline([0 0]);
hline.Color = 'k';

saveas(gcf,'IVcurve.png');


%% Hyperpolarization parameters
% Calculate the hypolarized steady state and the sag

%1) Identify hyperpolarizing pulses...
[~,pulseIdentity] = find(DifCurrents <0);

for i = 1:length(pulseIdentity) % For every hyperpolarizing pulse
    
%2) The hyperpolarized steady state
hyperpolsteadystate(i) = median(responses{1,i}((pulseStart+8000):pulseEnd)); % Look for the steady state at the end of the response to the current injection, for the hyperpolarizing sweeps
hyperpolsteadystatelocation{i} = (pulseStart+8999)+find(responses{1,i}((pulseStart+9000):pulseEnd) == hyperpolsteadystate(i));

% figure,
% plot(responses{1,i}) % plot the responses
% hold on
% plot([(pulseStart+8000),pulseEnd],[hyperpolsteadystate(i),hyperpolsteadystate(i)],'r','linewidth',2) % Overlay the steady state

%3) The sag
minV(i) = (min(responses{1,i}(pulseStart:(pulseStart+5000)))); %take the sag as the minimum value in the voltage during the pulse, for the hyperpolarizing sweeps.
saglocation{i} =(pulseStart-1)+find(responses{1,i}(pulseStart:(pulseStart+5000)) == minV(i)); %find the location of that minimum as the point number for the whole trace

% hold on
% plot(saglocation{i},minV(i),'c','marker','*','MarkerSize', 12)
int_sag(i) = (hyperpolsteadystate(i)-minV(i)); %calculate the sag as the difference between the steady state and the min value
if int_sag(i)>0 %if there is a sag (i.e., if there is a value lower than the steady state)
    int_sag(i) = int_sag(i);
else
    int_sag(i) = 0;
end
sag(i) = int_sag(i)/minV(i);
end

for i = 1:length(pulseIdentity)
hyperpolV(i) = mode(round(responses{1,i}(pulseStart:pulseEnd))); %Subset responses to hyperpolarizing pulses
end
figure,
plot(hyperpolV,sag,'ro'); hold on
plot(hyperpolV,sag,'r');
title('Sag ratio as a function of voltage');
xlabel('Voltage response (mV)'); ylabel('Sag ratio');

saveas(gcf,'Sagratio.png');

%% first sweep with AP

% AP number in sweep with first AP
spikenumberfirsttrain = APnum(sweepnumberwithfirstAP);

%rheobase
rheobase = DifCurrents(sweepnumberwithfirstAP);
%% Firing rate

% 1) Absolute number of APs during the 1 s stimulation
maxAPnum = max(APnum);

% sweep with max. APnum
sweepwithmaxfreq = find(APnum == maxAPnum);
sweepwithmaxfreq = sweepwithmaxfreq(1);

% 2) Instantaneous firing frequency
InstFR = cell(1,length(peakLoc));
two_peaks = cell(1,length(peakLoc));
nonEmpty = (~cellfun(@isempty,peakLoc));

for i = 1:length(peakLoc)
    if nonEmpty(i) == 0 || numel(peakLoc{1,i}) == 1
        two_peaks{1,i} = 0;
        InstFR{1,i} = 0;
    else
        two_peaks{1,i} = peakLoc{i}(2) - peakLoc{i}(1);
        InstFR{1,i} = 2/(two_peaks{1,i}/10000);
    end
end

maxInstFR = max(cell2mat(InstFR));

% 3) Frequency taken between first and last spike of the train
totFR = cell(1,length(peakLoc));
all_peaks = cell(1,length(peakLoc));

for i = 1:length(peakLoc)
    if nonEmpty(i) == 0 || numel(peakLoc{1,i}) == 1
        all_peaks{1,i} = 0;
        totFR{1,i} = 0;
    else
        all_peaks{1,i} = peakLoc{i}(APnum(1,i)) - peakLoc{i}(1);
        totFR{1,i} = APnum(1,i)/(all_peaks{1,i}/10000);
    end
end

maxtotFR = max(cell2mat(totFR));

%% If curves with the three different frequencies

%figure, set(gcf,'units','points','position',[100,100,1000,600]); %if I run it in lab
figure, set(gcf,'units','points','position',[80,80,600,350]); %if I run it in my laptop
subplot(1,3,1);
plot(DifCurrents,APnum,'b^')
hold on
plot(DifCurrents,APnum,'b')
title('If curve');
xlim([0 600]); ylim([0 max(APnum)+1]);
ylabel('Total AP number');xlabel('Current (pA)');

subplot(1,3,2);
plot(DifCurrents,cell2mat(InstFR),'b^')
hold on
plot(DifCurrents,cell2mat(InstFR),'b')
title('If curve');
xlim([0 600]); ylim([0 maxInstFR+1]);
ylabel('Instantenous firing rate (spikes/s)');xlabel('Current (pA)');

subplot(1,3,3);
plot(DifCurrents,cell2mat(totFR),'b^')
hold on
plot(DifCurrents,cell2mat(totFR),'b')
title('If curve');
xlim([0 600]); ylim([0 maxtotFR+1]);
ylabel('Total firing rate (spikes/s)');xlabel('Current (pA)');

saveas(gcf,'IFcurves.png');

%% Resting properties

% Trying to automate this
% 1) Load the spreadsheet with all the info
[num,txt,parameters] = xlsread('Z:\MICROSCOPE\Melanie\ephys\Mel\experiment\everyCell\summary.xlsx',1);

% 2) Identify the row corresponding to the cell in question
cellID = CellPath(54:end); %save the cell name from the path
cellName = strcmp(parameters,cellID); %find the cell with the cell name in the parameters array
cellPar = cell(1,size(parameters,2));
for i = 1:size(parameters,2)
cellPar{1,i} = parameters{find(cellName,1),i};%save the row corresponding to that cell
end

% 3) Save the parameters
Date = cellPar{1,2};
Cage = cellPar{1,3};
mouseID = cellPar{1,4};
Rin = cellPar{1,5};
Rs = cellPar{1,6};
Cm = cellPar{1,7};
Vrest = cellPar{1,8};
Ihold = cellPar{1,9};
Temp = cellPar{1,10};
Rsfin = cellPar{1,11};
Rinfin = cellPar{1,12};
Cmfin = cellPar{1,13};
Vrestfin = cellPar{1,14};

%% Save outputs to struct

% Generate struct
cellProp = struct('InstFR',InstFR,'totFR',totFR,'sag',sag,'Date', Date, 'Cage', Cage, 'mouseID', mouseID, 'Rs', Rs,'Temp',Temp,'Rsfin',Rsfin,'Rinfind',Rinfin,'Cmfin',Cmfin,'Vrestfin',Vrestfin,'maxtotFR',maxtotFR,'maxInstFR',maxInstFR,'rheobase',rheobase,'pulses',correctedPulses,'responses',responses,'peakLoc',peakLoc,'peakMag',peakMag,'voltage',pulseV,'currents',DifCurrents,'APnum',APnum,'Cm',Cm,'Rin',Rin,'Vrest',Vrest,'Ihold',Ihold,'APthreshold',APthreshold,'Latency', latency, 'APamplitude',APamplitude, 'APhalfwidth',APhalfwidth, 'maxAPnum',maxAPnum);

% Save in the cell's folder
save('cellProp');
