% ANALYZING THE CC RESULTS

clear all;close all;

% choose directory of cell to be analyzed
CellPath = uigetdir('\\files.med.harvard.edu\Neurobio\MICROSCOPE\Melanie\ephys\Mel','Choose folder')

% set the current path to that directory
path = cd(CellPath);

%% extracting the data into one array

files = dir('*.mat'); %list the .mat files in the cwd 

vars = cell(size(files)); %initialize a cell the size of the number of files we have

for ii = 1:numel(files) %for each file
   fields{ii} = files(ii).name; %save the name of the files 
end

[cs,index] = sort_nat(fields);

for ii = 1:numel(files) %for each file
   fields{ii} = cs{ii}(1:end-4); %substract the .mat from the name
   vars{ii} = load(cs{ii}); %in each cell of the array 'vars' load the struct with the data from each file
   data{ii} = vars{ii,1}.(fields{ii}).data; %in a new cell array called 'data', store in each cell the info corresponding to the waveforms in each file. 
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

figure, set(gcf,'units','points','position',[100,100,1000,600]);
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
saveas(gcf,'AllTraces.svg');

%% Find peaks (APs) in the individual responses

% Find the APs with an automatic function
for ii = 1:size(responses,2)
 [peakLoc{ii},peakMag{ii}] = peakfinder(responses{ii},[],-10); % I added 10 as a threshold above which the data should be considered to look for peaks, cause otherwise it gives me weird peaks that are not there
end

% plot the data with the APs extracted and add or substract peaks that
% haven't been properly identified
% figure
% scrsz = get(0, 'ScreenSize');
% hf=figure(200); clf; set(gcf,'units','points','position',scrsz);
% hold on
% 
%     for count = 1:size (peakLoc,2)
%         NumTrial = count;
%         clf;
%         set(gcf,'Toolbar','figure','name',[ 'Trial ' num2str(count) ]);
%         
%         okMax = 0;
%         DelMax = 0;
%         AddMax = 0;
% 
%         bla = linspace(1,5,50000);
% %         indStart=Indices(NumeroCangrejo,count).indStart;
% %         indEnd=Indices(NumeroCangrejo,count).indEnd;
%         subplot(5,1,1:4)
%         hold on
%         plot(responses{count},'color','b','linewidth',2);
%         shg
%         hold on
%         locs = peakLoc{count};
%         mag = peakMag{count};
%         plot(peakLoc{count},peakMag{count},'ro');
%         
% %         for j=1:length(peakLoc)
% %             text(bla(locs(j)),responses{1,j},num2str(j),'FontSize',18)
% %         end
% 
%          %%%%%%%%%%%%
%         uicontrol('Style','pushbutton','String','Add peaks','CallBack','AddMax = 1;','Position',[scrsz(3)*1/10 scrsz(4)*1/10 120 60]);
%         uicontrol('Style','pushbutton','String','Delete peaks','CallBack','DelMax = 1;','Position',[scrsz(3)*5/10 scrsz(4)*1/10 120 60]);
%         uicontrol('Style','pushbutton','String','Done!','CallBack','okMax = 1','Position',[scrsz(3)*9/10 scrsz(4)*1/10 120 60]);
% 
%         while okMax == 0
%             figure(200)
%             title({'Waiting input',[' trial: ' num2str(count)]} )
%             
%              if DelMax ==1
%                 title('Del Max')
%                 inputTitle = char('Do you want to delete any peak?');
%                 prompt = {'Which ones? (write the numbers separated with a space):'};
%                 answer = inputdlg(prompt, inputTitle, 1);
%                 mov_to_add= str2num(answer{1});
%                 peakLoc{count}(mov_to_add)=[];
%                 
%                 DelMax=2;
%                 
%                 elseif DelMax == 2
%                 subplot(5,1,1:4)
%                 cla
%                 hold on
%                 
%                 plot(responses{count},'color','b','linewidth',2);
%                 shg
%                 hold on
%                 peakLoc = peakLoc{count};
%                 plot(peakLoc{count},peakMag{count},'ro');
% %                 for j=1:length(peakLoc)
% %                     text(bla(locs(j)+indStart-1),datosECGAll{NumeroCangrejo,NumTrial}(locs(j)+indStart-1),num2str(j),'FontSize',18)
% %                 end
%                 DelMax = 0;
%             end
%                 if AddMax == 1
%                 disp('hola mundo')
%                 clear dcm_obj sacc_add2
%                 dcm_obj = datacursormode(hf);
%                 sacc_add2 =  getCursorInfo(dcm_obj);
%                 ind_add=zeros(1, length(sacc_add2));
%                 for j=1:length(sacc_add2);
%                     ind_add(j)=sacc_add2(j).DataIndex(1);
%                 end
%                 if ind_add
%                     ind_add=sort(ind_add);
%                     startx=ind_add;
%                     AddMax=0;
%                     sacc_add.startx=(startx);
%                     [ind_start]=unique(sort([(locs) sacc_add.startx]));
%                     peakLoc{count}=ind_start;
%                     clear vector
%                 end
%                 AddMax=2;
%             elseif AddMax==2
%                 subplot(5,1,1:4)
%                 cla
%                 hold on
%                 plot(responses{count},'color','b','linewidth',2);
%                 shg
%                 hold on
%                 locs=peakLoc{count};
%                 plot(peakLoc{count},peakMag{count},'ro');
% %                 for j=1:length(locs)
% %                     text(bla(locs(j)+indStart-1),datosECGAll{NumeroCangrejo,NumTrial}(locs(j)+indStart-1),num2str(j),'FontSize',18)
% %                 end
%                 AddMax=0;
%             end
%             drawnow
%         end
%         
%         if ~exist('temp','dir')
%             mkdir temp
%         end
%         fileTemp=['./temp/temp' date num2str(round(rem(now,1))) 'trial' num2str(count)];
%         save(fileTemp)
%     end
%     
%     close all
%     
% pathNameSaeDatos = pwd; filesave='CheckedAPs.mat';
% if ~exist(fullfile(pathNameSaeDatos, filesave),'file')   
%     save(fullfile(pathNameSaeDatos, filesave))
% else    
%     inputTitle = filesave;
%     prompt = {['Data file already exists, change name (add.mat)']};
%     answer = inputdlg(prompt, inputTitle, 1);
%     newfilesave= answer{1};
%     save(fullfile(pathNameSaeDatos, newfilesave))
% end
    
for ii =1:size(responses,2)
figure,  set(gcf,'units','points','position',[100,100,1000,600]);
subplot(2,1,1)
plot(correctedPulses{ii},'k')
title('Current pulses'), ylabel('Current (pA)');
ylim([-400,1000]);

subplot(2,1,2)
plot(responses{ii},'k')
hold on
plot(peakLoc{ii},peakMag{ii},'ro')
title('Responses'), ylabel('Voltage (mV)'); xlabel('Time (s)');
ylim([-120,60]);
shg
end

%% AP analysis

close all;

pulseStart = 5000; % when in the protocol is our pulse starting (in points)
pulseEnd = 15000; % when is it ending
threshold = 10; % define a threshold to look for AP

% find trace with first AP
myData = cell2mat(responses);
myData2 = reshape(myData,[30000,length(pulses)]);
[row,col] = find(myData2(pulseStart:pulseEnd,:)>threshold); %look for every trace the points above threshold
sweepwithfirstAP = myData2(:,col(1)); % take the first column in X (= first sweep) with points above threshold
sweepnumberwithfirstAP = col(1); % which waveform number is the first with an AP?

% APthreshold
[row,col,v] = find(sweepwithfirstAP(pulseStart:pulseEnd)>threshold); % in the first sweep with an AP, find the points above threshold
AP = row+(pulseStart-1);
thresholdsearch  = sweepwithfirstAP((AP(1)-15):AP(1));
distance = (thresholdsearch(2:16)-thresholdsearch(1:15));
[row] = find(distance>3);
thresholdlocation=AP(1)-(16-row(1));
APthreshold = sweepwithfirstAP(thresholdlocation);

% latency to first AP
latency = rdivide((AP(1)-(pulseStart-1)),10);

% AP amplitude
APamplitude = minus(max(sweepwithfirstAP(thresholdlocation:thresholdlocation+10)),APthreshold);
[row] = find(round(sweepwithfirstAP(pulseStart:pulseEnd)) == round(APamplitude+APthreshold));
APlocation = row(1)+(pulseStart-1);

% AP 1/2 width
[row] = find(sweepwithfirstAP(thresholdlocation:(thresholdlocation+10))>(rdivide(APamplitude,2)+APthreshold));
P2 = (row(1))+thresholdlocation-1;
P1 = P2-1;
[row] = find(sweepwithfirstAP(APlocation:APlocation+20)<(rdivide(APamplitude,2)+APthreshold));
P4 = row(1)+APlocation-1;
P3 = P4-1;
yAPhalfwidth = APthreshold+rdivide(APamplitude,2);
APhalfwidth = (rdivide((yAPhalfwidth-sweepwithfirstAP(P3)),(sweepwithfirstAP(P4)-sweepwithfirstAP(P3)))+P3)-(rdivide((yAPhalfwidth-sweepwithfirstAP(P1)),(sweepwithfirstAP(P2)-sweepwithfirstAP(P1)))+P1);
xvectorforplot = [((rdivide((yAPhalfwidth-sweepwithfirstAP(P1)),(sweepwithfirstAP(P2)-sweepwithfirstAP(P1)))+P1)-thresholdlocation+11),((rdivide((yAPhalfwidth-sweepwithfirstAP(P3)),(sweepwithfirstAP(P4)-sweepwithfirstAP(P3)))+P3)-thresholdlocation+11)];
yvectorforplot = [yAPhalfwidth,yAPhalfwidth];
APhalfwidth = rdivide(APhalfwidth,10);

%AHP
   
 


    %% IV curve

for ii = 1:size(responses,2)
    pulseV(ii) = mode(round(responses{1,ii}(pulseStart:pulseEnd)));
    currents(ii) = correctedPulses{1,ii}(10300);
    DifCurrents(ii) = currents(ii) - correctedPulses{1,ii}(300); %add an extra correction, given that they start shifted from zero sometimes.
%make it read the difference instead.   
    firingRate(ii) = size(peakLoc{ii},2);
end


figure,  set(gcf,'units','points','position',[100,100,1000,600]);
subplot(1,2,1)
plot(DifCurrents,pulseV,'ro')
hold on
plot(DifCurrents,pulseV,'r')
title('IV curve');
xlim([-200 450]); ylim([-150 50]);
ylabel('V (mV)');xlabel('Current (pA)');

%If curve

subplot(1,2,2);
plot(DifCurrents,firingRate,'b^')
hold on
plot(DifCurrents,firingRate,'b')
title('If curve');
xlim([0 450]); ylim([0 max(firingRate)+1]);
ylabel('Firing rate (spikes/s)');xlabel('Current (pA)');

saveas(gcf,'IVandIFcurves.png');
saveas(gcf,'IVandIFcurves.svg');

%% first sweep with AP

% AP number in sweep with first AP
spikenumberfirsttrain = firingRate(sweepnumberwithfirstAP);

%rheobase
rheobase = DifCurrents(sweepnumberwithfirstAP);
%% maximal frequency
maxfreq = max(firingRate);

% sweep with max. freq.
sweepwithmaxfreq = find(firingRate == maxfreq);
sweepwithmaxfreq = sweepwithmaxfreq(1);
figure;
plot(firingRate,'o')
title('Spikes per sweep indicating sweep with max. freq.');
xlabel('Sweep number'); ylabel('Number of APs');
xlim([0 length(responses(1,:))]); ylim([0 max(firingRate)+1]);
hold on
plot(sweepwithmaxfreq,firingRate(sweepwithmaxfreq),'r','marker','*')

%% RC check

	% where the RC check occurs
	checkPulseSize = -200; % how big is the RC check
	checkPulseStart = 25000; % when does it start
	checkPulseEnd = 26000; % when does it end

    acqRate = 10; % points per ms

%     	if isempty(checkPulseStart)
% 					notPulse=[SR(1, pulseStart-10) SR(pulseEnd+150, acqLen-1)]; 
% 					rPeak=NaN;
% 					rEnd=NaN;
% 					tau=NaN;
% 				else
% 					if checkPulseStart>pulseStart % the RC check comes late
% 						notPulse=[SR(myData2,1, pulseStart-10) SR(myData2,pulseEnd+150, checkPulseStart-10)]; 
% 					else
% 						notPulse=[SR(myData2,1, checkPulseStart-1) SR(myData2,checkPulseEnd+50, pulseStart-10) SR(myData2,pulseEnd+150, acqLen-1)]; 
% 					end
% 					[rPeak, rEnd, tau]=csCurrentClampPulseAnalysis(SR(myData2,checkPulseStart, checkPulseEnd)- mean(notPulse), acqRate, checkPulseSize);
% 				end
                
    
    % Check this out to calculate the Rm
    
%     for i = 1:length(responses)
%     notPulse1{i} = [responses{1,i}(1:pulseStart-10)];
%     notPulse2{i} = [responses{1,i}(pulseEnd+150:end-1)];
%     end
%     
%     restMean = mean(notPulse);
%     pulseRm = 1000*(pulseV-restMean)/currents;

    for i = 1:length(responses)
[rPeak(i), rEnd(i), tau(i)] = csCurrentClampPulseAnalysis(mode(round(responses{1,i}(checkPulseStart:checkPulseEnd))), acqRate, checkPulseSize);   
    end 

    figure,
a4rPeak = subplot(1, 3, 1);
			title(a4rPeak, ['RC R']);
			xlabel(a4rPeak,'acq'); 
			ylabel(a4rPeak,'MO');
            ylim([0 700]);
			hold on
			a4rPulse=subplot(1, 3, 2);
			title(a4rPulse, ['Pulse R']);
			xlabel(a4rPulse, 'acq') 
			hold on
			a4tau=subplot(1, 3, 3);
			title(a4tau, ['RC Tau']);
			xlabel(a4tau,'acq') 
			ylabel(a4tau,'ms')
			hold on
            
				plot(a4rPeak, rPeak, 'go')
				plot(a4rPeak, rEnd, 'gx')
				%plot(a4rPulse, goodTraces, newCell.pulseRm(goodTraces), 'go')		
				plot(a4tau, tau, 'go')               
    
%% save image with sweep with first AP and sweep with max freq

figure,  set(gcf,'units','points','position',[100,100,1000,600]);
subplot(2,1,1)
plot(myData2(:,sweepwithmaxfreq),'k');axis([0 length(myData2) (max(myData2(:,sweepwithmaxfreq))-150) (max(myData2(:,sweepwithmaxfreq))+10)])
ylabel('Voltage (mV)');
hold on
text(17000,-20,strcat('Max AP number = ', num2str(maxfreq)),'FontSize',12);
title('Sweep with max AP number','Fontsize',12,'FontName','Calibri'); ylabel('Voltage (mV)');
subplot(2,1,2)
hold on
text(17000,-20,strcat('Rheobase = ', num2str(round(rheobase)),' pA'),'FontSize',12);
plot(myData2(:,sweepnumberwithfirstAP),'k');axis([0 length(myData2) (max(myData2(:,sweepwithmaxfreq))-150) (max(myData2(:,sweepwithmaxfreq))+10)])
ylabel('Voltage (mV)');
title('Sweep with first AP','Fontsize',12,'FontName','Calibri'); ylabel('Voltage (mV)');

saveas(gcf,'MaxFreqAndRheobase.png');
saveas(gcf,'MaxFreqAndRheobase.png');

    
%% Vrest

% Given that I am using the Ihold, it doesn't make sense to read the Vrest,
% but rather we open a window to input the value of Vrest without holding
% current, and/or the value of injected current

prompt = {'Enter Vrest:','Enter Ihold:'};
title = 'Cell resting properties';
dims = [1 35];
answer = inputdlg(prompt,title,dims);

Vrest = answer{1};
Ihold = answer{2};

%% Quality control

% 	% values for QC inclusion of individual sweeps
% 	maxRestSD=100;
%  	minRm=0;
%  	maxRm=1000;
% 	maxVm=2000;
% 	minVm=-2000;


%% Save output to excel table
A={'output',Vrest,Ihold,APthreshold,APamplitude,APhalfwidth,spikenumberfirsttrain,maxfreq,sweepwithmaxfreq,latency,sweepnumberwithfirstAP,rheobase};
xlswrite(strcat(CellPath,'table.xls'),A,1);

%% Save outputs to struct

% Generate struct
cellProp = struct('rheobase',rheobase,'pulses',correctedPulses,'responses',responses,'peakLoc',peakLoc,'peakMag',peakMag,'voltage',pulseV,'currents',DifCurrents,'firingRate',firingRate,'Vrest',Vrest,'Ihold',Ihold,'APthreshold',APthreshold,'Latency', latency, 'APamplitude',APamplitude, 'APhalfwidth',APhalfwidth, 'maxfreq',maxfreq);

% Save in the cell's folder
save('cellProp');
