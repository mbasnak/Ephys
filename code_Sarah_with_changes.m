% Sarah's code to analyze the ephys experiments

clear all; close all;

%(in this code, you need to change lines 4 through 12 if necessary)

outputcell = 'A0'; % I presume this is the channel that reads the cell's response
file = 'example'; % rename file!
hyperpolsweep = 1;
pAincrease = 20;% pA increase per sweep! (in my protocol, this varies, so I'll have to change it)
currentonset = 5001; % point where the current pulse starts
threshold = 10; % set the threshold to define an AP
firstcurrentinjection = -200; % value of the first current injection
sweepfirstposcurrent = 11; % number of sweeps with neg and 0 current
directory='\\files.med.harvard.edu\Neurobio\MICROSCOPE\Melanie\ephys\Mel\experiment\cage8\20180517\other\cell1\'; % directory you read from
contents = dir(strcat(directory,'AD0_*')); % list the files in the directory that correspond to the cell's responses

% use the following lines to make sure the files are opened in order
for i = 1:size(contents,1)
    names{i} = contents(i).name;
end

[cs,index] = sort_nat(names);

%% Plot responses

figure;
X = []; % generate an empty array
for i = 1:length(cs); % for every response trace
    filename = cs{i};
    f = load(strcat(directory,filename)); % open the waveform file
    response = struct2cell(f); response = response{1};
    response = struct2cell(response); response = response{1}(1,1:30000); % save the response values. I'm saving up to three seconds, because the pulses have different durations
    plot(response); hold on; % plot them
    X=[X;response]; % fill the empty array with a response in each row
end
X = X'; %transpose to matrix
%X=X(:,[15:30]); choose sweeps if you are only analyzing a certain number
title('All traces'); xlabel('Time (points)'); ylabel('Voltage (mV)');
%% Calculate parameters

currentend = 10000+currentonset-1; % Point where the current injection ends
equalhundredpA = rdivide(100,pAincrease); % I don't know what this division is for, and the rdivide isn't necessary
thresholdpersistentfiring = 0;

%in case that not all sweeps shall be analyzed
%X=X(:,1:24);

%V resting membrane potential for the first trace
Vrestingmembrane = min(X(200:(currentonset-3),1)); % Calculate Vrest (but I am injecting Ihold, so this will read about -70 mV always) for the first trace
Vrestingmembranemean = mean(X(200:(currentonset-3),1));
Vrestingmembranelocation = find(X(1:(currentonset-3),1) == Vrestingmembrane); % find the locations before the current pulse where Vrest is worth that
Y = X;
while X(Vrestingmembranelocation(1)+2,1)>(Vrestingmembrane+0.11);
    Y(Vrestingmembranelocation(1),1) = 100;
    Vrestingmembrane = min(Y(200:(currentonset-3),1));
    Vrestingmembranelocation = find(Y(1:(currentonset-3),1) == Vrestingmembrane);
end

% Plot first trace with resting membrane
figure
plot(X(:,1))
title('first trace with resting membrane potential, sag, input resistance, rebound')
hold on
plot([1,(currentonset-1)],[Vrestingmembrane,Vrestingmembrane],'r','linewidth',2)
plot([1,(currentonset-1)],[Vrestingmembranemean,Vrestingmembranemean],'g','linewidth',2)

% add input resistance
hyperpolsteadystate = min(X((currentonset+9000):currentend,hyperpolsweep)); % Look for the steady state at the end of the response to the current injection, for the hyperpolarizing sweeps
hyperpolsteadystatelocation = (currentonset+8999)+find(X((currentonset+9000):currentend,hyperpolsweep) == hyperpolsteadystate);
Y=X;
while X(hyperpolsteadystatelocation(1)-2,hyperpolsweep)>(hyperpolsteadystate+0.12);
    Y(hyperpolsteadystatelocation(1),hyperpolsweep)=100;
    hyperpolsteadystate=min(Y((currentonset+9000):currentend,hyperpolsweep));
    hyperpolsteadystatelocation=(currentonset+8999)+find(Y((currentonset+9000):currentend,hyperpolsweep)==hyperpolsteadystate);
end
inputresistance = rdivide((hyperpolsteadystate-Vrestingmembrane),firstcurrentinjection*10^-3);
plot([(currentonset+9000),currentend],[hyperpolsteadystate,hyperpolsteadystate],'r','linewidth',2)

% Add sag
%(The sag is a dip in the voltage. Is it always present? What ions
%contribute to it?)

%rebound and depolarizing sag only proper if no EPSPs present

sag = (min(X(currentonset:(currentonset+5000),hyperpolsweep)));
saglocation =(currentonset-1)+find(X(currentonset:(currentonset+5000),hyperpolsweep) == sag);
Y=X;
while X(saglocation(1)+2,hyperpolsweep)>(sag+0.11);
    Y(saglocation(1),hyperpolsweep) = 100;
    sag = min(Y(currentonset:(currentonset+5000),hyperpolsweep));
    saglocation = (currentonset-1)+find(Y(currentonset:(currentonset+5000),hyperpolsweep) == sag);
end
hold on
plot(saglocation,sag,'c','marker','*','MarkerSize', 12)
sag = hyperpolsteadystate-sag;
if sag>0
    sag = sag;
else
    sag = 0;
end
sag = -sag/(hyperpolsteadystate-Vrestingmembrane)*100;


% add rebound
Vrestingmax = max(X(1:(currentonset-3),hyperpolsweep));
Vrestingmaxlocation = find(X(1:(currentonset-3),hyperpolsweep)==Vrestingmax);
Y=X;
while X(Vrestingmaxlocation(1)+2,hyperpolsweep)<(Vrestingmax-0.11);
    Y(Vrestingmaxlocation(1),hyperpolsweep) = Vrestingmembrane;
    Vrestingmax = max(Y(1:(currentonset-3),hyperpolsweep));
    Vrestingmaxlocation = find(Y(1:(currentonset-3),hyperpolsweep) == Vrestingmax);
end
    
rebound = max(X(currentend:(currentonset+12000),hyperpolsweep));
reboundlocation = currentend-1+find(X(currentend:(currentonset+12000),hyperpolsweep)== rebound);
Y=X;
while X(reboundlocation(1)+2,hyperpolsweep)<(rebound-0.11);
        Y(reboundlocation(1),hyperpolsweep) = -100;
        rebound = max(Y(currentend:(currentonset+11999),hyperpolsweep));
        reboundlocation = currentend-1+find(Y(currentend:(currentonset+11999),hyperpolsweep) == rebound);
end


if max(X(currentonset:(currentonset+12000),hyperpolsweep))>threshold;
    rebound = [];
else
    plot((reboundlocation(1)),rebound,'g','marker','*','MarkerSize', 12)
    hold on
   % plot(Vrestingmaxlocation(1),Vrestingmax,'g','marker','*','MarkerSize',
   % 12) % why plotting this?
    rebound = rebound-Vrestingmax;
end



%% AP analysis

%find trace with first AP
[row,col] = find(X(currentonset:currentend,:)>threshold); %look for every trace the points above threshold
sweepwithfirstAP = X(:,col(1)); % take the first column in X (= first sweep) with points above threshold
sweepnumberwithfirstAP = col(1); % which waveform number is the first with an AP?

%APthreshold
[row,col,v] = find(sweepwithfirstAP(currentonset:currentend)>threshold); % in the first sweep with an AP, find the points above threshold
AP = row+(currentonset-1);
thresholdsearch  = sweepwithfirstAP((AP(1)-15):AP(1));
distance = (thresholdsearch(2:16)-thresholdsearch(1:15));
[row] = find(distance>3);
thresholdlocation=AP(1)-(16-row(1));
APthreshold = sweepwithfirstAP(thresholdlocation);

%latency to first AP
latency=rdivide((AP(1)-(currentonset-1)),10)

%AP amplitude
APamplitude = minus(max(sweepwithfirstAP(thresholdlocation:thresholdlocation+10)),APthreshold);
[row] = find(sweepwithfirstAP(currentonset:currentend) > APamplitude+APthreshold);
APlocation = row(1)+(currentonset-1);

%AP 1/2 width
[row]=find(sweepwithfirstAP(thresholdlocation:(thresholdlocation+10))>(rdivide(APamplitude,2)+APthreshold));
P2=(row(1))+thresholdlocation-1;
P1=P2-1;
[row]=find(sweepwithfirstAP(APlocation:APlocation+20)<(rdivide(APamplitude,2)+APthreshold));
P4=row(1)+APlocation-1;
P3=P4-1;
yAPhalfwidth=APthreshold+rdivide(APamplitude,2);
APhalfwidth=(rdivide((yAPhalfwidth-sweepwithfirstAP(P3)),(sweepwithfirstAP(P4)-sweepwithfirstAP(P3)))+P3)-(rdivide((yAPhalfwidth-sweepwithfirstAP(P1)),(sweepwithfirstAP(P2)-sweepwithfirstAP(P1)))+P1);
xvectorforplot=[((rdivide((yAPhalfwidth-sweepwithfirstAP(P1)),(sweepwithfirstAP(P2)-sweepwithfirstAP(P1)))+P1)-thresholdlocation+11),((rdivide((yAPhalfwidth-sweepwithfirstAP(P3)),(sweepwithfirstAP(P4)-sweepwithfirstAP(P3)))+P3)-thresholdlocation+11)];
yvectorforplot=[yAPhalfwidth,yAPhalfwidth];
APhalfwidth=rdivide(APhalfwidth,10);
% 
% %AHP amplitude
% if APlocation<(currentend-100)
%     beginofsecondAPorendofcurrent=find(sweepwithfirstAP((APlocation+100):currentend)<APthreshold);
%     beginofsecondAPorendofcurrent=find(sweepwithfirstAP(beginofsecondAPorendofcurrent(1)+APlocation+99:currentend)>APthreshold);
%     if sum(beginofsecondAPorendofcurrent)>0;
%         beginofsecondAPorendofcurrent=beginofsecondAPorendofcurrent(1)+APlocation-1;
%     else
%         beginofsecondAPorendofcurrent=currentend;
%     end
%     AHPamplitude=APthreshold-(min(sweepwithfirstAP(APlocation:beginofsecondAPorendofcurrent)))
%     xAHPpeak=find(sweepwithfirstAP(thresholdlocation:beginofsecondAPorendofcurrent)==(-AHPamplitude+APthreshold));
%     xAHPpeak=xAHPpeak(1)+thresholdlocation-1;
%     %AHP time to peak
%     AHPtimetopeak=rdivide(xAHPpeak-APlocation,10);
%     
%     %first AHP (in case that largest AHP is not first AHP)
%     [r]=findpeaks(-(sweepwithfirstAP(APlocation:xAHPpeak)));
%     if sum(r)>0;
%         firstAHPamplitude=APthreshold+r(1)
%         xfirstAHP=find(sweepwithfirstAP(APlocation:xAHPpeak)==(-r(1)));
%         xfirstAHP=xfirstAHP(1)+APlocation-1;
%     else
%         firstAHPamplitude=AHPamplitude
%         xfirstAHP=xAHPpeak
%     end
%     %DAP and second AHP
%     amppotentialDAPhyperpollist=[];
%     if beginofsecondAPorendofcurrent-xfirstAHP>=1
%         potentialDAPdepol=findpeaks(sweepwithfirstAP((xfirstAHP):beginofsecondAPorendofcurrent));
%         for i=1:length(potentialDAPdepol)
%             xpotentialDAPdepol=find(sweepwithfirstAP(xfirstAHP:beginofsecondAPorendofcurrent)==potentialDAPdepol(i));
%             amppotentialDAPhyperpol=potentialDAPdepol(i)-(min(sweepwithfirstAP((xpotentialDAPdepol(1)+xfirstAHP-1):beginofsecondAPorendofcurrent)));
%             if amppotentialDAPhyperpol>0.25;
%                 amppotentialDAPhyperpol=amppotentialDAPhyperpol;
%             else
%                 amppotentialDAPhyperpol=0;
%             end
%             amppotentialDAPhyperpollist=[amppotentialDAPhyperpollist,amppotentialDAPhyperpol];
%         end
%         amppotentialDAPhyperpollist=find(amppotentialDAPhyperpollist);
%         if sum(amppotentialDAPhyperpollist)==0;
%             DAPdepolamplitude=0
%             DAPhyperpolamplitude=0
%             figure
%             plot(sweepwithfirstAP((thresholdlocation-10):(thresholdlocation+500)))
%         else
%             amppotentialDAPhyperpollist=amppotentialDAPhyperpollist(1);
%             xpotentialDAPdepol=find(sweepwithfirstAP(xfirstAHP:beginofsecondAPorendofcurrent)==potentialDAPdepol(amppotentialDAPhyperpollist));
%             xpotentialDAPdepol=xpotentialDAPdepol(1)+xfirstAHP-1;
%             yDAPhyperpol=min(sweepwithfirstAP(xpotentialDAPdepol:beginofsecondAPorendofcurrent));
%             xDAPhyperpol=find(sweepwithfirstAP(xpotentialDAPdepol:beginofsecondAPorendofcurrent)==yDAPhyperpol);
%             xDAPhyperpol=xDAPhyperpol(1)+xpotentialDAPdepol-1;
%             yDAPdepol=max(sweepwithfirstAP(xfirstAHP:xDAPhyperpol));
%             xDAPdepol=find(sweepwithfirstAP(xfirstAHP:xDAPhyperpol)==yDAPdepol);
%             xDAPdepol=xDAPdepol(1)+xfirstAHP-1;
%             if xDAPhyperpol<xfirstAHP+1000; %increased from 200 to 1000 for ca1 pyramidal cells
%                 if xDAPdepol<xfirstAHP+1000;%increased from 100 to 1000 for ca1 pyramidal cells
%                     DAPdepolamplitude=yDAPdepol-(APthreshold-firstAHPamplitude)
%                     DAPhyperpolamplitude=yDAPdepol-yDAPhyperpol
%                     %plot AP and AHP
%                     figure
%                     plot(sweepwithfirstAP((thresholdlocation-10):(xDAPhyperpol+50)))
%                     hold on
%                     plot((xDAPhyperpol-thresholdlocation+11),(yDAPhyperpol),'y','marker','*')
%                     hold on
%                     plot((xDAPdepol-thresholdlocation+11),(yDAPdepol),'g','marker','.')
%                 else
%                     DAPdepolamplitude=0
%                     DAPhyperpolamplitude=0
%                     figure
%                     plot(sweepwithfirstAP((thresholdlocation-10):(thresholdlocation+500)))
%                 end
%             else
%                 DAPdepolamplitude=0
%                 DAPhyperpolamplitude=0
%                 figure
%                 plot(sweepwithfirstAP((thresholdlocation-10):(thresholdlocation+500)))
%                 
%             end
%         end
%     else
%         DAPdepolamplitude=0
%         DAPhyperpolamplitude=0
%         figure
%         plot(sweepwithfirstAP((thresholdlocation-10):(thresholdlocation+500)))
%     end
%     hold on
%     title('sweep with first AP indicating APthreshold, AP, first AHP, DAP and second AHP')
%     plot(11,APthreshold,'r','marker','.') 
%     plot(APlocation-thresholdlocation+11,APamplitude+APthreshold,'g','marker','.')
%     plot(xvectorforplot,yvectorforplot,'r','linewidth',3)
%     plot((xAHPpeak-thresholdlocation+11),(APthreshold-AHPamplitude),'k','marker','.')
%     plot((xfirstAHP-thresholdlocation+11),(APthreshold-firstAHPamplitude),'y','marker','*')
%     
%     figure
%     plot(sweepwithfirstAP)
%     title('sweep with first AP indicating max. AHP and latency')
%     hold on
%     plot((xAHPpeak),(APthreshold-AHPamplitude),'y','marker','.')
%     plot([1,((latency*10)+(currentonset-1))],[10,10],'g')
%     
% else
%     AHPamplitude=NaN
%     xAHPpeak=NaN
%     figure
%     plot(sweepwithfirstAP)
%     title('sweep with first AP indicating max. AHP and latency')
%     hold on
%     plot([1,((latency*10)+(currentonset-1))],[10,10],'g')
%     DAPdepolamplitude=NaN
%     DAPhyperpolamplitude=NaN
%     AHPtimetopeak=NaN
%     firstAHPamplitude=NaN
%     xfirstAHP=NaN
%     
% end

%%  spike frequency in all sweeps
R = X;
for numsweep=1:length(R(1,:)); % for every sweep
    for numsamplingpoint = 1:length(R); % for each data point
        if R(numsamplingpoint,numsweep)>threshold; % if the data point is above threshold
            R(numsamplingpoint,numsweep) = R(numsamplingpoint,numsweep); % keep the actual value
        else
            R(numsamplingpoint,numsweep) = threshold; %otherwise set it to threshold
        end
    end
end
spikecount=[];
for numsweep=1:length(R(1,:)); % for every sweep
    [peaks]=findpeaks(R(currentonset:currentend,numsweep)); % find peaks
    spikecount=[spikecount,length(peaks)];
end
figure;
plot(spikecount,'o')
title('spikes per sweep indicating sweep with max. freq.');
xlabel('Sweep number'); ylabel('Number of APs');
xlim([0 length(R(1,:))]); ylim([0 max(spikecount)+1]);

%% 

%AP number in sweep with first AP
spikenumberfirsttrain = spikecount(sweepnumberwithfirstAP);

%maximal frequency
maxfreq = max(spikecount);

%sweep with max. freq.
sweepwithmaxfreq = find(spikecount==maxfreq);
sweepwithmaxfreq = sweepwithmaxfreq(1);
hold on
plot(sweepwithmaxfreq,spikecount(sweepwithmaxfreq),'r','marker','*')

% spikecount first 50 ms after the current injection
spikecount50=[];
for numsweep=1:length(R(1,:));
    [peaks]=findpeaks(R(currentonset:(currentonset+500),numsweep));
    spikecount50=[spikecount50,length(peaks)];
end
figure;
plot(spikecount50,'o')
title('Spikes per sweep during first 50 ms of current injection indicating sweep with max. freq.')
xlabel('Sweep number'); ylabel('Number of APs');
xlim([0 length(R(1,:))]); ylim([0 max(spikecount50)+1]);

%% mean ISI in all sweeps
R = X;
for numsweep = 1:length(R(1,:)); % for every sweep
    for numsamplingpoint = 1:length(R); % for each data point
        if R(numsamplingpoint,numsweep)>threshold;
            R(numsamplingpoint,numsweep) = R(numsamplingpoint,numsweep);
        else
            R(numsamplingpoint,numsweep)=(threshold);
        end
    end
end
meanISIs=[];
for numsweep=1:length(R(1,:));
    [row,col]=findpeaks(R(currentonset:currentend,numsweep));
    ISImean = rdivide(mean(col(2:length(col))-col(1:(length(col)-1))),10);
    meanISIs=[meanISIs,ISImean];
end
figure
plot(meanISIs,'c*')
title('Mean ISI per sweep')
xlabel('Sweep number'); ylabel('Interspike interval (ms)');

%% 

% Initial freq (calculated from ISI of sweep with max freq)
[row,col]=findpeaks(R(currentonset:currentend,sweepwithmaxfreq));
initialISI=rdivide((col(2)-col(1)),10)
figure
plot(X(:,sweepwithmaxfreq))
title('initial/middle/final ISI')
hold on
plot(col(1)+currentonset,row(1),'r','marker','*')
plot(col(2)+currentonset,row(2),'r','marker','*')


%freq at 200 ms 
[row,col]=findpeaks(R(currentonset:currentend,sweepwithmaxfreq));
[r]=find(col>(currentonset+2000));
if length(r)>1;
    middleISI=rdivide((col(r(2))-col(r(1))),10)
    hold on
    plot(col(r(1))+(currentonset-1),row(r(1)),'g','marker','*')
    plot(col(r(2))+(currentonset-1),row(r(2)),'g','marker','*')
else
    middleISI=NaN
end


%final freq.
[row,col]=findpeaks(R(currentonset:currentend,sweepwithmaxfreq));
[r]=find(col>(8000));
if length(r)>1;
    finalISI=rdivide((col(length(col))-col(length(col)-1)),10)
    hold on
    plot(col(length(row)-1)+(currentonset-1),row(length(row)-1),'y','marker','*')
    plot(col(length(row))+(currentonset-1),row(length(row)),'y','marker','*')
else
    finalISI=NaN
end

%% 

%total adaptation
totaladaptation=100-((1/finalISI)*100/(1/initialISI));

%persistent firing frequency in all sweeps
%(for neurons that keep firing after the current injection)
R=X;
for numsweep=1:length(R(1,:));
    for numsamplingpoint=1:length(R);
        if R(numsamplingpoint,numsweep)>thresholdpersistentfiring;
            R(numsamplingpoint,numsweep)=R(numsamplingpoint,numsweep);
        else
            R(numsamplingpoint,numsweep)=threshold;
        end
    end
end
persistentfiring=[];
for numsweep=1:length(R(1,:));
    [row]=findpeaks(R((currentend+15):length(R(:,1)),numsweep));
    persistentfiring=[persistentfiring,length(row)];
end
figure;
plot(persistentfiring)
title('spikes per sweep after current injection')

%% ratio 1st/2nd ISI in first sweep with >2 spikes
R=X;
for numsweep=1:length(R(1,:));
    for numsamplingpoint=1:length(R);
        if R(numsamplingpoint,numsweep)>threshold;
            R(numsamplingpoint,numsweep)=R(numsamplingpoint,numsweep);
        else
            R(numsamplingpoint,numsweep)=threshold;
        end
    end
end
if sweepnumberwithfirstAP+equalhundredpA<=sweepwithmaxfreq; %this might be wrong because of the "equalhundred pA"
    [row,col]=findpeaks(R(currentonset:currentend,(sweepnumberwithfirstAP+equalhundredpA)));
    if length(col)<=2
        ratiofirstsecondISI=NaN;
    else
        ratiofirstsecondISI=rdivide(rdivide((col(2)-col(1)),10),(rdivide((col(3)-col(2)),10)));
        figure
        plot(X(:,(sweepnumberwithfirstAP+equalhundredpA)))
        title('ratio 1st/2nd ISI')
        hold on
        plot(col(1)+(currentonset-1),row(1),'r','marker','*')
        plot(col(2)+(currentonset-1),row(2),'r','marker','*')
        plot(col(3)+(currentonset-1),row(3),'r','marker','*')
    end
else
    ratiofirstsecondISI=NaN
end
%% 

%rheobase
rheobase = ((sweepnumberwithfirstAP-sweepfirstposcurrent)*pAincrease);

%saturation current
currentformaxfreq = ((sweepwithmaxfreq-sweepfirstposcurrent)*pAincrease); % this is not right because pA increase is not always the same


%% Save output image and excel file

%save image with first sweep, sweep with first AP and sweep with max freq
figure,
subplot('Position', [0.1,0.74,0.6325,0.21],'XTickLabel',{'1000','2000','3000'},'XTick',[10001 20001 30001])
hold on
title('Sweep with max firing rate','Fontsize',10,'FontName','Calibri'); ylabel('Voltage (mV)');
plot(X(:,sweepwithmaxfreq),'k');axis([0 length(X) (max(X(:,sweepwithmaxfreq))-150) (max(X(:,sweepwithmaxfreq))+10)])
subplot('Position', [0.1,0.41,0.6325,0.21],'XTickLabel',{'1000','2000','3000'},'XTick',[10001 20001 30001])
hold on
title('Sweep with first AP','Fontsize',10,'FontName','Calibri'); ylabel('Voltage (mV)');
plot(X(:,sweepnumberwithfirstAP),'k');axis([0 length(X) (max(X(:,sweepwithmaxfreq))-150) (max(X(:,sweepwithmaxfreq))+10)])
subplot('Position', [0.1,0.09,0.6325,0.21],'XTickLabel',{'1000','2000','3000'},'XTick',[10001 20001 30001])
hold on
title('First sweep','Fontsize',10,'FontName','Calibri'); ylabel('Voltage (mV)');
plot(X(:,1),'k');axis([0 length(X) (max(X(:,sweepwithmaxfreq))-200) (max(X(:,sweepwithmaxfreq))-40)])
xlabel('Time (ms)');
print('-depsc',strcat(directory,'figure.eps'))

save(strcat(directory,'spikecount'),'spikecount');

%output to excel table
%A={file,Vrestingmembranemean,inputresistance,sag,rebound,APthreshold,APamplitude,APhalfwidth,AHPamplitude,AHPtimetopeak,firstAHPamplitude,DAPdepolamplitude,DAPhyperpolamplitude,spikenumberfirsttrain,maxfreq,sweepwithmaxfreq,initialISI,middleISI,finalISI,totaladaptation,latency,ratiofirstsecondISI,sweepnumberwithfirstAP,rheobase,currentformaxfreq};
A={file,Vrestingmembranemean,inputresistance,sag,rebound,APthreshold,APamplitude,APhalfwidth,spikenumberfirsttrain,maxfreq,sweepwithmaxfreq,initialISI,middleISI,finalISI,totaladaptation,latency,ratiofirstsecondISI,sweepnumberwithfirstAP,rheobase,currentformaxfreq};
xlswrite(strcat(directory,'table.xls'),A,1);