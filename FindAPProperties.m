function [APlocation,APpeak,APthreshold,latency,APthrough,APhalfwidth,APamplitude,sweepnumberwithfirstAP] = FindAPProperties(pulseStart,pulseEnd,responses,threshold,peakLoc,peakMag)
%AP analysis

% Find trace with first AP
myData = cell2mat(responses); %change the responses to a double (for now a vector)
myData2 = reshape(myData,[30000,length(responses)]); %change the vector of responses into a matrix where each column is one trace
[~,SweepsWithAP] = find(myData2(pulseStart:pulseEnd,:)>threshold); %look for every trace the points above threshold
sweepwithfirstAP = myData2(:,SweepsWithAP(1)); % take the first column in myData2 (= first sweep) with points above threshold
sweepnumberwithfirstAP = SweepsWithAP(1); % which waveform number is the first with an AP?

% % Interpolate linearly datapoints in the waveform
% InterpolatedAP = interp1(sweepwithfirstAP,1:0.1:length(sweepwithfirstAP));
% [APlocation,APpeak] = peakfinder(InterpolatedAP,[],-10); 

%Save location and values for the first AP of the first trace with APs
APlocation = peakLoc{1,sweepnumberwithfirstAP}(1);
APpeak = peakMag{1,sweepnumberwithfirstAP}(1);

% % APthreshold
% % Calculated as the maximum of the second derivative
% responseAroundFirstAP = sweepwithfirstAP(APlocation-200:APlocation+200); %keep data points around the first AP
% secDerFirstAP = diff(responseAroundFirstAP,2); %second derivative
% %secDerFirstAP = [0;0;secDerFirstAP]; %shifting the vector by 2 to correct the lag for taking the second derivative
% [MagMaxDiff,IdxMaxDiff] = max(secDerFirstAP);
% APthreshold = responseAroundFirstAP(IdxMaxDiff); %the threshold is the value that the Voltage takes in the position of the max of the second derivative
% thresholdlocation = APlocation-(201-IdxMaxDiff);
% plot((thresholdlocation ),APthreshold,'*')
%agregar gui para elegir el umbral manual

figure, plot(sweepwithfirstAP);xlim([APlocation-200,APlocation+200]);
[thresholdlocation,APthreshold] = ginput(1);
waitforbuttonpress;
thresholdlocation = round(thresholdlocation);

% latency to first AP
latency = (thresholdlocation-pulseStart)/10; %the latency is the position of the threshold since the pulse started, divided by 10 to get it in ms.

% AP amplitude
APamplitude= APpeak-APthreshold;

% AP through
%(minimum value of the membrane potential between the peak and the next AP)
if numel(peakLoc{1,sweepnumberwithfirstAP})>1 % if there are at least two peaks, take the minimum between them
    APthrough = min(sweepwithfirstAP(APlocation:peakLoc{1,sweepnumberwithfirstAP}(2)));
else % otherwise take the minimum in the next few ms
    APthrough = min(sweepwithfirstAP(APlocation:APlocation+100));
end
[APthroughLocation] = find(sweepwithfirstAP == APthrough);
APthroughloc = find(APthroughLocation > APlocation);
APthroughLocation = APthroughLocation(APthroughloc(1));


% The half-width of the AP is the width (in time units) of the trace at
% half height. Here we considered the spike amplitude from the threshold,
% so we will take the half height as the half amplitude (but we might want
% to recompute things).

% 1) Calculate the half-height (or in this case half-amplitude) of the
% action potential
halfAmplitude = APthreshold + APamplitude/2;

% 2) Look for the closest points at either side of the peak that match that height
halfAmpInLocation = knnsearch(sweepwithfirstAP(APlocation-20:APlocation),halfAmplitude);
halfAmpEndLocation = knnsearch(sweepwithfirstAP(APlocation:APlocation+50),halfAmplitude);

% 3) Define the width as the time interval between them
APhalfwidth = (20-halfAmpInLocation + halfAmpEndLocation)/10;


% Plot everything

% Plot it and add stuff to it as the analysis progresses
figure, plot(sweepwithfirstAP);xlim([5000,9500]); ylim([APthrough-5,APpeak+3]);
xlabel('Time (points)');ylabel('Voltage (mV)'); 
title('First trace with an AP (stimulation time only)');
hold on

plot((thresholdlocation ),APthreshold,'*')
plot([pulseStart,pulseStart+(latency*10)],[APthreshold,APthreshold]);
plot([APlocation,APlocation],[APthreshold,APpeak],'g');
plot(APthroughLocation,APthrough,'bo')
plot([APlocation-(20-halfAmpInLocation),APlocation+halfAmpEndLocation],[halfAmplitude,halfAmplitude])

end