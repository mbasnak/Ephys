function [peakLoc,peakMag] = myPeakFinder(correctedPulses,responses,plotMe,monitor)

%the way to call this function is myPeakFinder(correctedPulses,responses,plotMe,monitor)
%the plotMe input will either show me the plots or not. It is either true
%or false
%the monitor input will change the size of the figure acording to where I
%am. It is either 1 (for the lab monitor) or any other number (for the
%laptop)

for ii = 1:size(responses,2)
 [peakLoc{ii},peakMag{ii}] = peakfinder(responses{ii},[],-15); % I added 10 as a threshold above which the data should be considered to look for peaks, cause otherwise it gives me weird peaks that are not there
end

% plot the data with the APs extracted and add or substract peaks that
% haven't been properly identified
% (agregar corrector de picos de Vero)

if plotMe
    
if monitor == 1
figure, set(gcf,'units','points','position',[100,100,1000,600]); %if I run it in lab
else
figure, set(gcf,'units','points','position',[80,80,600,350]); %if I run it in my laptop
end

for ii =1:size(responses,2)
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
hold off
disp('Press enter for next plot');
waitforbuttonpress
end

end
end