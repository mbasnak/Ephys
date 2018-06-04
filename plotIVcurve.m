function [pulseV,DifCurrents] = plotIVcurve(responses,pulseStart,pulseEnd,correctedPulses,monitor)
% This function saves the voltage and current data from the cell
% information and gives back the IV curve. The inputs it need are the
% responses, the pulseStart, the pulseEnd, the correctedPulses, and the
% monitor. If monitor is set to 1, it is lab, otherwise it is laptop

for ii = 1:size(responses,2)
    pulseV(ii) = median(responses{1,ii}(pulseStart:pulseEnd)); %calculate the voltage of the plateau, using the median of the round value of V during the pulse
    currents(ii) = correctedPulses{1,ii}(10300); %measure the current during a specific point of the pulse given
    DifCurrents(ii) = currents(ii) - correctedPulses{1,ii}(300); %add an extra correction, given that they start shifted from zero sometimes. Make it read the difference instead.   
end

if monitor == 1
figure, set(gcf,'units','points','position',[100,100,1000,600]); %if I run it in lab
else
figure, set(gcf,'units','points','position',[80,80,600,350]); %if I run it in my laptop
end

plot(pulseV,DifCurrents,'ro')
hold on
plot(pulseV,DifCurrents,'r')
title('IV curve'); ylim([-250 600]); xlim([-120 0]);
xlabel('V (mV)');ylabel('Current (pA)');
hline = refline([0 0]);
hline.Color = 'k';
saveas(gcf,'IVcurve.png');

end