function [] = plotIVcurve(pulseV,DifCurrents,monitor)
% This function gives back the IV curve. The inputs it need are the
% responses, the pulseStart, the pulseEnd, the correctedPulses, and the
% monitor. If monitor is set to 1, it is lab, otherwise it is laptop
a=1

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