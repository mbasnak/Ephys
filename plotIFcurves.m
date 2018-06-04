function [] = plotIFcurves(DifCurrents,APnum,InstFR,maxInstFR,totFR,maxtotFR,monitor)

if monitor == 1
figure, set(gcf,'units','points','position',[100,100,1000,600]); %if I run it in lab
else
figure, set(gcf,'units','points','position',[80,80,600,350]); %if I run it in my laptop
end

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

end