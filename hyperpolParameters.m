function [hyperpolsteadystate,sag] = hyperpolParameters(DifCurrents,responses,pulseStart,pulseEnd)

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

end