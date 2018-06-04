function [responses,correctedPulses] = DefinePulsesAndResponses(data)
%% Define the responses and the pulses

IndexResponses = find(contains(data(2,:),'AD0')); %find the data corresponding to the responses
IndexPulses = find(contains(data(2,:),'AD4')); %find the data corresponding to the pulses

responses = data(1,IndexResponses); %save as responses the data for the A0 files. 
pulses = data(1,IndexPulses); %save as pulses the data for the A4 files. 

% Correct the pulses so that they are in pA
correctedPulses = cellfun(@(x) x*2000,pulses,'un',0);

for ii = 1:size(responses,2)
    pulseV(ii) = median(responses{1,ii}(5000:1500)); %calculate the voltage of the plateau, using the median of the round value of V during the pulse
    currents(ii) = correctedPulses{1,ii}(10300); %measure the current during a specific point of the pulse given
    DifCurrents(ii) = currents(ii) - correctedPulses{1,ii}(300); %add an extra correction, given that they start shifted from zero sometimes. Make it read the difference instead.   
end

end