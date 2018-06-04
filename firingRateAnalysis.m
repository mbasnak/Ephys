function [maxAPnum,InstFR,totFR,maxInstFR,maxtotFR] = firingRateAnalysis(APnum,peakLoc)
% This function calculates the instantaneous and total FR of the sweeps.
% The inputs it needs are the APnum and the peakLoc

% 1) Absolute number of APs during the 1 s stimulation
maxAPnum = max(APnum);

% sweep with max. APnum
sweepwithmaxfreq = find(APnum == maxAPnum);
sweepwithmaxfreq = sweepwithmaxfreq(1);

% 2) Instantaneous firing frequency
InstFR = cell(1,length(peakLoc)); % initialize an empty cell array with as many cells as sweeps
two_peaks = cell(1,length(peakLoc)); % do that again with a different name
nonEmpty = (~cellfun(@isempty,peakLoc)); % return a logical with 1s in the cells corresponding to sweeps that have APs

for i = 1:length(peakLoc) % for every sweep
    if nonEmpty(i) == 0 || numel(peakLoc{1,i}) == 1 % if there are no APs or there is only 1 AP
        two_peaks{1,i} = 0; % then set the corresponding cel in two_peaks to 0
        InstFR{1,i} = 0; % and the one in InstFR
    else % if there are more than 1 AP
        two_peaks{1,i} = peakLoc{i}(2) - peakLoc{i}(1); % then two_peaks is the difference in position of the first two peaks
        InstFR{1,i} = 2/(two_peaks{1,i}/10000); % calculate the instantaneous FR (in Hz)
    end
end

maxInstFR = max(cell2mat(InstFR)); % save the maximum Inst FR

% 3) Frequency taken between first and last spike of the train
totFR = cell(1,length(peakLoc));
all_peaks = cell(1,length(peakLoc));

for i = 1:length(peakLoc)
    if nonEmpty(i) == 0 || numel(peakLoc{1,i}) == 1
        all_peaks{1,i} = 0;
        totFR{1,i} = 0;
    else % If the sweep has more than 1 AP
        all_peaks{1,i} = peakLoc{i}(APnum(1,i)) - peakLoc{i}(1); % Calculate the difference in time frames between the first and the last AP of the sweep
        totFR{1,i} = APnum(1,i)/(all_peaks{1,i}/10000); % use that difference to calculate the total firing rate in Hz
    end
end

maxtotFR = max(cell2mat(totFR)); % Save the max total FR

end