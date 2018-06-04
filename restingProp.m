function [cellProp] = restingProp(CellPath,APthrough,InstFR,totFR,sag,maxtotFR,maxInstFR,rheobase,correctedPulses,responses,peakLoc,peakMag,pulseV,DifCurrents,APnum,APthreshold,latency,APamplitude,APhalfwidth,maxAPnum)
% Load Excel spreadsheet with different session and cell parameters and
% save the information

% 1) Load the spreadsheet with all the info
%[num,txt,parameters] = xlsread('Z:\MICROSCOPE\Melanie\ephys\Mel\experiment\everyCell\summary.xlsx',1);
[num,txt,parameters] = xlsread('D:\Doctorado\rotations\Bernardo\Ephys\everyCell\summary.xlsx',1);

% 2) Identify the row corresponding to the cell in question
%cellID = CellPath(54:end); %save the cell name from the path
cellID = CellPath(49:end); %save the cell name from the path
cellName = strcmp(parameters,cellID); %find the cell with the cell name in the parameters array
cellPar = cell(1,size(parameters,2));
for i = 1:size(parameters,2)
cellPar{1,i} = parameters{find(cellName,1),i}; %save the row corresponding to that cell
end

% 3) Save the parameters
Date = cellPar{1,2};
Cage = cellPar{1,3};
mouseID = cellPar{1,4};
Rin = cellPar{1,5};
Rs = cellPar{1,6};
Cm = cellPar{1,7};
Vrest = cellPar{1,8};
Ihold = cellPar{1,9};
Temp = cellPar{1,10};
Rsfin = cellPar{1,11};
Rinfin = cellPar{1,12};
Cmfin = cellPar{1,13};
Vrestfin = cellPar{1,14};

% Generate struct
cellProp = struct('APthrough',APthrough,'InstFR',InstFR,'totFR',totFR,'sag',sag,'Date', Date, 'Cage', Cage, 'mouseID', mouseID, 'Rs', Rs,'Temp',Temp,'Rsfin',Rsfin,'Rinfin',Rinfin,'Cmfin',Cmfin,'Vrestfin',Vrestfin,'maxtotFR',maxtotFR,'maxInstFR',maxInstFR,'rheobase',rheobase,'pulses',correctedPulses,'responses',responses,'peakLoc',peakLoc,'peakMag',peakMag,'voltage',pulseV,'currents',DifCurrents,'APnum',APnum,'Cm',Cm,'Rin',Rin,'Vrest',Vrest,'Ihold',Ihold,'APthreshold',APthreshold,'Latency', latency, 'APamplitude',APamplitude, 'APhalfwidth',APhalfwidth, 'maxAPnum',maxAPnum);

cd(CellPath);
% Save in the cell's folder
save('cellProp');

end