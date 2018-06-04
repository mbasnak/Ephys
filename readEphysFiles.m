function [data,CellPath] = readEphysFiles(trimdata)

% choose directory of cell to be analyzed
CellPath = uigetdir() 
path = cd(CellPath);

files = dir(fullfile(CellPath,'*.mat')); %list the .mat files in the cwd
vars = cell(size(files)); %initialize a cell the size of the number of files we have

for ii = 1:numel(files) %for each file
   fields{ii} = files(ii).name; %save the name of the files 
end

[OrderedFileNames,index] = sort_nat(fields); %order the files by name in ascending order using a helper function
waveformFiles = strfind(OrderedFileNames,'AD'); %find the files corresponding to waveforms
OrderedFileNames = OrderedFileNames(1:sum(cell2mat(waveformFiles))); % keep only the names of the files corresponding to waveforms

for ii = 1:numel(OrderedFileNames) %for each file
   vars{ii} = load(OrderedFileNames{ii}); %in each cell of the array 'vars' load the struct with the data from each file
   fieldsID{ii} = regexprep(OrderedFileNames{ii},'.mat',''); %substract the .mat from the name to use later on for reading the struct names
   data{1,ii} = vars{ii,1}.(fieldsID{ii}).data; %in a new cell array called 'data', store in each cell the info corresponding to the waveforms in each file. 
   data{2,ii} = OrderedFileNames{ii};
end

if trimdata %if I give this input as true, cut the responses and pulses and make them 3 s
for i=1:size(data,2)
    data{1,i} = data{1,i}(1:30000);
end
end


end