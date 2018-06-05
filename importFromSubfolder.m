function [Data] = importFromSubfolder(CellPath)
%% Import the structs with relevant info from all the folders in the path

% Get list of all subfolders in the path
allSubFolders = genpath(CellPath);
% Parse into a cell array.
remain = allSubFolders;
listOfFolderNames = {}; % Initialize an empty cell array
while true
    % The strtok function returns in the output "singleSubFolder" the text in 
    % the input "remain" contained within the delimiter ";" until "remain"
    % is empty because it has already recovered all the folder names.
	[singleSubFolder, remain] = strtok(remain, ';');
	if isempty(singleSubFolder)
		break;
    end
    % It moves those folders into listOfFolderNames before singleSubFolder
    % gets empty
	listOfFolderNames = [listOfFolderNames singleSubFolder];
end
numberOfFolders = length(listOfFolderNames);

cellData = {};
for i= 2:length(listOfFolderNames) % For every folder in the path except from the one making the path ("everyCell")
    cellData{i-1} = load(strcat(listOfFolderNames{i}, '\cellProp.mat')); % Load the file with all the relevant info into a cell array 
end

% Extract specifically the relevant struct
for i = 1:size(cellData,2)
    Data{i} = cellData{1,i}.cellProp;
end

end