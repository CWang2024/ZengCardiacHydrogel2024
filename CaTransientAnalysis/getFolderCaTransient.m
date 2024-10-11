clc;
clear;
close all;

Folder = uigetdir(); % Full path of the folder extracted
if isempty(Folder)
    disp("User cancelled the program.");
end
Files = dir(Folder);
idx = strfind(Folder, '/'); % Find all occurrences of '/'
if ~isempty(idx)
    lastSlash = idx(end); % Get the position of the last '/'
    temp = Folder(lastSlash+1:end);
    % Extract only the folder name
else
    disp("Please check for folder picking.");
    return;
end

FileList=[];
for n=1:length(Files)
    if(contains(Files(n,1).name, '.mp4')||contains(Files(n,1).name, '.avi')||contains(Files(n,1).name, '.AVI'))
        FileList(end+1)=n;
    end
end

Files=Files(FileList);
T=struct2table(Files);
%T=natsortrows(T);
Files=table2struct(T);

number_files=size(Files,1);

% Create an output folder
Output_folder_name=strcat(temp, '_Output');

if exist(Output_folder_name)
    disp("Output folder already existed, scanning existing folder ... ... ");
    excelFiles = dir(fullfile(Output_folder_name, '*.xlsx'));
    Files=Files(length(excelFiles)+1: end);
    for i=1:length(Files)
        %Functions to read each file and write the excel, should take input
        %fullFileName, i, k, plot the data, take response, write the data.
        %path= '/Users/clarkwang/Documents/MATLAB/AP&Ca analysis/BatchCaAnalysis/videoCa analysis/para-10x.mp4';
        %filename='para-10x.mp4';
        %Output_folder_name='ZQJCa_Output';
        path=strcat(Files(i).folder, '/', Files(i).name);
        disp("Current file name: "+Files(i).name);
        get1Ca(path, Files(i).name, Output_folder_name);
    end
    disp('Conatenating excel sheets ... ...');
    concatenateAllVideoParameters(Output_folder_name, 'AllVideoParameters.xlsx', 1);
    concatenateAllVideoParameters(Output_folder_name, 'AllVideoParameters.xlsx', 'Raw dF');
    concatenateAllVideoParameters(Output_folder_name, 'AllVideoParameters.xlsx', 'Filtered dF');
    concatenateAllVideoParameters(Output_folder_name, 'AllVideoParameters.xlsx', 'Filtered dFF0');
    concatenateAllVideoParameters(Output_folder_name, 'AllVideoParameters.xlsx', 'AvgTrace');
    concatenateAllVideoParameters(Output_folder_name, 'AllVideoParameters.xlsx', 'dFF0AvgTrace');
else
    mkdir(Output_folder_name);
    %fullFileName = fullfile(Output_folder_name, baseFileName);
    cellInt=1;
    %%
    for i=1:length(Files)
        %Functions to read each file and write the excel, should take input
        %fullFileName, i, k, plot the data, take response, write the data.
        %path= '/Users/clarkwang/Documents/MATLAB/AP&Ca analysis/BatchCaAnalysis/videoCa analysis/para-10x.mp4';
        %filename='para-10x.mp4';
        %Output_folder_name='ZQJCa_Output';
        path=strcat(Files(i).folder, '/', Files(i).name);
        disp("Current file name: "+Files(i).name);
        get1Ca(path, Files(i).name, Output_folder_name);
    end
    %%
    disp('Conatenating excel sheets ... ...');
    concatenateAllVideoParameters(Output_folder_name, 'AllVideoParameters.xlsx', 1);
    concatenateAllVideoParameters(Output_folder_name, 'AllVideoParameters.xlsx', 'Raw dF');
    concatenateAllVideoParameters(Output_folder_name, 'AllVideoParameters.xlsx', 'Filtered dF');
    concatenateAllVideoParameters(Output_folder_name, 'AllVideoParameters.xlsx', 'Filtered dFF0');
    concatenateAllVideoParameters(Output_folder_name, 'AllVideoParameters.xlsx', 'AvgTrace');
    concatenateAllVideoParameters(Output_folder_name, 'AllVideoParameters.xlsx', 'dFF0AvgTrace');
end
%}
disp("Analysis of all files finished!");

