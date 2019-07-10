my_init;
%% Load files
addpath('foam_data');
firstLine   = 6;
numColumns  = 3;
numFiles    = 10;

dataset = 'D'; % 'C'; %

for iFile=1:numFiles
    fileName = ['FOAM',num2str(iFile),dataset,'.txt']; % read all as strings
    D = textread(fileName, '%s','delimiter','\n'); % read full file as strings
    V = cellfun(@(s) sscanf(s,'%f').',D, 'un', 0); % read each string in D as numbers
    data_foam{iFile} = vertcat(V{:}); % concentrate arrays
    fileName = [num2str(iFile),dataset];
    fileData = data_foam{iFile};
    save(fileName, 'fileData');
end
fileAll = ['all_data_',dataset]
save(fileAll,'data_foam');
%% Plots
figure;
id = 1;
for iFile=1:5
    subplot(5,2,id)
    plot(data_foam{iFile}(:,1),data_foam{iFile}(:,2)); hold on;
    names(iFile) = {[dataset,num2str(iFile)]};
    title(names(iFile));
    ylabel('Load')
    id = id + 1;
    subplot(5,2,id)
    plot(data_foam{iFile}(:,1),data_foam{iFile}(:,3)); hold on;
    title(names(iFile));
    ylabel('Disp')
    id = id + 1;
end

figure;
id = 1;
for iFile=6:10
    subplot(5,2,id)
    plot(data_foam{iFile}(:,1),data_foam{iFile}(:,2)); hold on;
    names(iFile) = {[dataset,num2str(iFile)]};
    title(names(iFile));
    ylabel('Load')
    id = id + 1;
    subplot(5,2,id)
    plot(data_foam{iFile}(:,1),data_foam{iFile}(:,3)); hold on;
    title(names(iFile));
    ylabel('Disp')
    id = id + 1;
end
