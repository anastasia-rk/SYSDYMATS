my_init;
%% Load files
folderName = 'foam_2010';
addpath(folderName);
firstLine   = 6;
numColumns  = 3;
numFiles    = 10;

dataset = 'C'; % 'D'; %  

for iFile=1:numFiles
    fileName = ['FOAM',num2str(iFile),dataset,'.txt']; % read all as strings
    D = textread(fileName, '%s','delimiter','\n'); % read full file as strings
    V = cellfun(@(s) sscanf(s,'%f').',D, 'un', 0); % read each string in D as numbers
    data_foam{iFile} = vertcat(V{:}); % concentrate arrays
    fileName = [folderName,'/',num2str(iFile),dataset];
    fileData = data_foam{iFile};
    save(fileName, 'fileData');
end
% fileAll = ['all_data_',dataset]
% save(fileAll,'data_foam');
%% Plots
indSample = [1000:1500];
figure;
id = 1;
for iFile=1:numFiles
    subplot(5,2,id)
    plot(indSample,data_foam{iFile}(indSample,2),'LineWidth',0.5); hold on;           %data_foam{iFile}(indSample,1),
    names(iFile) = {[dataset,num2str(iFile), ' load, kN']};
    xlim([1000 1500]);
    xlabel('Sample index')
    ylabel(names(iFile));
    id = id + 1;
end
%%
tikzName = ['data_',dataset,'_outputs.tikz'];
cleanfigure;
matlab2tikz(tikzName, 'showInfo', false,'parseStrings',false,'standalone', ...
            false, 'height', '18cm', 'width','12cm','checkForUpdates',false);

%%
indSample = [1:2:1501];
figure;
plot(indSample,data_foam{iFile}(indSample,3),'LineWidth',0.5); hold on;           %data_foam{iFile}(indSample,1),
xlim([1 1500]);
% ylim([-7 -4.999]);
xlabel('Sample index')
ylabel('$\Delta x$, mm')

%%
tikzName = ['data_',dataset,'_input.tikz'];
cleanfigure;
matlab2tikz(tikzName, 'showInfo', false,'parseStrings',false,'standalone', ...
            false, 'height', '5cm', 'width','12cm','checkForUpdates',false);
