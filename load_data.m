my_init;
%% Load files
addpath('foam_data');
firstLine   = 6;
numColumns  = 3;
numFiles    = 10;

dataset =  'C'; % 'D'; %

for iFile=1:numFiles
    fileName = ['FOAM',num2str(iFile),dataset,'.txt']; % read all as strings
    D = textread(fileName, '%s','delimiter','\n'); % read full file as strings
    V = cellfun(@(s) sscanf(s,'%f').',D, 'un', 0); % read each string in D as numbers
    data_foam{iFile} = vertcat(V{:}); % concentrate arrays
%     fileName = [num2str(iFile),dataset];
%     fileData = data_foam{iFile};
%     save(fileName, 'fileData');
end
% fileAll = ['all_data_',dataset]
% save(fileAll,'data_foam');
%% Plots
indSample = [1000:1500];
figure;
id = 1;
for iFile=1:5
    subplot(5,2,id)
    plot(indSample,data_foam{iFile}(indSample,3),'LineWidth',0.5); hold on;           %data_foam{iFile}(indSample,1),
    names(iFile) = {[dataset,num2str(iFile)]};
    xlim([1000 1500]);
    ylim([-7 -4.999]);
%     title(names(iFile));
    xlabel('Time, sec')
    ylabel('$\Delta x$, mm')
    id = id + 1;
    subplot(5,2,id)
    plot(indSample,data_foam{iFile}(indSample,2),'LineWidth',0.5); hold on;           %data_foam{iFile}(indSample,1),
%     title(names(iFile));
    xlim([1000 1500]);
    ylim([-60 0]);
    xlabel(' Time, sec')
    ylabel(' Load, kN')
    id = id + 1;
end
%%
tikzName = ['data_',dataset,'_1to5.tikz'];
cleanfigure;
matlab2tikz(tikzName, 'showInfo', false,'parseStrings',false,'standalone', ...
            false, 'height', '18cm', 'width','12cm','checkForUpdates',false);

%%
figure;
id = 1;
for iFile=6:10
    subplot(5,2,id)
    plot(indSample,data_foam{iFile}(indSample,3),'LineWidth',0.5); hold on;           %data_foam{iFile}(indSample,1),
    names(iFile) = {[dataset,num2str(iFile)]};
    xlim([1000 1500]);
    ylim([-7 -4.999]);
%     title(names(iFile));
    xlabel('Time, sec')
    ylabel('$\Delta x$, mm')
    id = id + 1;
    subplot(5,2,id)
    plot(indSample,data_foam{iFile}(indSample,2),'LineWidth',0.5); hold on;           %data_foam{iFile}(indSample,1),
%     title(names(iFile));
    xlim([1000 1500]);
    ylim([-405 0]);
    xlabel('Time, sec')
    ylabel('Load, kN')
    id = id + 1;
end
%%
tikzName = ['data_',dataset,'_6to10.tikz'];
cleanfigure;
matlab2tikz(tikzName, 'showInfo', false,'parseStrings',false,'standalone', ...
            false, 'height', '18cm', 'width','12cm','checkForUpdates',false);
