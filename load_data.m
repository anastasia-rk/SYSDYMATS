my_init;
%% Load files
addpath('foam_data');
firstLine   = 6;
numColumns  = 3;
numFiles    = 10;

dataset = 'D'; % 'C'; %

% switch dataset
%     case 'C'
%         % paste procedure from below
%     case 'D'
%          % paste procedure from below
% end

for iFile=1:numFiles
    fileName = ['FOAM',num2str(iFile),dataset,'.txt']; % read all as strings
    D = textread(fileName, '%s','delimiter','\n'); % read full file as strings
    V = cellfun(@(s) sscanf(s,'%f').',D, 'un', 0); % read each string in D as numbers
    data_foam{iFile} = vertcat(V{:}); % concentrate arrays
end

%% Plots
figure;
for iFile=1:5
    plot(data_foam{iFile}(:,1),data_foam{iFile}(:,2)); hold on;
    names(iFile) = {[dataset,num2str(iFile)]};
end
title('Loads');
legend(names(1:5));

figure;
for iFile=1:5
    plot(data_foam{iFile}(:,1),data_foam{iFile}(:,3)); hold on;
end
title('Disp')
legend(names(1:5));

figure;
for iFile=6:10
    plot(data_foam{iFile}(:,1),data_foam{iFile}(:,2)); hold on;
    names(iFile) = {[dataset,num2str(iFile)]};
end
title('Loads')
legend(names(6:10));

figure;
for iFile=6:10
    plot(data_foam{iFile}(:,1),data_foam{iFile}(:,3)); hold on;
end
title('Disp')
legend(names(6:10));
