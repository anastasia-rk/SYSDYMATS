my_init;
%% Load files
addpath('foam_data');
firstLine = 6;
num_columns = 3;

dataset = 'C'; %'D';

% switch dataset
%     case 'C'
%         % paste procedure from below
%     case 'D'
%          % paste procedure from below
% end

for iFile=1:10
    fileName = ['FOAM',num2str(iFile),dataset,'.txt']; % 
    D = textread(fileName, '%s','delimiter','\n'); % read full file as strings
    V = cellfun(@(s) sscanf(s,'%f').',D, 'un', 0); % read each string in D as numbers
    data_foam{iFile} = vertcat(V{:}); % concentrate arrays
end
