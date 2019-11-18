function [folderName] = make_folder(name,varargin)
    % To create a new folder with the folderName indicative of the settings
    % The first input variable should be a string array containing settings
    % names
    % The rest of the input variables are the values of corresponding
    % settings
    for i=2:nargin-1
        name = [name,'_',varargin{1}{i-1},'_',num2str(varargin{i})];

    end
    folderName = name;
    if ~exist(folderName, 'dir')
       mkdir(folderName)
    end
end