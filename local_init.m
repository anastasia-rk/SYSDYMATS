%% Add mylibs to and subdirectories to path in local directory
% Requires that mylibs is located in the parent folder.
LibFolder = '../mylibs';
addpath(LibFolder)
f = dir(LibFolder);
git = 4;
N = length(f)-git;
for iDir=1:N
    if f(iDir+git).isdir
       addpath([LibFolder,'/',f(iDir+git).name]);
    end
end
my_init;