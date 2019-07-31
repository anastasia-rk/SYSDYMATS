my_init;
%% Set up global parameters
dataset = 'C'; % 'D'; %                                                     % name of dataset
iFile   = 1;                                                                % id of the sample
% Length of input and output lags
n_u     = 2;                                                                % input signal lag length
n_y     = 3;                                                                % output signal lag length
d       = n_u + n_y;                                                        % size of input vector x
lambda  = 4;                                                                % order of polynomial
a       = sym('a',[1 d]);                                                   % associated symbolic vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Prelimianries
for iFile=1:K
%% Upload data
clear Input Output 
fileName = [num2str(iFile),dataset];
load(fileName);
Input  = fileData(:,2);
Output = fileData(:,3);
T = length(Input); % length of the observation sequence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create batch of input vectors
diff = n_u - n_y;                                                           % difference between lags
switch sign(diff)
    case 0
        disp('zero')
        t_0 = n_u;
    case -1
        disp('negative')
        t_0 = n_y;
    case 1
        disp('positive')
        t_0 = n_u;
end

iNarx = 0;                                                                  % batch index of the input vector in AR model
for t=t_0:T
    iNarx = iNarx + 1;
    x_narx{iFile}(:,iNarx) = [Input(t-n_u+1:t,1);Output(t-n_y+1:t,1)];      % NARX input
    y_narx{iFile}(:,iNarx) = Output(t,1);                                   % NARX output
end
nNarx = iNarx;                                                              % Length of NARX input batch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create sum index permutations
indeces{1} = [1:d]';
for iLambda=2:lambda                                                        % iLambda - order of polynomial term
    initial_m = permn(1:d,iLambda);                                         % get all permutations with repetition
    M = initial_m;                                                          % M - temporary set of permutations
    for j=iLambda:-1:2
        ind = find(M(:,j)>=M(:,j-1));                                       % sort out only increasing indeces
        M = M(ind,:);                                                       % update set
        clear ind
    end
    indeces{iLambda} = M;                                                   % all index combinations of order iLambda
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create batch of regressors
% Computing all polynomial terms
iTerm = 0;
for iLambda = 1:lambda
    for k = 1:size(indeces{iLambda},1)
        iTerm = iTerm + 1;  
        for iNarx = 1:nNarx 
            term{iFile}(iNarx,iTerm) = regressor(x_narx{iFile}(:,iNarx),indeces{iLambda}(k,:));         % compute the regressor
        end
        symb_term{iFile,iTerm} = a(indeces{iLambda}(k,:));                  % dictionary of regressors
    end
end
nTerms = iTerm;                                                              % total number of regressors
end                                                                         % loop over files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Select significant terms for all datasets
% Compute correlation coefficients
for iTerm = 1:nTerms                                                        % over all polynomial terms in the dictionary
    AEER{iTerm} = zeros(nTerm,1);                                           % placeholder for AERR criteria
    for iFile=1:K                                                           % over all datasets
        for jTerm = 1:nTerms                                                % over all polynomial terms in the dictionary
            cf(iFile,jTerm) = cor_sqr(y_narx{iFile},term{iFile}(:,jTerm));
            AEER{iTerm}(jTerm) = AEER{iTerm}(jTerm) + cf(iFile,jTerm);
        end
    end
    [AEER_max,iMax] = max(AEER{iTerm});
    S{iTerm} = iMax;
    phi{iTerm} = term{iFile}(:,iMax);
    
end

%




%% Set up and simulate the model
% 
% narx = @(u,y,n_u,n_y,Theta) Theta(1:n_u)*u(:,end-n_u+1:end)' + Theta(n_u+1:n_u+n_y)*y(:,end-n_y+1:end)';
% 
% theta = ones(1,10);