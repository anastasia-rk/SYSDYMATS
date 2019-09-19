my_init;
%% Set up global parameters
dataset =    'C'; %   'D'; %                                                % name of dataset
iFile   = 1;                                                                % id of the sample
K       = 10;                                                               % number of datasets
% Length of input and output lags
n_u     = 5;                                                                % input signal lag length
n_y     = 5;                                                                % output signal lag length
d       = n_y + n_u;                                                        % size of input vector x
lambda  = 5;                                                                % order of polynomial
a       = sym('a',[1 d]);                                                   % associated symbolic vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Identify difference in lag
diff = n_u - n_y;                                                           % difference between lags
switch sign(diff)
    case 0
        disp('zero')
        t_0 = n_u+1;
    case -1
        disp('negative')
        t_0 = n_y+1;
    case 1
        disp('positive')
        t_0 = n_u+1;
end
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
for iFile=1:K
disp(['Dataset_',num2str(iFile)])
%% Upload data
clear Input Output 
fileName = [num2str(iFile),dataset];
load(fileName);
Input  = fileData(:,2);
Output = fileData(:,3);
T = length(Input); % length of the observation sequence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create the batch of input vectors
iNarx = 0;                                                                  % batch index of the input vector in AR model
for t=t_0:T
    iNarx = iNarx + 1;
    x_narx(:,iNarx) = [Output(t-n_y:t-1,1);Input(t-n_u:t-1,1)];      % NARX input
end
nNarx = iNarx;                                                              % length of NARX input batch
y_narx(:,:) = Output(t_0:end);                                       % NARX output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create the batch of regressors and compute all polynomial terms
iTerm = 0;
if (iFile > 1)
    for iLambda = 1:lambda                                                  % iLambda - order of polynomial term
        for k = 1:size(indeces{iLambda},1)                                  % get index combination for the regressor in the sum
            iTerm = iTerm + 1;                                              % term index in the polynomial
            for iNarx = 1:nNarx                                             % going throught the full NARX batch
                term(iNarx,iTerm) = regressor(x_narx(:,iNarx),...
                                              indeces{iLambda}(k,:));       % compute the regressor (numeric)
            end
        end    
    end
else
    for iLambda = 1:lambda                                                  % iLambda - order of polynomial term
        for k = 1:size(indeces{iLambda},1)                                  % get index combination for the regressor in the sum
            iTerm = iTerm + 1;                                              % term index in the polynomial
            for iNarx = 1:nNarx                                             % going throught the full NARX batch
                term(iNarx,iTerm) = regressor(x_narx(:,iNarx),...
                                              indeces{iLambda}(k,:));       % compute the regressor (numeric)
            end
            symb_term{iTerm} = a(indeces{iLambda}(k,:));                    % dictionary of regressors (symbolic)
        end    
    end

end
disp('Dictionary complete')
fileName = ['dict_',dataset,num2str(iFile),'.mat'];
save(fileName, 'term','symb_term','x_narx','y_narx','nNarx','-v7.3');
clear term0 x_narx y_narx nNarx
end                                                                         % end loop over files
nTerms = iTerm;                                                             % total number of regressors in the polynomial
dict_terms = [1:nTerms];                                                    % dictionary of all terms
fileMeta = ['Meta_',dataset];
save(fileMeta, 'nTerms','dict_terms','indeces','lambda','n_y','n_u','K');   % save metadata