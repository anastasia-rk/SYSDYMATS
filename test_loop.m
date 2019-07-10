my_init;

d = 4;                                  % Size of input vector x
lambda = 3;                             % order of polynomial terms
times = 1;                              % time indeces
x = [1:4]';                             % example input vector
a = sym('a',[1 length(x)]);             % associated symbolic vector
% ind{1:lambda} = zeros(1,d);
%% Index combinations in polynomial NARX of order lambda
indeces{1} = [1:d]';
for iLambda=2:lambda                    % iLambda - order of polynomial term
    initial_m = permn(1:d,iLambda);     % get all permutations with repetition
    M = initial_m;                      % M - temp set
    for j=iLambda:-1:2
        ind = find(M(:,j)>=M(:,j-1));   % sort out only increasing indeces
        M = M(ind,:);                   % update set
        clear ind
    end
    indeces{iLambda} = M;               % all index combinations of order iLambda
end

for t = times
    % Computing all polynomial terms
    iTerm = 0;
    for iLambda = 1:lambda
        for k = 1:size(indeces{iLambda},1)
           iTerm = iTerm + 1; 
           term(1,iTerm) = regressor(x(:,t),indeces{iLambda}(k,:));
           symb_term{iTerm} = a(indeces{iLambda}(k,:));
        end
    end
    nTerms = iTerms;    % total number of regressors
end