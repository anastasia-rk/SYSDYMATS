my_init;
%% Set up and simulate the model

narx = @(u,y,n_u,n_y,Theta) Theta(1:n_u)*u(:,end-n_u+1:end)' + Theta(n_u+1:n_u+n_y)*y(:,end-n_y+1:end)';

theta = ones(1,10);
u = [1:10];
y = [1:10];
n_u = 3;
n_y = 2;
%% Create batch of regressors
d = n_u + n_y;                                                              % size of input vector x
lambda = 4;                                                                 % order of polynomial
times = [1];                                                                % batch over times
for t = times
    x(:,t) = [u(t-n_u:t);y(t-n_y:t)];                                       % create input vector 
    % Computing all polynomial terms
    iTerm = 0;
    for iLambda = 1:lambda
        for k = 1:size(indeces{iLambda},1)
           iTerm = iTerm + 1;  
           term(1,iTerm) = regressor(x(:,t),indeces{iLambda}(k,:));         % compute the regressor
           symb_term{iTerm} = a(indeces{iLambda}(k,:));                     % dictionary of regressors
        end
    end
    nTerms = iTerms;                                                        % total number of regressors
end
