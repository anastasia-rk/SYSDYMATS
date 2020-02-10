function [result] = BIC(r_n,nData,nTerms)
% Input
%   nData  - length of the dataset;
%   nTerms - number of model terms;
%   r_n    - associated model residual;
% Output
%   result - BIC of the profiled likelhihood
%   BIC - Bayesian Information Criteria
    norm_r = norm(r_n);
    result = log(nData)*nTerms + nData*log(norm_r/nData);
end