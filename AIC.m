function [result] = AIC(r_n,nData,nTerms)
% Input
%   nData  - length of the dataset;
%   nTerms - number of model terms;
%   r_n    - associated model residual;
% Output
%   result - AIC of the profiled likelhihood
%   AIC    - Akaike Information Criteria
    norm_r = norm(r_n);
    result = 2*nTerms + nData*log(norm_r/nData);
end