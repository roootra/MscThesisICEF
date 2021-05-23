function models = calculateStructShocks( models, nVars )
% CALCULATESTRUCTSHOCKS Calculates the structural shocks
% for a VAR model.
%
%   Structural shocks = (impact Matrix)^-1 * reduced form residuals
%   et = (PQ)^-1 * ut
%   where, P is the lower triangular Cholesky factor and Q is a orthogonal
%   matrix with QQ'=I.
%
%   The output is saved in the identifiedModel-Structure with the following
%   order:
%   - Col 1: Shock 1
%   - Col 2: Shock 2
%   - Col k: Shock k
%   where k = {number of variables in the VAR}
%
%   The rows correspond to the time period.

for ii = 1:size(models,2)
    models(ii).structShocks = (inv(models(ii).orthP)*models(ii).res')';
end
    
end