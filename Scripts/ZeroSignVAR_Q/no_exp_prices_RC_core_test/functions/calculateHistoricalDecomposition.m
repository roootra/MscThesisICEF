function models = calculateHistoricalDecomposition( models, nVars, y )
% CALCULATEHISTORICALDECOMPOSITION Calculates the historical decomposition
% for a VAR model.
%
%   The output is ordered as follows:
%   - Col 1: Contribution of shock 1 to variable 1
%   - Col 2: Contribution of shock 2 to variable 1
%   ...
%   - Col (i-1)*k+j: Response of variable i to variable j, where
%   i,j(el){1,...k} and k = {number of variables in the VAR}
%
%   The rows correspond to the time period.



% stack the phis and shocks of all identified models in 3-dimensional array
%
% order of tempPhiStack: phi0 to phi max horizons with dimension
% (nVars x nVars)in columns; models are stacked page-wise; hence overall
% dimension (nVars x nVars*Horizont x nIdentifiedModels).
%
% order of phi per model (nVars*(nVars*periods+1):
% first column shows the effects of shock 1 to variable 1 to nVars
% at h=0
% second column shows effects of shocks 2 to variable 1 to nVars
% at h=0
% phis for each model are vertically stacked (size(models,2)*nVars*(nVars*periods)

% order of tempShockStack: shocks are in columns, shocks in t=1..T in rows;
% and models are stacked page-wise; hence overall dimension 
% (Time x nVars x nIdentifiedModels).
tempPhiStack = cat(3,models.phiStruct);
tempShockStack = cat(3,models.structShocks);

HistDecomp = zeros(size(models(1).res,1),nVars^2,size(models,2));
for jj = 1:size(models(1).res,1) % jj from first to last observation period
    % the loop calculates the HD for the whole observation period; however,
    % starting and end point can actually be chosen arbitrarily. values
    % are the same (maybe interesting if calculation for whole period may
    % be too slow).
    
    horizon = jj-1; % jj-1 for whole horizon; specify if only shorter 
    % horizon should be considered.
    hh = min(horizon, jj-1);
    % select shocks from jj-horizon to jj and transpose the matrix
    tempShocks = permute(tempShockStack(jj-hh:jj,:,:),[2 1 3]);
    % reorder such that final matrix goes from t=jj to t=jj-1 to t=1
    tempShocks = flipdim(tempShocks,2);
    % create a row vector, which can then be duplicated
    tempShocks = reshape(tempShocks, 1 , nVars*(hh+1), []);
    % create repeat entries for element-wise multiplication
    tempShocks = repmat(tempShocks, [nVars 1 1]);
    
    % multiply shocks with phis
    % same order as phis:
    % first column   - effect of shock 1 at t=jj on variables 1 to nVars
    % second column  - effect of shock 2 at t=jj on variables 1 to nVars
    % nVars+1 column - effect of shock 1 at t=jj-1 on variables 1 to nVars
    tempMatrix = tempPhiStack(:,1:nVars*(hh+1),:).*tempShocks;
    
    % Historical decomposition is the sum across jj periods
    % reshape for summation across the periods jj and reshape back that the
    % output can be similarly ordered as the IRFs
    tempMatrix = reshape(tempMatrix, nVars^2, (hh+1), []);
    tempSum = sum(tempMatrix,2);
    tempSum = reshape(tempSum, nVars, nVars, []);

    % save HistDecomp of period jj
    % change order according to IRFs:
    % first column   - effect of shock 1 on variable 1
    % second column  - effect of shock 2 on variable 1
    % nVars+1 column - effect of shock 1 on variable 2
    HistDecomp(jj,:,:) = reshape(permute(tempSum,[2 1 3]), 1, nVars^2, []);
        
end


% Reshape HistDecomp and match with model structure; generate base
% projection of each model within the same loop (loop works relatively
% fast)
for ii = 1:size(models,2)
    models(ii).HistDecomp=HistDecomp(:,:,ii);

    % Generates the base projection of the SVAR
    for kk = 1:nVars
        models(ii).BaseProject(:,kk) = y(models(1).nLags+1:end,kk);
        models(ii).BaseProject(:,kk) = models(ii).BaseProject(:,kk)-sum(models(ii).HistDecomp(:,(kk-1)*nVars+1:kk*nVars),2);
    end
end
    
end