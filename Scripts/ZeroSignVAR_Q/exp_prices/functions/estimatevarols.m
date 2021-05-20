function model = estimatevarols(y,nLags,useConstant,useTrend)
% ESTIMATEVAROLS calculates the OLS coefficients of the reduced form VAR.
% The observation period reduces according to the lag-length. If the 
% data-matrix y includes missing values the estimation-function forces 
% the data-matrix to be a balenced panel, by producing first the lag-matrix
% x and subsequently deleting all rows with missing values from both 
% matrices.
%
% The function returns a structure containing the following elements:
%   - model.A ... the reduced form parameter matrix.
%   - model.vc ... the variance-covariance matrix.
%   - model.nLags ... the number of lags in the model.
%   - model.res ... the residuals of the model.
%   - model.orthP ... an orthogonal matrix P with P*P' = vc.

    x = lagmatrix(y,1:nLags);

    if useConstant == 1
        x = [ones(length(x),1), x];
    end
    if useTrend == 1
        x = [(1:length(x))', x];
    end
    
    y(any(isnan(x),2),:) = [];
    x(any(isnan(x),2),:) = [];
    x(any(isnan(y),2),:) = [];
    y(any(isnan(y),2),:) = [];
    
    bar = x\y;    

    % Residuen ausrechnen:
    res = y - x*bar;

    % und die Kovarianzmatrix (nach EViews Methode, d.h. nicht ML sondern dof corrected):
    vc = res'*res/(size(y,1)-size(x,2));

    % create A matrix, which has the following structure
    %   e.g. for a VAR with 3 variables and 2 lags:
    %   Y1t=a11,1*Y1t-1+a11,2*Y1t-2+a12,1*Y2t-1+...+a13,2*Y3t-2+u1  
    %   Y2t=a21,1*Y1t-1+a21,2*Y1t-2+a22,1*Y2t-1+...+a23,2*Y3t-2+u1 
    %   Y3t=a31,1*Y1t-1+a31,2*Y1t-2+a32,1*Y2t-1+...+a33,2*Y3t-2+u1 
    %   
    %       a11,1   a12,1   a13,1   |   a11,2   a12,2   ... |
    %   A = a21,1   ...             |   a21,1               |
    %       a31,1           a33,1   |   a31,1           ... |
    %
    % if there is a constant/trend it has to be dropped
    if useConstant == 1 && useTrend == 1
        A = bar(3:end,:)';
    elseif useConstant == 1 || useTrend == 1
        A = bar(2:end,:)';
    else
        A = bar';
    end

    % standard errors and t-statistics of coefficients
    
    seA = kron(inv(x'*x), vc);
    seA = sqrt(diag(seA));
    
    vecA = reshape(bar',size(bar,1)*size(bar,2),1);
    
    tstatA = zeros(size(bar,1)*size(bar,2),1);
    for ii = 1:size(bar,1)*size(bar,2)
        tstatA(ii) = vecA(ii)/seA(ii);
    end
    
    tstatA = reshape(tstatA,size(bar,2),size(bar,1));
    
    if useConstant == 1
        tstatA = tstatA(:,2:end);
    end
    
    seA = reshape(seA,size(bar,2),size(bar,1));
    
    if useConstant == 1
        seA = seA(:,2:end);
    end
    

    model.A = A;
    model.Ase = seA;
    model.Atstat = tstatA;
    model.vc = vc;
    model.nLags = nLags;
    model.res = res;
    model.orthP = chol(vc,'lower');
	
end