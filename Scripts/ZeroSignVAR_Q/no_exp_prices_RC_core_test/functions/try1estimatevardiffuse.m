function models = try1estimatevardiffuse(y,nLags,useConstant,nDrawsFromBvar)

    nObs = length(y);
    n = size(y,2);                  % number of variables in the VAR
    k = n*nLags+useConstant;      	% number of coefficients each equation

    x = lagmatrix(y,1:nLags);

    x(any(isnan(x),2),:) = [];
    
    if useConstant == 1
        x = [ones(length(x),1), x];
    end

    y = y((nLags+1):end,:);
    

    bar = x\y;

    % Residuen ausrechnen:
    res = y - x*bar;

    Sigma = cov(res);

    nu = nObs - k*n;
    S=Sigma;
    B=bar;

    % draws from the posterior
    % ??? Canova uses OLS estimate as first draw.
    betaDraw = B;
    if useConstant == 1
        A = betaDraw(2:end,:)';
    else
        A = betaDraw';
    end        

    models(1).A = A;
    % ??? Canova uses the OLS Sigma for each model!
    models(1).vc = Sigma;
    models(1).nLags = nLags;
    models(1).res = res;
    models(1).orthP = chol(Sigma,'lower');
    for ii=2:nDrawsFromBvar
        % draw from an IW for Sigma
        PS = real(sqrtm(inv(nu*S))); 
        u = randn(n,nu);
        sigmaDraw(:,:) = inv(PS*u*u'*PS');
        % draw from N for beta, conditional on the draw for Sigma 
        N = kron(inv(sigmaDraw),x'*x);
        P=chol(inv(N));
        BBBet=reshape(B,n*(n*nLags+useConstant),1)+ P*randn(n*(n*nLags+useConstant),1);
        betaDraw = reshape(BBBet,n*nLags+1,n);
        
        % Residuen ausrechnen:
        res = y - x*betaDraw;
        
        if useConstant == 1
            A = betaDraw(2:end,:)';
        else
            A = betaDraw';
        end        
        
        models(ii).A = A;
        % ??? Canova uses the OLS Sigma for each model!
        models(ii).vc = sigmaDraw;
        models(ii).nLags = nLags;
        models(ii).res = res;
        models(ii).orthP = chol(sigmaDraw,'lower');
    end