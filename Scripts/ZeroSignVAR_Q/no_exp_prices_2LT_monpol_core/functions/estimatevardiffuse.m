function models = estimatevardiffuse(y,p,cons,bVARdraws,trend)

[t,k] = size(y);

x = lagmatrix(y,1:p);

if cons == 1
    x = [ones(length(x),1), x];
end
if trend == 1
    x = [(1:length(x))', x];
end

y(any(isnan(x),2),:) = [];
x(any(isnan(x),2),:) = [];
x(any(isnan(y),2),:) = [];
y(any(isnan(y),2),:) = [];
    
bar = x\y;
res = y - x*bar;
vc = cov(res);

nu=t;
N=x'*x;
B=bar;
S=vc;

for ii=1:bVARdraws
    % Draw variances from the posterior
    Sigma=wishrnd(inv(S)/nu,nu);
    Sigma=inv(Sigma);

    % Draw VAR parameters from the posterior
    vecB=mvnrnd(reshape(B,k*(cons+trend+k*p),1),kron(Sigma,inv(N)),1)';
    Bet=reshape(vecB,cons+trend+k*p,k);

    
    
    res = y - x*Bet;

    if cons == 1 && trend == 1
        A = Bet(3:end,:)';
    elseif cons == 1 || trend == 1
        A = Bet(2:end,:)';
    else
        A = Bet';
    end

    models(ii).A = A;
    % ??? Canova uses the OLS Sigma for each model!
    models(ii).vc = Sigma;
    models(ii).nLags = p;
    models(ii).res = res;
    models(ii).orthP = chol(Sigma,'lower');
end
