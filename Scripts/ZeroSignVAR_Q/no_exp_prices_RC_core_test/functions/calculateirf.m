function [ irf , fA , phiStruct ] = calculateirf(model,q)
% CALCULATEIRF  Calculate the impulse responses and the stacked impulse
% responses for a VAR model.
%
%   The impulse responses are ordered as follows:
%   - Col 1: Response of variable 1 to variable 1
%   - Col 2: Response of variable 1 to variable 2
%   ...
%   - Col (i-1)*k+j: Response of variable i to variable j, where
%   i,j(el){1,...k} and k = {number of variables in the VAR}
%
%   The rows correspond to the horizon of the response.

    p = model.nLags;
    P = model.orthP;
    A = model.A;
    k = size(A,1);
    
	% compute Moving Average coefficients
	phi = eye(k);
	for ii = 1:q
		sum = 0;
		for jj = 1:min([p,ii])
			sum = sum + phi(:,1+(ii-jj)*k:1+(ii-jj)*k+(k-1))*A(:,1+(ii-(ii-jj)-1)*k:1+(ii-(ii-jj)-1)*k+(k-1));
		end
		phi = [phi, sum];
	end

	% transform Phi
	phis = zeros(k,(q+1)*k);
	for ii = 1:(q+1)
		phis(:,(ii-1)*k+1:(ii-1)*k+k) = phi(:,(ii-1)*k+1:(ii-1)*k+k)*P; % 1 std. impulse 
		% phis(:,(ii-1)*k+1:(ii-1)*k+k) = phi(:,(ii-1)*k+1:(ii-1)*k+k)*W; % unit impulse
	end

	impulse_chol = zeros(q+1,k^2);
	impulse_chol_unit = zeros(q+1,k^2);
	for ii = 1:k
		impulse_chol(:,ii:k:end) = phis(:,ii:k:end)';
		% impulse_chol_unit(:,ii:k:end) = phis(:,ii:k:end)'*(1/D(ii,ii));
    end
	
    ident = eye(size(phis,2));
    fA = zeros(k*(q+1),k);
    for ii = 1:(q+1)
        fA(1+(ii-1)*k:1+(ii-1)*k+(k-1),:) = phis*ident(1+(ii-1)*k:1+(ii-1)*k+(k-1),:)';
    end
    
	irf = impulse_chol;

    phiStruct = phis;
end