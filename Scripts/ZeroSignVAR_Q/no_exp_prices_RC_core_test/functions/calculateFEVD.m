function models = calculateFEVD( models , nVars , opt )
    eyek = eye(nVars);
    for ii = 1:size(models,2)
        models(ii).fevd = zeros(opt.nImpulseHorizon+1,nVars^2);
        for jj = 1:nVars
            for kk = 1:nVars
                tempSum = 0;
                tempSum2 = 0;
                for ll = 1:(opt.nImpulseHorizon+1)
                    tempSum = tempSum + power(eyek(:,jj)'*models(ii).phiStruct(:,(ll-1)*nVars+1:(ll-1)*nVars+nVars)*eyek(:,kk),2);
                    tempSum2 = tempSum2 + models(ii).phiStruct(:,(ll-1)*nVars+1:(ll-1)*nVars+nVars)*models(ii).phiStruct(:,(ll-1)*nVars+1:(ll-1)*nVars+nVars)';
                    models(ii).fevd(ll,(jj-1)*nVars+1+kk-1) = (tempSum/tempSum2(jj,jj))*100;
                end
            end
        end
    end
end