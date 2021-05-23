function [ output , ctmModel ] = generateOutput( y, identifiedModels, S, Z, nModelsFound, nVars, opt )

    % Stack the impulse responses into one matrix for further calculations.
    % This matrix has the following structure:
    % Columns 1:nVars^2 contain the contemporaneous impulse responses.
    % The next nVars^2 columns contain the impulse responses at horizon 1.
    % ...
    % Each row contains the impulse responses from one identified model.
    temp.irf = cat(2,[identifiedModels.irf]);
    temp.irfStack = zeros(nModelsFound,(nVars^2)*(opt.nImpulseHorizon+1));
    for ii = 1:(opt.nImpulseHorizon+1)
        for jj = 1:nVars^2
            temp.irfStack(:,(ii-1)*(nVars^2)+jj) = temp.irf(ii,jj:nVars^2:end)';
        end
    end

    % Calculate the pointwise median and the pointwise 68 percentile band
    % around the median.
    temp.irfMedian = median(temp.irfStack,1);
    temp.irf95 = prctile(temp.irfStack,95,1);
    temp.irf84 = prctile(temp.irfStack,84,1);
    temp.irf16 = prctile(temp.irfStack,16,1);
    temp.irf05 = prctile(temp.irfStack,5,1);

    % Generate selection matrix to multiply sqDistanceFromMedian with.
    selectCTM = eye((opt.nCTMHorizon+1)*(nVars^2));
    if opt.narrowCTMSearch == 1
        for ii = 1:nVars
            if ~(sum(any(S(:,:,ii)))>1 || any(any(Z(:,:,ii))))
                selectCTM(:,ii:nVars:end) = 0;
            end
        end   
    end

    % For each model take the impulse responses up to the horizon specified
    % in the options and calculate the difference to the pointwise median.
    % Divide this difference by the pointwise standard deviation across
    % the models. Square the result (alternatively the absolute values could
    % be used), calculate sums over the models and select the model with 
    % the smallest sum.
    
    error = (temp.irfStack(:,1:(opt.nCTMHorizon+1)*(nVars^2)) - repmat(temp.irfMedian(:,1:(opt.nCTMHorizon+1)*(nVars^2)),nModelsFound,1))./repmat(std(temp.irfStack(:,1:(opt.nCTMHorizon+1)*nVars^2)),nModelsFound,1);

    % if 0 pointwise standard deviation are present (e.g. by
    % Cholesky) replace NaN with 0 after the division
    
    error(isnan(error)) = 0 ;
    
    sqDistanceFromMedian = power(error,2)*selectCTM;
    totalSqDistanceFromMedian = sum(sqDistanceFromMedian')';
    [~,closestImpulse] = min(totalSqDistanceFromMedian);
    ctmModel = identifiedModels(closestImpulse);
    output.irfCTM = ctmModel.irf(1:opt.nImpulseHorizon+1,:);

    % Bring the calculated percentiles back to the regular form
    for ii = 1:nVars^2
        output.irfMedian(:,ii) = temp.irfMedian(ii:nVars^2:end)';
        output.irf95(:,ii) = temp.irf95(ii:nVars^2:end)';
        output.irf84(:,ii) = temp.irf84(ii:nVars^2:end)';
        output.irf16(:,ii) = temp.irf16(ii:nVars^2:end)';
        output.irf05(:,ii) = temp.irf05(ii:nVars^2:end)';
    end
    
    %% FEVD
    if opt.isFevdTable == 1    
        
        identifiedModels = calculateFEVD ( identifiedModels , nVars , opt );
        
        temp.fevd = cat(2,[identifiedModels.fevd]);
        temp.fevdStack = zeros(nModelsFound,(nVars^2)*(opt.nImpulseHorizon+1));
        for ii = 1:(opt.nImpulseHorizon+1)
            for jj = 1:nVars^2
                temp.fevdStack(:,(ii-1)*(nVars^2)+jj) = temp.fevd(ii,jj:nVars^2:end)';
            end
        end

        % Calculate the pointwise median and the pointwise 68 percentile band
        % around the median.
        temp.fevdMedian = median(temp.fevdStack,1);
        temp.fevd95 = prctile(temp.fevdStack,95,1);
        temp.fevd84 = prctile(temp.fevdStack,84,1);
        temp.fevd16 = prctile(temp.fevdStack,16,1);
        temp.fevd05 = prctile(temp.fevdStack,5,1);

        % Bring the calculated percentiles back to the regular form
        for ii = 1:nVars^2
            output.fevdMedian(:,ii) = temp.fevdMedian(ii:nVars^2:end)';
            output.fevd95(:,ii) = temp.fevd95(ii:nVars^2:end)';
            output.fevd84(:,ii) = temp.fevd84(ii:nVars^2:end)';
            output.fevd16(:,ii) = temp.fevd16(ii:nVars^2:end)';
            output.fevd05(:,ii) = temp.fevd05(ii:nVars^2:end)';
        end        
         
    end
    
    % CTM - Forecast Error Variance Decomposition
    if opt.isFevdTable == 1
        ctmModel = identifiedModels(closestImpulse);
        output.fevdCTM = ctmModel.fevd;
    elseif opt.isFevdPlots == 1
        eyek = eye(nVars);
        output.ctmFevd = zeros(opt.nImpulseHorizon+1,nVars^2);
        for jj = 1:nVars
            for kk = 1:nVars
                tempSum = 0;
                tempSum2 = 0;
                for ii = 1:(opt.nImpulseHorizon+1)
                    tempSum = tempSum + power(eyek(:,jj)'*ctmModel.phiStruct(:,(ii-1)*nVars+1:(ii-1)*nVars+nVars)*eyek(:,kk),2);
                    tempSum2 = tempSum2 + ctmModel.phiStruct(:,(ii-1)*nVars+1:(ii-1)*nVars+nVars)*ctmModel.phiStruct(:,(ii-1)*nVars+1:(ii-1)*nVars+nVars)';
                    output.ctmFevd(ii,(jj-1)*nVars+1+kk-1) = (tempSum/tempSum2(jj,jj))*100;
                end
            end
        end
    end

    if opt.isStructShockPlots || opt.isStructShockTable
        % Structural shocks
        temp.structShock = cat(2,[identifiedModels.structShocks]);
        temp.structShockStack = zeros(nModelsFound,nVars*size(identifiedModels(1).res,1));
        for ii = 1:nVars
            temp.structShockStack(:,1+(ii-1)*size(temp.structShock,1):ii*size(temp.structShock,1)) = temp.structShock(:,ii:nVars:end)';
        end
        temp.structShockMedian = median(temp.structShockStack,1);
        temp.structShock95 = prctile(temp.structShockStack,95,1);
        temp.structShock84 = prctile(temp.structShockStack,84,1);
        temp.structShock16 = prctile(temp.structShockStack,16,1);
        temp.structShock05 = prctile(temp.structShockStack,5,1);   
        output.structShockMedian = zeros(size(identifiedModels(1).structShocks));
        output.structShock95 = output.structShockMedian;
        output.structShock84 = output.structShockMedian;
        output.structShock16 = output.structShockMedian;
        output.structShock05 = output.structShockMedian;
        for ii = 1:nVars
            output.structShockMedian(:,ii) = temp.structShockMedian(:,1+size(identifiedModels(1).res,1)*(ii-1):size(identifiedModels(1).res,1)*ii)';
            output.structShock95(:,ii) = temp.structShock95(:,1+size(identifiedModels(1).res,1)*(ii-1):size(identifiedModels(1).res,1)*ii)';
            output.structShock84(:,ii) = temp.structShock84(:,1+size(identifiedModels(1).res,1)*(ii-1):size(identifiedModels(1).res,1)*ii)';
            output.structShock16(:,ii) = temp.structShock16(:,1+size(identifiedModels(1).res,1)*(ii-1):size(identifiedModels(1).res,1)*ii)';
            output.structShock05(:,ii) = temp.structShock05(:,1+size(identifiedModels(1).res,1)*(ii-1):size(identifiedModels(1).res,1)*ii)';
        end
        
        x = lagmatrix(y,1:opt.nLags);
        temp2 = any(isnan(x),2);
        temp2(any(isnan(y),2)) = 1;
        temp2 = temp2(opt.nLags+1:end,:);
        if any(temp2)
            structShockMedian = NaN(length(temp2), nVars);
            structShock95 = structShockMedian;
            structShock84 = structShockMedian;
            structShock16 = structShockMedian;
            structShock05 = structShockMedian;
            temp3 = [];
            for ii = 1:length(temp2)
                if temp2(ii)==0
                    structShockMedian(ii,:)=output.structShockMedian(ii-length(temp3),:);
                    structShock95(ii,:)=output.structShock95(ii-length(temp3),:);
                    structShock84(ii,:)=output.structShock84(ii-length(temp3),:);
                    structShock16(ii,:)=output.structShock16(ii-length(temp3),:);
                    structShock05(ii,:)=output.structShock05(ii-length(temp3),:);
                else
                    temp3 = [temp3, 1];
                end
            end      
            output.structShockMedian = structShockMedian;
            output.structShock95 = structShock95;
            output.structShock84 = structShock84;
            output.structShock16 = structShock16;
            output.structShock05 = structShock05;
        end
    end
    
    if opt.isHistDecompTable || opt.isHistDecompPlots
        % Historical Decomposition
        temp.HistDecomp = cat(2,[identifiedModels.HistDecomp]);
        temp.HistDecompStack = zeros(nModelsFound,nVars*nVars*size(identifiedModels(1).res,1));
        for ii = 1:nVars*nVars
            temp.HistDecompStack(:,1+(ii-1)*size(temp.HistDecomp,1):ii*size(temp.HistDecomp,1)) = temp.HistDecomp(:,ii:nVars*nVars:end)';
        end
        temp.HistDecompMedian = median(temp.HistDecompStack,1);
        temp.HistDecomp95 = prctile(temp.HistDecompStack,95,1);
        temp.HistDecomp84 = prctile(temp.HistDecompStack,84,1);
        temp.HistDecomp16 = prctile(temp.HistDecompStack,16,1);
        temp.HistDecomp05 = prctile(temp.HistDecompStack,5,1);   
        output.HistDecompMedian = zeros(size(identifiedModels(1).HistDecomp));
        output.HistDecomp95 = output.HistDecompMedian;
        output.HistDecomp84 = output.HistDecompMedian;
        output.HistDecomp16 = output.HistDecompMedian;
        output.HistDecomp05 = output.HistDecompMedian;
        for ii = 1:nVars*nVars
            output.HistDecompMedian(:,ii) = temp.HistDecompMedian(:,1+size(identifiedModels(1).res,1)*(ii-1):size(identifiedModels(1).res,1)*ii)';
            output.HistDecomp95(:,ii) = temp.HistDecomp95(:,1+size(identifiedModels(1).res,1)*(ii-1):size(identifiedModels(1).res,1)*ii)';
            output.HistDecomp84(:,ii) = temp.HistDecomp84(:,1+size(identifiedModels(1).res,1)*(ii-1):size(identifiedModels(1).res,1)*ii)';
            output.HistDecomp16(:,ii) = temp.HistDecomp16(:,1+size(identifiedModels(1).res,1)*(ii-1):size(identifiedModels(1).res,1)*ii)';
            output.HistDecomp05(:,ii) = temp.HistDecomp05(:,1+size(identifiedModels(1).res,1)*(ii-1):size(identifiedModels(1).res,1)*ii)';
        end
        
        x = lagmatrix(y,1:opt.nLags);
        temp2 = any(isnan(x),2);
        temp2(any(isnan(y),2)) = 1;
        temp2 = temp2(opt.nLags+1:end,:);
        if any(temp2)
            HistDecompMedian = NaN(length(temp2), nVars*nVars);
            HistDecomp95 = HistDecompMedian;
            HistDecomp84 = HistDecompMedian;
            HistDecomp16 = HistDecompMedian;
            HistDecomp05 = HistDecompMedian;
            temp3 = [];
            for ii = 1:length(temp2)
                if temp2(ii)==0
                    HistDecompMedian(ii,:)=output.HistDecompMedian(ii-length(temp3),:);
                    HistDecomp95(ii,:)=output.HistDecomp95(ii-length(temp3),:);
                    HistDecomp84(ii,:)=output.HistDecomp84(ii-length(temp3),:);
                    HistDecomp16(ii,:)=output.HistDecomp16(ii-length(temp3),:);
                    HistDecomp05(ii,:)=output.HistDecomp05(ii-length(temp3),:);
                else
                    temp3 = [temp3, 1];
                end
            end      
            output.HistDecompMedian = HistDecompMedian;
            output.HistDecomp95 = HistDecomp95;
            output.HistDecomp84 = HistDecomp84;
            output.HistDecomp16 = HistDecomp16;
            output.HistDecomp05 = HistDecomp05;
        end

        % Median Base Projection (for plotting histDecomp required)
        temp.BaseProjectMedian = cat(2,[identifiedModels.BaseProject]);
        temp.BaseProjectMedianStack = zeros(nModelsFound,nVars*size(identifiedModels(1).BaseProject,1));
        for ii = 1:nVars
            temp.BaseProjectMedianStack(:,1+(ii-1)*size(temp.BaseProjectMedian,1):ii*size(temp.BaseProjectMedian,1)) = temp.BaseProjectMedian(:,ii:nVars:end)';
        end
        temp.BaseProjectMedian = median(temp.BaseProjectMedianStack,1);   
        output.BaseProjectMedian = zeros(size(identifiedModels(1).BaseProject));
        for ii = 1:nVars
            output.BaseProjectMedian(:,ii) = temp.BaseProjectMedian(:,1+size(identifiedModels(1).BaseProject,1)*(ii-1):size(identifiedModels(1).BaseProject,1)*ii)';
        end
        
        x = lagmatrix(y,1:opt.nLags);
        temp2 = any(isnan(x),2);
        temp2(any(isnan(y),2)) = 1;
        temp2 = temp2(opt.nLags+1:end,:);
        if any(temp2)
            BaseProjectMedian = NaN(length(temp2), nVars);
            temp3 = [];
            for ii = 1:length(temp2)
                if temp2(ii)==0
                    BaseProjectMedian(ii,:)=output.BaseProjectMedian(ii-length(temp3),:);
                else
                    temp3 = [temp3, 1];
                end
            end      
            output.BaseProjectMedian = BaseProjectMedian;
        end
    end    
    
end