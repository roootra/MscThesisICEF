% Sign Restrictions Search
%
% -----------------------------------------------------------------------
% Explanation:
% ===========
%	
% -----------------------------------------------------------------------

function [ identifiedModels , nModelsFound, nModelsChecked ] = signrestrictionsearch(models, S, SN, Z, nVars, hasZeroRestrictions, nModelsFound, nModelsChecked, opt)

% initialize progress bar
totalruns=0; ii=0; jj=0; h = waitbar(0,sprintf('model draw: %d, transformation: %d, models found: %d',ii,jj,nModelsFound),'CreateCancelBtn','setappdata(gcbf,''Cancel'',1)');
% BUG: waitbar does not work with keepAllValid == 0


% Determine first the search order: we search recursively starting always 
% with the structural shock that reveals the highest number of sign 
% restrictions; shocks are only searched among strucural shocks with the 
% same amount of zero restrictions
searchOrder=zeros(nVars);
searchNumber=zeros(nVars,1);
for kk = 1:nVars
    if any(any(Z(:,:,kk)))
    % 2. If yes, search only among structural shocks with
    % exactly the same order of zeros
        sameZeroShocks=min(repmat(any(Z(:,:,kk)),nVars,1)==reshape(any(Z(:,:,:)),[],nVars)',[],2);
    % Count how many sign restrictions each shock has;
    % start to check structural shocks with the hightest
    % amount of sign restrictions
        signShocks=sum(reshape(any(S(:,:,:)),[],nVars))';
        [~,tempOrder]=sort((1+signShocks).*sameZeroShocks,'descend');
        searchOrder(1:size(find(sameZeroShocks==1),1),kk)=tempOrder(1:size(find(sameZeroShocks==1),1));
        searchNumber(kk)=size(find(sameZeroShocks==1),1);
    % 3. If no, search only among structural shockw without
    % zero restrictions
    else
        noZeroShocks=min(reshape(any(Z(:,:,:)),[],nVars)'==0,[],2);
    % Order sign restriction shocks; always start with the
    % shocks with the highest numboer of sign restirctions
        signShocks=sum(reshape(any(S(:,:,:)),[],nVars))';
        [~,tempOrder]=sort((1+signShocks).*noZeroShocks,'descend');
        searchOrder(1:size(find(noZeroShocks==1),1),kk)=tempOrder(1:size(find(noZeroShocks==1),1));
        searchNumber(kk)=size(find(noZeroShocks==1),1);
    end
end

% These two loops search for models fulfilling the sign restrictions.
% The outer loop draws a model from the posterior or goes through all
% models in the posterior (depending on the option chosen).
if opt.drawFromPosterior == 1
    nModels = opt.nModelDraws;
else
    nModels = size(models,2);
end
% The inner loop draws orthogonal matrices to transform the model and
% checks the zero / sign restrictions.
for ii = 1:nModels

    if opt.drawFromPosterior == 1
        % If the option to draw a model from the posterior was selected
        % draw one of the models using the uniform distribution.
        modelPick = randi(size(models,2));
    else
        % If the option to loop through all models was selected simply
        % use the current position in the loop.
        modelPick = ii;
    end
    
    % Pick out the model according to modelPick.
	currentModel = models(modelPick);
    
    % Initialize the while-loop (see below).
    jj = 1;	
    while jj <= opt.nTransformationsPerDraw
            % Update the progress bar to reflect the current transformation
            % and the number of models already found.
            if opt.keepAllValid
                waitbar(totalruns/(nModels*opt.nTransformationsPerDraw),h,sprintf('model draw: %d, transformation: %d, models found: %d',ii,jj,nModelsFound));
                if getappdata(h, 'Cancel'); break; end
            else
                waitbar(ii/nModels,h,sprintf('model draw: %d, transformation: %d, models found: %d',ii,jj,nModelsFound));
                if getappdata(h, 'Cancel'); break; end
            end
            % Copy the model to be transformed to another instance. This is
            % necessary to go through different transformations of the same
            % model.
            transformedModel = currentModel;
            % Initialize the algorithm by drawing a random normal kxk 
            % matrix. 
            X = randn(nVars);
            % Check whether there are any zero restrictions.
            if hasZeroRestrictions              
                % Initialize the Q matrix.
                Q = zeros(nVars);                
            else              
                % Compute the QR decomposition of the random matrix X
                [Q,r]=qr(X);
                for ll=1:nVars
                    if r(ll,ll)<0
                        Q(:,ll)=-Q(:,ll); % reflection matrix
                    end
                end
            end            
            
            % Go through the shocks to be identified (maximum would be
            % equal to the number of variables = k). isSignsFulfilled will
            % be turned 0 if a restriction is not fulfilled and the loop
            % ends.
            % kk = orthogonal shock
			kk = 1; isSignsFulfilled = 1;
            % Generate identifier variable whether one of the structural 
            % shocks is already identified
            isShockIdent=zeros(nVars,1);
            % Initialize new Q matrix (needed to updated order of
            % orthogonal shocks)
            newQ = zeros(nVars);
            % Initialize number of identified shocks
            isNrShockIdent = 0;

            while kk <= nVars && isSignsFulfilled == 1
                if hasZeroRestrictions                  
                    % Construct R. The max()-part ensures that in the first
                    % run Q(kk-1) is the zero vector (as defined in Arias et 
					% al.) taken from the initialization of Q.
                    if any(any(Z(:,:,kk)))
                        % this shock has zero restrictions.
                        % stack them above the Q-matrix to generate R.
                        ZZ = Z(any(Z(:,:,kk),2),:,kk);
                        R = [ ZZ*transformedModel.fA ; Q(:,1:max(1,kk-1))'];
                    else
                        % this shock does not have zero restrictions.
                        % use only the Q-matrix for R.
                        R = Q(:,1:max(1,kk-1))';
                    end
                    % Find a matrix N whose columns form an orthonormal 
                    % basis for the null space of R. 
                    N = null(R);
                    if isempty(N)
                        error('Unable to find a matrix whose columns form an orthonormal basis for the null space of R.')
                    end
                    % basis for the null space of R.
                    % Find a column Q(:,kk) that is the projection of 
                    % X(:,kk) onto the null space of R normalized to be of 
                    % unit length (orthonormal). We then have 
                    % norm(Q:,kk) == 1 and R*Q(:,kk) == 0.
                    Q(:,kk) = N*N'*X(:,kk)/norm(N*N'*X(:,kk));
                    if ~opt.checkNegativeQ && Q(kk,kk)<0
                        Q(:,kk) = Q(:,kk)*(-1);
                    end
                end

                % Check sign restrictions or zero restrictions only if
                % debug-mode is switched off
                if opt.isNoTransform
                    newQ = eye(nVars); 
                else                                    
                    for strkk = searchOrder(1:searchNumber(kk),kk)'
                        
                        % obtain random orthogonal fA matrix
                        orthfA = transformedModel.fA*Q(:,kk);
                        
                        % Cumulate IRFs if applicable
                        if max(opt.isIRFcum) > 0
                            tempIRF = reshape(orthfA,nVars,opt.nMaxSignHorizon+1);
                            for cc = opt.isIRFcum
                                tempIRF(cc,:) = cumsum(tempIRF(cc,:),2);
                            end
                            orthfA = reshape(tempIRF,[],1);
                        end
                        
                        % If shock has more than one sign restriction
                        % (otherweise it is jsut a normalization
                        % restriction), check whether sign restrictions
                        % are fulfilled (zero restrictions are fulfilled
                        % by construction)
                        if sum(any(S(:,:,strkk))) > 1 && min(S(:,:,strkk)*orthfA-SN(:,:,strkk)*ones(size(S(:,:,strkk),2),1)) >= 0
                            % sign restrictions are fulfilled --> check
                            % whether sign restrictons were already
                            % fulfilled before; update identifier variable 
                            isShockIdent(strkk) = isShockIdent(strkk)+1;
                            if isShockIdent(strkk)==1
                                % First time sign restrictions are 
                                % fulfilled --> update new Q matrix and
                                % break for loop
                                newQ(:,strkk) = Q(:,kk);
                                isNrShockIdent = isNrShockIdent+1;
                            end
                            % Sign restrictions were already fulfilled
                            % by other orthogonal shock kk --> do not
                            % update Q
                            % if sign restrictions were fulfilled the first
                            % time, break stops for-loop and we go to next
                            % orthogonal shock; if sign restictions were
                            % already fulfilled before, break starts search
                            % of next model (since isNrShockIdent is not 
                            % updated)
                            break
                        elseif opt.checkNegativeQ && sum(any(S(:,:,strkk))) > 1 && max(S(:,:,strkk)*orthfA+SN(:,:,strkk)*ones(size(S(:,:,strkk),2),1)) <= 0
                            % Check opposite sign restrictions
                            isShockIdent(strkk) = isShockIdent(strkk)+1;
                            if isShockIdent(strkk)==1
                                newQ(:,strkk) = -Q(:,kk);
                                isNrShockIdent = isNrShockIdent+1;
                            end
                            break
                        % If shock has no sign restrictions or only one
                        % sign restriction (i.e. a normalization
                        % restriction) it is categorized as a residual
                        % shock. Residual shocks are ordered randomly
                        % and are not exactly assigned (meaning that
                        % residual shock 1 may also be consistent with
                        % normalization of residual shock 2 and vice versa)
                        elseif sum(any(S(:,:,strkk))) <= 1
                            % If shock fulfills restrictions the first
                            % time order Q matrix accordingly, otherwise go
                            % on with loop
                            isShockIdent(strkk) = isShockIdent(strkk)+1;
                            % If opposit sign of normalization is fulfilled
                            % multiply Q-column with -1
                            tempSign=1;
                            if max(S(:,:,strkk)*orthfA+SN(:,:,strkk)*ones(size(S(:,:,strkk),2),1)) <= 0
                                tempSign=-1;
                            end
                            if isShockIdent(strkk)==1
                                newQ(:,strkk) = Q(:,kk)*tempSign;
                                isNrShockIdent = isNrShockIdent+1;
                                break
                            end
                        end
                    end
                    if isNrShockIdent < kk
                        % The sign restrictions are not fulfilled or fulfilled 
                        % twice for one structural shock. Jump to the next 
                        % model.
                        isSignsFulfilled = 0;
                        break
                    end
                end
                
                if kk == nVars
					% All signs are fulfilled.
                    % Increase the count of models found.
					nModelsFound = nModelsFound + 1;
                    % Replace original Q with differently ordered Q matrix
                    Q = newQ;
                    % Transform the decomposition of the VC matrix from 
                    % PP' = vc to PQQ'P' = vc.
					transformedModel.orthP = currentModel.orthP*Q;
                    % Calculate the impulse responeses based on the
                    % transformed decomposition.
                    if opt.isHistDecompTable || opt.isHistDecompPlots
                        [ transformedModel.irf , ~ , transformedModel.phiStruct ] = calculateirf(transformedModel,size(currentModel.res,1));
                    else
                        [ transformedModel.irf , ~ , transformedModel.phiStruct ] = calculateirf(transformedModel,opt.nImpulseHorizon);
                    end
                    
                    % Cumulate IRFs if applicable
                    if max(opt.isIRFcum) > 0
                        for cc = opt.isIRFcum
                            transformedModel.irf(:,(cc-1)*nVars+1:cc*nVars) = cumsum(transformedModel.irf(:,(cc-1)*nVars+1:cc*nVars),1);
                        end
                    end          
                    
                    % Delete the stacked impulse responses as they are only
                    % required within the SR algorithm.
                    transformedModel = rmfield(transformedModel,'fA');
                    % Put the model into the array of models fulfilling the
                    % sign restrictions.
					identifiedModels(nModelsFound) = transformedModel;
                end
                kk = kk + 1;
            end
     
            nModelsChecked = nModelsChecked + 1;

            if opt.keepAllValid == 0 && isSignsFulfilled == 1
				% Only keep the first transformation that fits the sign 
                % restrictions for this model.
                break
            end
            
            if nModelsFound >= opt.nMatches
                % Stop while-loop if the requested number of models has been found
                break
            end
        % Totalruns is only needed for the progress bar.
        jj = jj + 1; totalruns = totalruns + 1;
    end
    
    jj = 1;
    % Update the progress bar with the new model draw and the number of 
    % models already found.
    if opt.keepAllValid
        waitbar(totalruns/(nModels*opt.nTransformationsPerDraw),h,sprintf('model draw: %d, transformation: %d, models found: %d',ii,jj,nModelsFound));
        if getappdata(h, 'Cancel'); break; end
    else
        waitbar(ii/nModels,h,sprintf('model draw: %d, transformation: %d, models found: %d',ii,jj,nModelsFound));
        if getappdata(h, 'Cancel'); break; end
    end
    
    if nModelsFound >= opt.nMatches
        % Stop for-loop if the requested number of models has been found
        
        break
    end
    
end
delete(h);

if ~nModelsFound
    error('No model fitting the restrictions was found.')
end

end