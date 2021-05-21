% SIGNRESTRICTIONS
%
% opt.
%   isFevdTable:        Calculate the FEVD for all models (1 = Yes, 0 = No) 
%                       (Default = 0)
%   narrowCTMSearch:    When looking for the closest-to-median impulse 
%                       response, use only the impulse responses associated 
%                       with a shock with zero or sign restrictions. 
%                       (1 = Yes, 0 = No, use all impulse responses)
%                       (Default = 0)
%   drawFromPosterior:  draw from posterior distribution when working through
%                       sign restrictions (1 = Yes, 0 = work through all draws)
%                       (Default = 1)
%   checkNegativeQ:     if a sign restriction is fullfilled for the negative 
%                       of a Q column transform this column of Q to it's negative.
%                       (1 = Yes, 0 = No) (Default = 1)
%   isNoTransform:      No transformations (for debugging) (1 = Yes, 0 = No)
%                       (Default = 0)
%   isImpulsePlots:     This option plots the impulse responses (default = 0).
%   isPlotCTM:          If isImpulsePlots is turned on, the CTM will be plotted by default.
%                       By setting isPlotCTM to 0, the impulse plots are provided with and
%                       without the CTM.

% Create folders if they don't exist already.
if ~exist(opt.modelPath,'dir')
  mkdir(opt.modelPath);
end

if ~exist(strcat(opt.modelPath, opt.modelName),'dir')
  mkdir(strcat(opt.modelPath, opt.modelName));
end

if ~exist(strcat(opt.modelPath, opt.modelName, '/imgs'),'dir')
  mkdir(strcat(opt.modelPath, opt.modelName, '/imgs'));
end

if ~exist(strcat(opt.modelPath, opt.modelName, '/tables'),'dir')
  mkdir(strcat(opt.modelPath, opt.modelName, '/tables'));
end

% assign standard values to variables that have not been assigned
if ~exist('SN','var')
    SN = zeros(size(S));
end

% assign standard values to options that have not been assigned
if ~isfield(opt, 'runNumber')
    opt.runNumber = 1;
end
if ~isfield(opt, 'hasConstant')
    opt.hasConstant = 1;
end
if ~isfield(opt, 'hasTrend')
    opt.hasTrend = 0;
end
if ~isfield(opt,'nImpulseHorizon')
    opt.nImpulseHorizon = 12;
end
if ~isfield(opt, 'nDrawsFromBvar')
    opt.nDrawsFromBvar = 1000;
end
if ~isfield(opt, 'nTransformationsPerDraw')
    opt.nTransformationsPerDraw = 1000;
end
if ~isfield(opt,'isFevdTable')
    opt.isFevdTable = 0;
end
if ~isfield(opt,'isFevdTexTable')
    opt.isFevdTexTable = 0;
end
if ~isfield(opt,'nCTMHorizon')
    opt.nCTMHorizon = opt.nImpulseHorizon;
end
if ~isfield(opt,'keepAllValid')
    opt.keepAllValid = 0;
end
if ~isfield(opt,'narrowCTMSearch')
    opt.narrowCTMSearch = 0;
end
if ~isfield(opt,'drawFromPosterior')
    opt.drawFromPosterior = 0;
end
if ~isfield(opt,'nModelDraws')
    opt.nModelDraws = 1000;
end
if ~isfield(opt,'checkNegativeQ')
    opt.checkNegativeQ = 1;
end
if ~isfield(opt,'isNoTransform')
    opt.isNoTransform = 0;
end
if ~isfield(opt,'isIRFtable')
    opt.isIRFtable = 0;
end
if ~isfield(opt,'isSaveResults')
    opt.isSaveResults = 0;
end
if ~isfield(opt,'nImpulseHorizon')
    opt.nImpulseHorizon = 12;
end
if ~isfield(opt,'isImpulsePlots')
    opt.isImpulsePlots = 1;
end
if ~isfield(opt,'isPlotCTM')
    opt.isPlotCTM = 0;
end
if ~isfield(opt,'nPlotRandomModels')
    opt.nPlotRandomModels = 0;
end
if ~isfield(opt,'isStructShockTable')
    opt.isStructShockTable = 0;
end
if ~isfield(opt,'nMatches')
    opt.nMatches = 1000000;
end
if ~isfield(opt,'isHistDecompTable')
    opt.isHistDecompTable = 0;
end
if ~isfield(opt,'isHistDecompPlots')
    opt.isHistDecompPlots = 0;
end
if ~isfield(opt,'isFevdPlots')
    opt.isFevdPlots = 0;
end
if ~isfield(opt,'isStructShockPlots')
    opt.isStructShockPlots = 0;
end
if ~isfield(opt,'isPlotTitle')
    opt.isPlotTitle = 0;
end
if ~isfield(opt,'isIRFcum')
    opt.isIRFcum = 0;
end

% consistency checks  

if opt.isFevdPlots == 2 || opt.isFevdPlots == 3 || opt.isFevdTexTable == 1
    opt.isFevdTable = 1; % calculations for FevdTable are required for these FevdPlots
end

% The algorithm differs depending whether or not zero restrictions are 
% used. So a boolean is created to allow to quickly vary the algorithm.
if any(any(any(Z)))
    hasZeroRestrictions = 1;
else
    hasZeroRestrictions = 0;
end

%% initialize random number generator
% The random number generator is initialized using the runNumber. This
% ensures that in different runs different random numbers are used.
rng(opt.runNumber); 

%% Estimation
%

% -----------------------------------------------------------------------
% Explanation:
% ===========
%   Each estimation method returns an array of 'nDrawsFromBvar' models 
%   (in the case of analytical solutions the array will only contain one 
%   model). Each model is a structure containing the following elements:
%   - model.A ... the reduced form parameter matrix.
%   - model.vc ... the variance-covariance matrix.
%   - model.nLags ... the number of lags in the model.
%   - model.res ... the residuals of the model.
%   - model.orthP ... an orthogonal matrix P with P*P' = vc.
% -----------------------------------------------------------------------

switch opt.estimationMethod
    case 'OLS'
        models = estimatevarols(y, opt.nLags, opt.hasConstant, opt.hasTrend);
    case 'diffuse'
        models = estimatevardiffuse(y, opt.nLags, opt.hasConstant, opt.nDrawsFromBvar, opt.hasTrend);
    otherwise
        error('Invalid estimation method specified.')
end

%}

if opt.runNumber == 1
	% If this is the first run, the impulse responses are generated for the
	% horizon needed to check the sign restrictions. To check the sign 
	% restrictions the impulse responses are needed in their stacked form
	% which is returned as the second object by calculateirf().
    for ii = 1:size(models,2)
        [  ~ , models(ii).fA , ~ ] = calculateirf(models(ii),opt.nMaxSignHorizon);
    end
	% variable that is to count the models found which fulfill imposed restrictions.
	nModelsFound = 0;
    nModelsChecked = 0;
	% If the option to run through all models drawn from the posterior was selected the
	% variable nDraw is initialized. This variable will increase by one each run until
	% all models have been ran through.
	if opt.drawFromPosterior == 0
		opt.nDraw = 1;
	end
	% Vector collecting all past runs. Needed to check whether a run is double.
    pastRuns = 1;
else
	% If this is not the first run, load previous results.
    load(strcat(opt.modelPath,opt.modelName,'/results'));
    if any(pastRuns == opt.runNumber)
        error('This run number was already used. Please choose another run number.')
    end
	% Add current run to the vector of past runs.
    pastRuns = [ pastRuns ; opt.runNumber ];
    clear opt.runNumber
end

%% Sign restrictions search
%
[ identifiedModels, nModelsFound, nModelsChecked ] = signrestrictionsearch(models, S, SN, Z, nVars, hasZeroRestrictions, nModelsFound, nModelsChecked, opt);
%}

%% Save results for later use or to continue with a further run
%
if opt.isSaveResults
    save(strcat(opt.modelPath, opt.modelName, '/results'));
end
%}

%% Estimation Output
%

identifiedModels = calculateStructShocks( identifiedModels, nVars);

if opt.isHistDecompTable || opt.isHistDecompPlots
    identifiedModels = calculateHistoricalDecomposition( identifiedModels, nVars, y);
end

[ output , ctmModel ] = generateOutput( y, identifiedModels, S, Z, nModelsFound, nVars, opt );

plotOutput( y, identifiedModels, output, nVars, t, nModelsFound, opt, ctmModel );

saveModelInfos( identifiedModels, S, Z, nModelsFound, nModelsChecked, nVars, opt );

%}