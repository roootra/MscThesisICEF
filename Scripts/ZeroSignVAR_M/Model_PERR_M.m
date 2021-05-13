%% Startup

clear all; close all; clc; %dbstop if error; 
rng('default');
addpath 'functions';


%% General program settings
% Choose a unique model name.
opt.modelName = 'Model_PERR_M';     

% Specify the path to the folder, in which all result files should be saved.
opt.modelPath = strcat(pwd,'/');

%% VAR specification and other setup 

% Define labels for the variables in the same order as they enter the VAR.
opt.lVars = {'oil_USD_mom', 'imp_price_mom', 'neer_mom','intrate','gdp_per_cap_mom', 'cpi_all_mom'};

% Define labels for the identified shocks. The number of shock labels has 
% to be equal to the number of variables (i.e. label also residual shocks).
opt.lShocks = {'Global persistent shock', 'Import price shock', 'Ex. rate shock', 'Monetary shock', 'Aggregate supply shock', 'Aggregate demand shock'};

% % Define start and end date (MM-DD-YYYY).
opt.startDate = '03-01-2005';
opt.endDate = '08-01-2020';

% Define the lag order.
opt.nLags = 4;   

% Choose the estimation method ('OLS' or 'diffuse'), corresponding to 
% either ordinary least square or Bayesian estimations.
opt.estimationMethod = 'diffuse';

% Maximum horizon of sign restrictions (0 = only contemporary, 
% 1 = contemporary + horizon 1, etc.; Default = 0).
opt.nMaxSignHorizon = 1;

%% Data Setup
%

% Read in the data here. Assign the final data to be used in the analysis 
%   to y. 

%y = xlsread('data.xls','A1:D161');
y = xlsread('data_sign_and_zero_mon.xlsx', 'A1:F186')
[t,nVars] = size(y);


%% Sign Restriction Setup


S = zeros(2,nVars*(opt.nMaxSignHorizon + 1),nVars);
Z = zeros(2,nVars*(opt.nMaxSignHorizon + 1),nVars);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SIGN RESTRICTION DEFINITION
% -----------------------------------------------------------------------
% Explanation:
% ===========
% The sign and zero restrictions are implemented using selection matrices 
%   (S for sign restrictions and Z for zero restrictions).
%   S and Z are three dimensional arrays. 
%   - The first dimension is the row of the selection matrix. For each 
%   shock it should be incremented by one for each restriction imposed 
%   to identify the shock starting at 1 for the first restriction.
%   - The second dimension is the column of the selection matrix. It
%   selects the variable on which the restriction should be set as well as
%   the horizon. For example the second variable at horizon 0 is selected
%   with 2, the second variable at horizon 1 is selected with 2 + nVars, 
%   the second variable at horizon 2 is selected with 2 + 2*nVars, and the 
%   second variable at horizon h is selected with 2 + h*nVars, where nVars 
%   is the number of variables in the system.
%   - The third dimension is the number of the shock. For the first shock
%   you want to identify pick 1, for the second 2, and so on.
%   - The matrix SN is used to impose a restriction which demands the
%   impulse response to be above or below zero by a certain amount, e.g.,
%   to restrict the contemporaneous response of variable 2 to shock 6 to be
%   above 2 in addition to the sign restriction S(1,2,6) = 1 the
%   restriction SN(1,2,6) = 2 would be set.
% Examples:
% ========
%   First shock: zero restrictions on variables 1 and 2 at horizon 0
%       Z(1,1,1) = 1;
%       Z(2,2,1) = 1;
%   Second shock: negative reaction of variable 1 up to maximum horizon
%     for ii = 1:(nMaxSignHorizon+1) 
%         S(ii,1+nVars*(ii-1),2) = -1;
%     end
%   second shock: positive reaction of variable 4 at horizons 0 and 1
%       S(nMaxSignHorizon+2,4,2);
%       S(nMaxSignHorizon+3,4+nVars,2);
%       Note: As we have already nMaxSignHorizon+1 restrictions imposed
%       for the second shock we have to start at nMaxSignHorizon+2.
%   For shocks without restrictions (residual shocks) nothing has to be
%   done, but they may be normalized.
% -----------------------------------------------------------------------


%%%Russia PERR short-run%%%

%Oil shock
S(1:4,1,1) = 1 %oil (normalized)
%S(2,2,1) = 1 %imp.price
S(3:4,3,1) = -1 %NEER
%S(4,4,1) = 1 %int.rate
%Z(1,5,1) = 1 %gdp
%Z(2,6,1) = 1 %cpi

%Import price shock
Z(1,1,2) = 1 %oil 
S(1,2,2) = 1 %imp.price (normalized)
%S(3,3,2) = -1 %NEER
%S(4,4,2) = 1 %int.rate
%Z(1,5,2) = 1 %gdp
S(2,6,2) = 1 %cpi

%Exchange rate shock
Z(1,1,3) = 0 %oil
%S(2,2,1) = 1 %imp.price
S(1,3,3) = 1 %NEER (normalized)
%S(4,4,1) = 1 %int.rate
%Z(1,5,1) = 1 %gdp
S(2,6,3) = 1 %cpi

%Interest rate shock
Z(1,1,4) = 0 %oil 
%S(2,2,1) = 1 %imp.price
%S(3,3,1) = -1 %NEER
S(1,4,4) = 1 %int.rate (normalized)
S(2,5,4) = -1 %gdp
S(3,6,4) = -1 %cpi

%Output shock
%Z(1,1,5) = 1 %oil 
%S(2,2,1) = 1 %imp.price
%S(3,3,1) = -1 %NEER
%S(4,4,1) = 1 %int.rate
S(1,5,5) = 1 %gdp (normalized)
S(2,6,5) = -1 %cpi ???

%Inflation exogenous shock
%Z(1,1,6) = 1 %oil (normalized)
%S(2,2,1) = 1 %imp.price
S(1,3,6) = 1 %NEER
S(2,4,6) = 1 %int.rate
S(3,5,6) = 1 %gdp
S(4,6,6) = 1 %cpi (normalized)

opt.isImpulsePlots = 0;
opt.isHistDecompPlots = 0;
opt.nModelDraws = 100;
opt.nDrawsFromBvar = 1e3;
opt.nTransformationsPerDraw = 1e3;
%opt.isNoTransform = 1;
signrestrictions

%%% IRF Calculation
irf_oil = sum(output.irfCTM(1:12,13))/sum(output.irfCTM(1:12,31));
irf_imp = sum(output.irfCTM(1:12,14))/sum(output.irfCTM(1:12,32));
irf_neer = sum(output.irfCTM(1:12,15))/sum(output.irfCTM(1:12,33));
irf_int = sum(output.irfCTM(1:12,16))/sum(output.irfCTM(1:12,34));
irf_gdp = sum(output.irfCTM(1:12,17))/sum(output.irfCTM(1:12,35));
irf_cpi = sum(output.irfCTM(1:12,18))/sum(output.irfCTM(1:12,36));

sprintf("PERRs:\tOil\tImp. price\tNEER\tInt. rate\tGDP\tCPI\n\t%0.5f\t%0.5f\t%0.5f\t%0.5f\t%0.5f\t%0.5f", irf_oil, irf_imp, irf_neer, irf_int, irf_gdp, irf_cpi)
%sprintf("\t%0.5f\t%0.5f\t%0.5f\t%0.5f\t%0.5f\t%0.5f", irf_oil, irf_imp, irf_neer, irf_int, irf_gdp, irf_cpi)