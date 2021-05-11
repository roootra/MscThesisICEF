%% Startup

clear all; close all; clc; %dbstop if error; 
rng('default');
addpath 'functions';


%% General program settings
% Choose a unique model name.
opt.modelName = 'Model_PERR_Q';     

% Specify the path to the folder, in which all result files should be saved.
opt.modelPath = strcat(pwd,'/');

%% VAR specification and other setup 

% Define labels for the variables in the same order as they enter the VAR.
opt.lVars = {'oil_USD_qoq', 'imp_price_qoq', 'neer_qoq','miacr_31'm'reserves_USD_qoq','broad_money_SA','gdp_nominal_qoq', 'cpi_all_qoq'};

% Define labels for the identified shocks. The number of shock labels has 
% to be equal to the number of variables (i.e. label also residual shocks).
opt.lShocks = {'Oil shock', 'Import price shock', 'Ex. rate shock', 'Monetary shock', 'Reserves shock', 'Money supply shock', 'Supply shock', 'Demand shock'};

% % Define start and end date (MM-DD-YYYY).
opt.startDate = '30-06-2005';
opt.endDate = '31-12-2020';

% Define the lag order.
opt.nLags = 1;   

opt.nDrawsFromBvar = 1e3;
opt.nTransformationsPerDraw = 1e3;

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
y = csvread('data_sign_and_zero.csv')
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
S(1,1,1) = 1 %oil self
S(2,2,1) = 1 %oil is costly => higher import prices %???
S(3,3,1) = -1 %oil is purchased => demand for rouble => lower NEER
%S(4,4,1) = 1 %oil price is higher => budget rule works => higher reserves
%Z(1,5,1) = 1
%Z(2,6,1) = 1
%S(1,7,1) = 
%S(1,8,1) = 

%Imp. price shock
Z(1,1,2) = 1 %import prices have no effect to oil
S(1,2,2) = 1 %imp. prices self
%S(2,3,2) = 1 %NEER?
Z(3,4,2) = 1 %reserves
Z(4,5,2) = 1 %money
%S(2,6,2) = 1 %int. rate?
%S(2,7,2) = 
S(2,8,2) = 1 %heats up cpi

%NEER shock
Z(1,1,3) = 1 %oil
Z(2,2,3) = 1 %imp. prices
S(1,3,3) = 1 %normalization
%S(3,4,3) = 1 %reserves? may be -
Z(3,5,3) = 1 %money
%S(3,6,3) = 1 %int rate should be higher to return rub investors
%S(3,7,3) = %gdp? + in LR???
S(2,8,3) = 1 %heats up prices

%Int. rate shock
Z(1,1,4) = 1
Z(2,2,4) = 1
S(1,3,4) = -1 %NEER
S(2,4,4) = 1 %normalization
%Z(3,5,4) = 1 %reserves
Z(4,6,4) = 1 %int rate
Z(5,7,4) = 1 %gdp
Z(5,8,4) = 1 %cpi

%Reserves shock
Z(1,1,5) = 1
Z(2,2,5) = 1
S(1,3,5) = 1 %more money => more rouble supply?
%Z(3,4,5) = 1 %reserves
S(2,5,5) = 1 %normalization
%S(3,6,5) = 1 %Int. rate
%S(4,7,5) = %in LR should be 0!!!
S(5,8,5) = 1 %money supply heats up prices

%Broad money shock
Z(1,1,6) = 1
Z(2,2,6) = 1
%S(1,3,6) = -1 %more attractive for foreign investors
Z(3,4,6) = 1 %reserves
S(2,5,6) = -1 %money?
S(3,6,6) = 1 %normalization
Z(3,7,6) = 1 %no contemporaneous reaction of GDP
S(4,8,6) = -1 %cools down prices, money are costly => less credits

%GDP shock
Z(1,1,7) = 1
Z(2,2,7) = 1
%S(1,3,7) = -1 %neer
%Z(3,4,7) = 1 %reserves
%S(2,5,7) = 1 %more output => more money???
%S(3,6,7) = 1 %int. rate
S(4,7,7) = 1 %normalization
%S(5,8,7) = 1 %cpi

%CPI shock
Z(1,1,8) = 1 
Z(2,2,8) = 1
%S(1,3,8) = -1 %NEER, may be 0
Z(3,4,8) = 0 %reserves
Z(4,5,8) = 0 %money
%S(2,6,8) = 1 %cb resists infl.
S(3,7,8) = 1 %Forbes 2020
S(4,8,8) = 1 %normalization
%{
%}

opt.isImpulsePlots = 0;
opt.nModelDraws = 500;
opt.isNoTransform = 1;
signrestrictions


