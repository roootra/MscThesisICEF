%% Startup

clear all; close all; clc; %dbstop if error; 
rng('default');
addpath 'functions';


%% General program settings
% Choose a unique model name.
opt.modelName = 'Model_PERR_Q_reduced';     

% Specify the path to the folder, in which all result files should be saved.
opt.modelPath = strcat(pwd,'/');

%% VAR specification and other setup 

% Define labels for the variables in the same order as they enter the VAR.
opt.lVars = {'oil_USD_qoq', 'imp_price_qoq', 'neer_qoq', 'miacr_31','gdp_nominal_qoq', 'cpi_all_qoq'};

% Define labels for the identified shocks. The number of shock labels has 
% to be equal to the number of variables (i.e. label also residual shocks).
opt.lShocks = {'Oil shock', 'Import price shock', 'Ex. rate shock', 'Monetary shock', 'Supply shock', 'Demand shock'};

% % Define start and end date (MM-DD-YYYY).
opt.startDate = '30-06-2005';
opt.endDate = '31-12-2020';

% Define the lag order.
opt.nLags = 1;   

opt.nDrawsFromBvar = 1e5;
opt.nTransformationsPerDraw = 1e3;

% Choose the estimation method ('OLS' or 'diffuse'), corresponding to 
% either ordinary least square or Bayesian estimations.
opt.estimationMethod = 'OLS';

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
S(2,2,1) = 1 %import prices
S(3,3,1) = -1 %NEER
%S(4,4,1) = %interest rate
%S(5,5,1) = %GDP
%S(6,6,1) = %CPI

%Imp. price shock
Z(1,1,2) = 0 %no effect to oil
S(2,2,2) = 1 %normalization
%S(2,3,1) = 1 %NEER
%S(4,4,1) = %interest rate
%S(5,5,1) = %GDP
S(6,6,2) = 1 %CPI

%NEER shock
Z(1,1,3) = 1 %no effect to oil
Z(2,2,3) = 1 %no effect to imp. prices
S(3,3,3) = 1 %normalization
%S(4,4,1) =  %int. rate
%S(5,5,1) = %GDP
S(6,6,3) = 1 %CPI

%Int. rate shock
Z(1,1,4) = 1 %no effect to oil
Z(2,2,4) = 1 %no effect to imp. prices
%S(1,3,4) = -1 %NEER
S(4,4,4) = 1 %normalization
S(5,5,4) = -1 %GDP ok
S(6,6,4) = -1 %CPI

%GDP shock
Z(1,1,5) = 1 %no effect to oil
Z(2,2,5) =1 %no effect to imp. prices
%S(1,3,4) = -1 %NEER
%S(2,4,4) = 1 %int. rate
S(5,5,5) = -1 %normalization
%S(4,6,4) = -1 %CPI

%CPI shock
Z(1,1,6) = 1 %no effect to oil
Z(2,2,6) = 1 %no effect to imp. prices
%S(1,3,6) = 1 %NEER
%S(2,4,6) = 1 %int. rate
%S(3,5,5) = -1 %GDP
S(6,6,6) = 1 %normalization

signrestrictions


