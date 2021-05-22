%% Startup

clear all; close all; clc; %dbstop if error; 
rng('default');
addpath 'functions';


%% General program settings
% Choose a unique model name.
opt.modelName = 'Model_PERR_Q_rc_rubusd';     

% Specify the path to the folder, in which all result files should be saved.
opt.modelPath = strcat(pwd,'/');

%% VAR specification and other setup 

% Define labels for the variables in the same order as they enter the VAR.
%opt.lVars = {'oil_USD_qoq', 'miacr_31', 'neer_qoq','gdp_real_qoq', 'cpi_all_qoq'};
opt.lVars = {'gdp_real_qoq', 'cpi_all_qoq', 'miacr_31','rubusd_qoq', 'oil_USD_qoq'}; %a-la Forbes et al. 2018

% Define labels for the identified shocks. The number of shock labels has 
% to be equal to the number of variables (i.e. label also residual shocks).
%opt.lShocks = {'Oil shock', 'Monetary shock', 'Ex. rate shock', 'Supply shock', 'Demand shock'};
opt.lShocks = {'Supply shock', 'Core demand shock', 'Monetary shock', 'Exchange rate shock', 'Oil price shock'};

% % Define start and end date (MM-DD-YYYY).
opt.startDate = '30-06-2005';
opt.endDate = '31-12-2020';

% Define the lag order.
opt.nLags = 1;   

opt.nDrawsFromBvar = 500;
opt.nTransformationsPerDraw = 1000;

% Choose the estimation method ('OLS' or 'diffuse'), corresponding to 
% either ordinary least square or Bayesian estimations.
opt.estimationMethod = 'diffuse';

% Maximum horizon of sign restrictions (0 = only contemporary, 
% 1 = contemporary + horizon 1, etc.; Default = 0).
opt.nMaxSignHorizon = 8;

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

%opt.lShocks = {'Oil shock', 'Monetary shock', 'Ex. rate shock', 'Supply shock', 'Demand shock'};

%%%Russia PERR short-run%%%

%opt.lShocks = {'Supply shock', 'Demand shock', 'Monetary shock', 'Exchange rate shock', 'Oil price shock'};

%Supply shock
S(1,1,1) = 1
S(4,6,1) = 1
Z(1,5 + nVars*0,1) = 1
%Z(2,5 + nVars*1,1) = 1
%Z(3,5 + nVars*2,1) = 1
%Z(4,5 + nVars*3,1) = 1
%NEW
S(2,2,1) = -1 %inflation response
S(3,7,1) = -1

%Demand shock
S(1,2,2) = 1
S(6,7,2) = 1
S(2,3,2) = 1 %interest rate response
S(3,8,2) = 1
Z(1,5 + nVars*0,2) = 1
%Z(2,5 + nVars*1,2) = 1
%Z(3,5 + nVars*2,2) = 1
%Z(4,5 + nVars*3,2) = 1
Z(2,1 + nVars*8,2) = 1
%NEW
S(4,1,2) = 1 %output response
S(5,6,2) = 1
S(6,4,2) = -1 %NEER response (this)
S(7,9,2) = -1

%Monetary shock
S(1,3,3) = 1
S(8,8,3) = 1
S(2,1,3) = -1 %output response
S(3,2,3) = -1 %inflation response
S(4,6,3) = -1 
S(5,7,3) = -1
Z(1,5 + nVars*0,3) = 1
%Z(2,5 + nVars*1,3) = 1
%Z(3,5 + nVars*2,3) = 1
%Z(4,5 + nVars*3,3) = 1
Z(2,1 + nVars*8,3) = 1
%NEW
S(6,4,3) = -1 %NEER response
S(7,9,3) = -1

%Ex. rate shock
S(1,4,4) = 1
S(4,9,4) = 1
S(2,2,4) = 1
S(3,7,4) = 1
Z(1,5 + nVars*0,4) = 1
%Z(2,5 + nVars*1,4) = 1
%Z(3,5 + nVars*2,4) = 1
%Z(4,5 + nVars*3,4) = 1
Z(2,1 + nVars*8,4) = 1
%NEW
S(4,3,4) = 1 %interest rate response
S(5,8,4) = 1

%Oil price shock
S(1,5,5) = 1
S(3,10,5) = 1
S(2,4,5) = -1



opt.isImpulsePlots = 1;
opt.isIRFtable = 1;
opt.isStructShockTable = 1;
opt.isSaveResults = 1;
opt.isFevdTable = 1;
opt.isFevdTexTable = 1;
opt.isHistDecompTable = 1;
opt.isPlotTitle = 1;
opt.isFevdPlots = 1;
opt.isStructShockPlots = 1;
opt.isHistDecompPlots = 1;

signrestrictions


%%% PERR Calculation
perr_gdp = sum(output.irfCTM(1:4,6))/sum(output.irfCTM(1:4,16));
perr_cpi = sum(output.irfCTM(1:4,7))/sum(output.irfCTM(1:4,17));
perr_int = sum(output.irfCTM(1:4,8))/sum(output.irfCTM(1:4,18));
perr_neer = sum(output.irfCTM(1:4,9))/sum(output.irfCTM(1:4,19));
perr_oil = sum(output.irfCTM(1:4,10))/sum(output.irfCTM(1:4,20));

sprintf("PERRs (4q):\n")
sprintf("Oil shock: %0.5f", perr_oil)
sprintf("Monetary shock: %0.5f", perr_int)
sprintf("Exchange rate shock: %0.5f", perr_neer)
sprintf("Output shock: %0.5f", perr_gdp)
sprintf("Inflation shock: %0.5f", perr_cpi)

perr_gdp_s = sum(output.irfCTM(1,6))/sum(output.irfCTM(1,16));
perr_cpi_s = sum(output.irfCTM(1,7))/sum(output.irfCTM(1,17));
perr_int_s = sum(output.irfCTM(1,8))/sum(output.irfCTM(1,18));
perr_neer_s = sum(output.irfCTM(1,9))/sum(output.irfCTM(1,19));
perr_oil_s = sum(output.irfCTM(1,10))/sum(output.irfCTM(1,20));

sprintf("PERRs (1q):\n")
sprintf("Oil shock: %0.5f", perr_oil_s)
sprintf("Monetary shock: %0.5f", perr_int_s)
sprintf("Exchange rate shock: %0.5f", perr_neer_s)
sprintf("Output shock: %0.5f", perr_gdp_s)
sprintf("Inflation shock: %0.5f", perr_cpi_s)