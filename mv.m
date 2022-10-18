
clear
clc

% Start with the default path
restoredefaultpath; 

% Path to the MATLAB asset pricing package
matlabPackagePath = 'D:\Published Repos\Moore-Velikov\'; 

% Path to the code for the paper
paperCodePath = 'D:\Published Repos\Moore-Velikov\Moore-Velikov'; 

% Add the relevant folders (with subfolders) to the path
addpath(genpath([matlabPackagePath, 'Data']))
addpath(genpath([matlabPackagePath, 'Functions']))
addpath(genpath([matlabPackagePath, 'Library Update']))
addpath(genpath([paperCodePath]))

% Navigate to the paper folder
cd(paperCodePath)

%% Add the directories

% Check if the /Data/ directory exists
if ~exist([pwd,'Data'], 'dir')
    mkdir(['Data'])
end

% Check if the /Results/ directory exists
if ~exist([pwd,'Results'], 'dir')
    mkdir(['Results'])
end

% Check if the /Figures/ directory exists
if ~exist([pwd,'Figures'], 'dir')
    mkdir(['Figures'])
end

% Make sure we add those to the path if we created them
addpath(genpath(pwd));

%% Start a log file

startLogFile([paperCodePath], 'mv_raps')
warning('off','all');

%% Make daily betas

run('make_daily_betas.m');

%% Make weighting function

run('make_weighting_function.m');

%% Make oil prices

run('make_oil_prices.m');

%% Make oil response forecast

run('make_oil_response_forecasts.m');

%% Make ad-hoc data/results

run('make_ad_hoc_data_results.m');

%% Make tables

run('make_tables.m');

%% Make figures

run('make_figures.m');

%% End the log file

diary('off');