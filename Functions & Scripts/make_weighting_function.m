%% Timekeeping

fprintf('Now working on the weighting functions. Run started @ %s.\n\n\n',char(datetime('now')));

%%
clear
clc

load FF10
load FF17
load FF49
load SIC
load ddates
load dret
load dff
load dcapmb 
load permno

% Make the 3-day CAR 
xret = dret - dcapmb .* repmat(dmkt,1,size(dret,2));
xret = (1+lag(xret,1,nan)).*(1+xret).*(1+lead(xret,1,nan))-1;
clear dret dcapmb 

% Store a few constants
nDays   = length(ddates);
nStocks = length(permno);
nObs    = nDays * nStocks;

% Crate a table with the CAR3s for all days 
xret_tab = array2table(reshape(xret, nObs, 1), 'VariableNames', {'xret'});
xret_tab.permno = reshape(repmat(permno', nDays, 1), nObs, 1);
xret_tab.ddates = reshape(repmat(ddates, 1, nStocks), nObs, 1);
isNanXret = isnan(xret_tab.xret);
xret_tab(isNanXret, :) = [];

% Initialize the structure
wghtFctns = struct;

% Loop over the six industry classifications
for i = 1:6
    % Choose the industry classification
    switch i
        case 1
            inds = FF10;
            wghtFctns(i).label = {'FF10'};
        case 2
            inds = FF17;
            wghtFctns(i).label = {'FF17'};            
        case 3 
            inds = FF49;
            wghtFctns(i).label = {'FF49'};
        case 4
            inds = floor(SIC/100);
            wghtFctns(i).label = {'SIC2d'};            
        case 5
            inds = floor(SIC/10);
            wghtFctns(i).label = {'SIC3d'};            
        case 6 
            inds = SIC;            
            wghtFctns(i).label = {'SIC4d'};            
    end

    % Timekeeping
    fprintf('Now working on %s. Run started at %s.\n\n',char(wghtFctns(i).label),char(datetime('now')));    
    
    % Get the wCAR3 matrix % store it in the structure
    wCAR3 = makeWghtFctn(inds, xret_tab);
    wghtFctns(i).wCAR3 = wCAR3;
    
    % Timekeeping
    fprintf('Done with %s at %s.\n\n\n',char(wghtFctns(i).label),char(datetime('now')));
    save Data\wghtFctns wghtFctns
    
end

save Data\wghtFctns wghtFctns

