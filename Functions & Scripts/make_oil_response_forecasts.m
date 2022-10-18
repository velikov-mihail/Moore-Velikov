%% Timekeeping

fprintf('Now working on the betas and oil response forecasts. Run started @ %s.\n\n\n',char(datetime('now')));


%% 

clear
clc

load wghtFctns
load wtiMonthlyPrices

horizon      = 12 * 3; % 12 quarters
min_quarters = 12; % 12 quarters minimum

varStruct = wghtFctns;
for i = 1:length(varStruct)
    varStruct(i).wtiMonthlyBetas = makeCommodityBetas(wghtFctns(i).wCAR3, wtiMonthlyPrices, horizon, min_quarters);
    varStruct(i).predictedCAR = makePredictedCAR(varStruct(i).wtiMonthlyBetas); 
end

save Results\varStruct varStruct -v7.3
