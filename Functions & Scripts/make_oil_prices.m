%% Timekeeping

fprintf('Now working on the oil prices. Run started @ %s.\n\n\n',char(datetime('now')));

%%
clear
clc

load dates
load ddates

% Series WTISPLC on FRED is constructed from OILPRICE and MCOILWTICO. See 
% https://fred.stlouisfed.org/series/WTISPLC:
% "This series was created by the Federal Reserve Bank of St. Louis to 
% expand the history of the monthly West Texas Intermediate oil price 
% series in FRED. We simply combined these two FRED series: 
% https://fred.stlouisfed.org/series/OILPRICE and 
% https://fred.stlouisfed.org/series/MCOILWTICO. From January 1946 through 
% July 2013, the series used is OILPRICE. From August 2013 to present, the 
% series used is MCOILWTICO."

% Get OILPRICE until 7/2013
fredStruct=getFredData('OILPRICE', '1970-01-01', '2013-07-01', 'lin', 'm', 'eop');
tempData = array2table(fredStruct.Data, 'VariableNames', [{'date'},'WTISPLC']);
tempData.date = datetime(tempData.date, 'ConvertFrom', 'datenum');

% Get MCOILWTICO starting in 8/2013:
fredStruct=getFredData('MCOILWTICO', '2013-08-01', [], 'lin', 'm', 'eop');
tempData2 = array2table(fredStruct.Data, 'VariableNames', [{'date'},'WTISPLC']);
tempData2.date = datetime(tempData2.date, 'ConvertFrom', 'datenum');

% Combine the two:
data = [tempData; tempData2];
        
% Initialize the structure 
wtiMonthlyPrices = struct;
% wtiMonthlyPrices.daily=nan(size(ddates));
wtiMonthlyPrices.monthly   = nan(size(dates));
wtiMonthlyPrices.quarterly = nan(size(dates));

yyyymm = 100 * year(data.date) + ...
               month(data.date);
           
% Intersect the dates           
[~, ia, ib] = intersect(dates, yyyymm);           

% s=find(dates==rawPrices(1,1));
wtiMonthlyPrices.monthly(ia) = data.WTISPLC(ib);

% Get the quarterly index
wtiMonthlyPrices.quarterly = wtiMonthlyPrices.monthly;
monthIndex = dates - 100*floor(dates/100);
isQuarterEnd = ismember(monthIndex, [3 6 9 12]);
wtiMonthlyPrices.quarterly(~isQuarterEnd) = nan;

% Get the changes and percent changes:
wtiMonthlyPrices.quarterlyChanges = wtiMonthlyPrices.quarterly - lag(wtiMonthlyPrices.quarterly,3,nan);
wtiMonthlyPrices.quarterlyReturns = wtiMonthlyPrices.quarterly ./ lag(wtiMonthlyPrices.quarterly,3,nan) - 1;

save Data/wtiMonthlyPrices wtiMonthlyPrices
