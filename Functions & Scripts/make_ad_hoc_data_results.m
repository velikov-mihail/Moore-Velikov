%% Timekeeping

fprintf('Now working on the ad-hoc results. Run started @ %s.\n\n\n',char(datetime('now')));

%% Make daily strategy portfolio returns for Figure 6

clear
clc

load ret
load me
load dates
load nyse
load ddates
load dme
load dret
load varStruct

% Determine the start date & number of portfolios
s = find(dates==197501);
nPtf = 5;


% Find the FF49 results
r = find(strcmp([varStruct.label],{'FF49'}));
oilChanges = varStruct(r).wtiMonthlyBetas.quarterlyChanges;
oilResponseForecast = varStruct(r).predictedCAR.Changes;


% Run the time-series regressions
ind = makeUnivSortInd(oilResponseForecast, nPtf, NYSE);
res = runUnivSort(ret, ind, dates, me, 'plotFigure', 0, ...
                                       'printResults', 0, ...
                                       'factorModel', 1);

s=find(dates==197401);

ptf_drets=nan(length(ddates), nPtf);

for i=find(floor(ddates/100)==dates(s),1):length(ddates)-1
    currMonth=find(dates==floor(ddates(i)/100));
    lastQuarter=currMonth-find(ismember(flipud(dates(currMonth-3:currMonth-1)-100*floor(dates(currMonth-3:currMonth-1)/100)),[3 6 9 12]));
    for j=1:nPtf        
        index=(ind(lastQuarter,:)==j) & isfinite(dret(i,:)) & isfinite(me(lastQuarter,:));
        ptf_drets(i,j)=sum(me(lastQuarter,index).*dret(i,index)/sum(me(lastQuarter,index)));                
    end
end

ptf_drets(:,nPtf+1)=ptf_drets(:,1)-ptf_drets(:,nPtf);


monthlyFromDailyReturns=d2wret(ptf_drets,eomflag);
res = runUnivSort(ret, ind, dates, me, 'plotFigure', 0, ...
                                       'printResults', 0, ...
                                       'factorModel', 1);

% for i=1:5
%     prt(nanols(res.pret(:,i),[ones(size(res.pret(:,1))) monthlyFromDailyReturns(:,i)]));
% end

save Data\ptf_drets ptf_drets

%% Make weighted turnover for Table 10

clear
clc

load FF49
load ddates
load dff
load permno
load dvol
load dshrout

dvol(dvol<=0) = nan;
dshrout(dshrout<=0) = nan;

dshrout = dshrout*1000;
xret = dvol./dshrout;
clear dvol dshrout 

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


% Get the wCAR3 matrix % store it in the structure
weightedTurnover = makeWghtFctn(FF49, xret_tab);
save Data\weightedTurnover weightedTurnover

%% Make recession subsample performance

clear
clc

load ret
load me
load dates
load nyse
load ff
load varStruct

% Determine the start date & number of portfolios
s = find(dates==197501);
nPtf = 5;


% Find the FF49 results
r = find(strcmp([varStruct.label],{'FF49'}));
oilChanges = varStruct(r).wtiMonthlyBetas.quarterlyChanges;
oilResponseForecast = varStruct(r).predictedCAR.Changes;


% Run the time-series regressions
ind = makeUnivSortInd(oilResponseForecast, nPtf, NYSE);
res = runUnivSort(ret, ind, dates, me, 'plotFigure', 0, ...
                                       'printResults', 0, ...
                                       'factorModel', 1);

% Account for the zero-oil-price-change quarters
pret = res.pret;
indZeroOilQtrChange = find(oilChanges == 0 & ...
                          dates>197412);
pret(indZeroOilQtrChange+1,end-1:end)=repmat(rf(indZeroOilQtrChange+1),1,2);
pret(indZeroOilQtrChange+2,end-1:end)=repmat(rf(indZeroOilQtrChange+2),1,2);
pret(indZeroOilQtrChange+3,end-1:end)=repmat(rf(indZeroOilQtrChange+3),1,2);
pret(indZeroOilQtrChange+1,1:end-2)=0;
pret(indZeroOilQtrChange+2,1:end-2)=0;
pret(indZeroOilQtrChange+3,1:end-2)=0;


% Get the GDP data
fredStruct = getFredData('USREC','1970-01-02','2021-12-31', 'lin', [],'eop');
fredData = array2table(fredStruct.Data, 'VariableNames',[{'date'},{'recInd'}]);
fredData.date = datetime(fredData.date, 'ConvertFrom', 'datenum');
fredData.date = fredData.date + calmonths(2);
data = fredData;
data.dates=100*year(data.date)+month(data.date);
[~,ia,ib]=intersect(dates,data.dates);
recInd = nan(size(dates));
recInd(ia) = data.recInd(ib);

strat_ret = pret(:,end);
strat_ret(dates<197501) = nan;

ind = (recInd == 1);
prt(nanols(strat_ret(ind),const(ind)))
prt(nanols(strat_ret(~ind),const(~ind)))

%%

clear
clc

load ret
load me
load dates
load nyse
load ff
load varStruct
load FF49

% Determine the start date & number of portfolios
s = find(dates==197501);
nPtf = 5;

% Find the FF49 results
r = find(strcmp([varStruct.label],{'FF49'}));
oilChanges = varStruct(r).wtiMonthlyBetas.Changes;
oilResponseForecast = varStruct(r).predictedCAR.Changes;

indRet = nan(49, 3);
indTRet = nan(49, 3);
indMean = nan(49, 1);
indStd = nan(49, 1);
indIQR = nan(49, 1);
indKurt = nan(49, 1);

for i = 1:49
    temp = oilChanges;
    temp(FF49~=i) = nan;
    indMean(i) = nanmean(nanmedian(temp(s:end,:),2));
    indStd(i) = nanmean(nanstd(temp(s:end,:),[],2));
    indIQR(i) = nanmean(iqr(temp(s:end,:),2));
    indKurt(i) = nanmean(kurtosis(temp(s:end,:),[],2));

    temp = oilResponseForecast;
    temp(FF49~=i) = nan;
    
    % Run the time-series regressions
    ind = makeUnivSortInd(temp, nPtf, NYSE);
    ind = 1*(ind==1) + 2*(ind==nPtf);
    res = runUnivSort(ret, ind, dates, me, 'plotFigure', 0, ...
                                           'printResults', 0, ...
                                           'factorModel', 1);

    % Account for the zero-oil-price-change quarters
    pret = res.pret;
    indZeroOilQtrChange = find(oilChanges == 0 & ...
                              dates>197412);
    pret(indZeroOilQtrChange+1,end-1:end)=repmat(rf(indZeroOilQtrChange+1),1,2);
    pret(indZeroOilQtrChange+2,end-1:end)=repmat(rf(indZeroOilQtrChange+2),1,2);
    pret(indZeroOilQtrChange+3,end-1:end)=repmat(rf(indZeroOilQtrChange+3),1,2);
    pret(indZeroOilQtrChange+1,1:end-2)=0;
    pret(indZeroOilQtrChange+2,1:end-2)=0;
    pret(indZeroOilQtrChange+3,1:end-2)=0;
    
    pret(dates<197501,:)=nan;
    for j = 1:size(pret,2)
        tempRes = nanols(pret(:,j), const);
        res(i).xret(j) = tempRes.beta;
        res(i).txret(j) = tempRes.tstat;
    end
    
    indRet(i,:) = res(i).xret;
    indTRet(i,:) = res(i).txret;    
end

x = indKurt;
y = indTRet(:,3);
z = indMean;

h1 = scatter(x, y, 'filled');
h2 = lsline;
res = nanols(y,[ones(size(y)) x]);
% prt(res)

bubblechart(x, y, z);

barh(indRet(:,3));
set(gca,'ytick',1:49);
set(gca,'yticklabel',FF49Names.shortName)

barh(indKurt);
set(gca,'ytick',1:49);
set(gca,'yticklabel',FF49Names.shortName)

barh(indStd);
set(gca,'ytick',1:49);
set(gca,'yticklabel',FF49Names.shortName)

barh(indTRet(:,3));
set(gca,'ytick',1:49);
set(gca,'yticklabel',FF49Names.shortName)
