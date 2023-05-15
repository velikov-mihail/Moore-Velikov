%% Start a diary (will have all the table results)
clear
clc

diary Results/MooreVelikovTablesOutput.txt

fprintf('\n\n\n\nTable printing started @ %s\n\n\n',char(datetime('now')));

%% Print Table 1 - portfolio characteristics

fprintf('\n\n\nTable 1 output:\n\n\n');

clear
clc

load ret
load me
load dates
load nyse
load varStruct

r = find(strcmp([varStruct.label],{'FF49'}));

wCAR3 =varStruct(r).wCAR3;
wtiBetas = varStruct(r).wtiMonthlyBetas.Changes;
oilResponseForecast = varStruct(r).predictedCAR.Changes;

ind = makeUnivSortInd(oilResponseForecast, 5, NYSE);
res = runUnivSort(ret, ind, dates, me, 'timePeriod', [197412 202112], ...
                                       'plotFigure', 0, ...
                                       'printResults', 0);

indFinite = find(sum(ind,2)>0);
index = zeros(size(ret));
for i = 1:length(indFinite)
    index(indFinite(i)+3,:)=ind(indFinite(i),:);        
end

a=[mean(res.ptfNumStocks(:,1:5), 1, 'omitnan'); mean(res.ptfMarketCap(:, 1:5), 1, 'omitnan')/1e6; 100*prt_char(wCAR3,ind); 100*prt_char(wtiBetas,ind); 100*prt_char(oilResponseForecast,ind);];
h={'\# of Firms','Market Capitalization','Industry-Weighted CAR3 ($wCAR3_{i,Q}$)','Oil Price Exposure ($\hat\beta_{i,Q}$)','Oil-Response Forecast ($\hat u_{i,Q}^\text{Oil}$)'};
mat2Tex(a, a, h, 2);

%% Table 2 - time-series regressions

fprintf('\n\n\nTable 2 output:\n\n\n');

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

printPtfResults(pret, s, 1)

%% Table 3 - Fama-MacBeth regressions

fprintf('\n\n\nTable 3 output:\n\n\n');

clear
clc

load ret
load me
load dates
load nyse
load ff
load varStruct
load bm
load R
load dates
load GP
load AT
load FinFirms

% Gross profitability
gp = GP./AT; 
gp(FinFirms == 1) = nan; 

% Asset growth
assetGrowth = -AT./lag(AT,12,nan);
assetGrowth(FinFirms==1) = nan;
inv = - assetGrowth;

% Find the FF49 results
r = find(strcmp([varStruct.label],{'FF49'}));
oilChanges = varStruct(r).wtiMonthlyBetas.quarterlyChanges;
oilResponseForecast = varStruct(r).predictedCAR.Changes;

% Create the rank variable
orp_rank = tiedrank(oilResponseForecast')'/1000;

% Account for the zero oil-price-change quarters
indZeroOilQtrChange = find(oilChanges == 0 & ...
                          dates>197412);
orp_rank(indZeroOilQtrChange,   :) = nan;
orp_rank(indZeroOilQtrChange+1, :) = nan;
orp_rank(indZeroOilQtrChange+2, :) = nan;

yret = ret;
yret(indZeroOilQtrChange+1, :) = nan;
yret(indZeroOilQtrChange+2, :) = nan;
yret(indZeroOilQtrChange+3, :) = nan;

% Run the regressions
res1 = runFamaMacBeth(100*yret, [orp_rank],                              dates, 'timePeriod', 197501, ...
                                                                                'weightMatrix', me, ...
                                                                                'neweyWestLags', 3, ...
                                                                                'printResults', 0);
res2 = runFamaMacBeth(100*yret, [log(me) log(bm) gp inv R ret],          dates, 'timePeriod', 197501, ...
                                                                                'weightMatrix', me, ...
                                                                                'neweyWestLags', 3, ...
                                                                                'printResults', 0);
res3 = runFamaMacBeth(100*yret, [orp_rank log(me) log(bm) gp inv R ret], dates, 'timePeriod', 197501, ...
                                                                                'weightMatrix', me, ...
                                                                                'neweyWestLags', 3, ...
                                                                                'printResults', 0);


a=[res1.bhat(2) nan           res3.bhat(2);
   nan          res2.bhat(2)  res3.bhat(3); 
   nan          res2.bhat(3)  res3.bhat(4);
   nan          res3.bhat(4)  res3.bhat(5);
   nan          res2.bhat(5)  res3.bhat(6);
   nan          res2.bhat(6)  res3.bhat(7);
   nan          res2.bhat(7)  res3.bhat(8);
   sum(isfinite(res1.beta(:,1))) sum(isfinite(res2.beta(:,1))) ...
   sum(isfinite(res3.beta(:,1)))];
   
tA=[res1.t(2)   nan        res3.t(2);
   nan          res2.t(2)  res3.t(3); 
   nan          res2.t(3)  res3.t(4);
   nan          res2.t(4)  res3.t(5);
   nan          res2.t(5)  res3.t(6);
   nan          res2.t(6)  res3.t(7);
   nan          res2.t(7)  res3.t(8);
   nan(1,3)];

h=[{'Rank ($\hat{u}_{i,Q}^{\DeltaP_{Oil}}$)'}, ...
   {'log(ME)'}, ...
   {'log(B/M)'}, ...
   {'GP/A'}, ...
   {'Investment'}, ...
   {'$r_{12,1}$'}, ...
   {'$r_{1,0}$'}, ...
   {'$n$'}];
a = [nan(8,1) a];
tA = [nan(8,1) tA];

mat2Tex(a,tA,h,2);

%% Table 4 - Controlling for other anomalies
 
fprintf('\n\n\nTable 4 output:\n\n\n');

clear
clc

load ret
load dates
load nyse
load me
load ff
load varStruct

% Find the FF49 results
r = find(strcmp([varStruct.label],{'FF49'}));
oilChanges = varStruct(r).wtiMonthlyBetas.quarterlyChanges;
oilResponseForecast = varStruct(r).predictedCAR.Changes;

[anoms23, labels23] = getAnomalySignals('novyMarxVelikovAnomalies.csv', 1, 2);

% Double Sort

a  = nan(6, size(anoms23, 3));
tA = nan(6, size(anoms23, 3));

for i=1:size(anoms23, 3)
    ind1 = makeUnivSortInd(FillMonths(anoms23(:,:,i)), 5);
    ind = cond_index_from_inds(ind1, FillMonths(oilResponseForecast), 5);

    index = 1*(ind==1 | ind==6 | ind==11 | ind==16 | ind==21) + ... 
            2*(ind==2 | ind==7 | ind==12 | ind==17 | ind==22) + ...
            3*(ind==3 | ind==8 | ind==13 | ind==18 | ind==23) + ...
            4*(ind==4 | ind==9 | ind==14 | ind==19 | ind==24) + ...
            5*(ind==5 | ind==10 | ind==15 | ind==20 | ind==25);
    s = find(dates==197501);
    
    res = runUnivSort(ret, index, dates,me, 'plotFigure', 0, ...
                                            'printResults', 0, ...
                                            'factorModel', 1);

    pret = res.pret;

    indZeroOilQtrChange = find(oilChanges == 0 & ...
                          dates>197412);
    pret(indZeroOilQtrChange+1,end-1:end)=repmat(rf(indZeroOilQtrChange+1),1,2);
    pret(indZeroOilQtrChange+2,end-1:end)=repmat(rf(indZeroOilQtrChange+2),1,2);
    pret(indZeroOilQtrChange+3,end-1:end)=repmat(rf(indZeroOilQtrChange+3),1,2);
    pret(indZeroOilQtrChange+1,1:end-2)=0;
    pret(indZeroOilQtrChange+2,1:end-2)=0;
    pret(indZeroOilQtrChange+3,1:end-2)=0;

    for j=1:size(pret,2)
        if j==size(pret,2)
            tempRes=nanols(pret(s:end,j), const(s:end));
        else
            tempRes=nanols(pret(s:end,j)-rf(s:end),const(s:end));
        end
            
        a(j,i)=tempRes.beta(1);
        tA(j,i)=tempRes.tstat(1);
    end
    
    
end


a(end,:)=a(end,:);
tA(end,:)=tA(end,:);

clc
h=[{'(L)'},{'(2)'},{'(3)'},{'(4)'},{'(H)'},{'(L-H)'}];
mat2Tex(a,tA,h,2);



%% Table 5 - robustness to strategy construction

fprintf('\n\n\nTable 5 output:\n\n\n');

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

fprintf('\n\nPanel A:\n\n');

% Run the time-series regressions
ind = makeUnivSortInd(oilResponseForecast, nPtf);
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

printPtfResults(pret, s)


fprintf('\n\nPanel B:\n\n');

% Run the time-series regressions
ind = makeUnivSortInd(oilResponseForecast, nPtf, NYSE);
res = runUnivSort(ret, ind, dates, me, 'weighting', 'e', ...
                                       'plotFigure', 0, ...
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


printPtfResults(pret, s)




fprintf('\n\nPanel C:\n\n');

nPtf = 10;

% Run the time-series regressions
ind = makeUnivSortInd(oilResponseForecast, nPtf);
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


printPtfResults(pret, s)


%% Table 6 - Subsamples 

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

% Panel A: oil price changes:
% Calculate the thresholds for the oil price changes
quarterlyChanges = (oilChanges);
quarterlyChanges(dates<197501) = nan;
quarterlyChanges(isnan(quarterlyChanges)) = [];
thresholds = prctile(quarterlyChanges, [30 70]);

% Get the three samples indicators
lagChanges=lag((oilChanges),3,nan);
lagChanges(dates<197501)=nan;
for i=find(ismember(dates-100*floor(dates/100),[3 6 9 12]), 1, 'first'):3:length(lagChanges)
    lagChanges(i-1)=lagChanges(i);
    lagChanges(i-2)=lagChanges(i);
end

qChOne   = lagChanges < thresholds(1);
qChTwo   = lagChanges >=  thresholds(1) & ...
           lagChanges < thresholds(2);
qChThree = lagChanges >= thresholds(2);


% Get the average returns 
res1 = nanols(pret(qChOne,   end), 100*const(qChOne));
res2 = nanols(pret(qChTwo,   end), 100*const(qChTwo));
res3 = nanols(pret(qChThree, end), 100*const(qChThree));

fprintf('Panel A: different total oil price changes:\n')
a  = 100 * [res1.beta res2.beta res3.beta];
tA =       [res1.tstat res2.tstat res3.tstat];
mat2Tex(a, tA, {'$r^e$'});

% Panel B: unexpected oil changes
tempData = readtable('ExpectedOilKilianData.xlsx', 'Sheet','Data');
tempData = tempData(:,{'Date','QUnexpectedOilChange'});
tempData.Properties.VariableNames = {'dates','uOilChange'};

monthIndex = dates-100*floor(dates/100);
quarterEndIndex = (monthIndex==3) | (monthIndex==6) | (monthIndex==9) | (monthIndex==12);


[~,ia,ib] = intersect(dates,tempData.dates);
uOilQuarterly=nan(size(dates));
uOilQuarterly(ia) = tempData.uOilChange(ib);
uOilQuarterly(~quarterEndIndex) = nan;

% Calculate the thresholds for the oil price changes
quarterlyChanges = uOilQuarterly;
quarterlyChanges(dates<197501) = nan;
quarterlyChanges(isnan(quarterlyChanges)) = [];
thresholds = prctile(quarterlyChanges, [30 70]);

% Get the three samples indicators
lagChanges=lag((oilChanges),3,nan);
lagChanges(dates<197501)=nan;
for i=find(ismember(dates-100*floor(dates/100),[3 6 9 12]), 1, 'first'):3:length(lagChanges)
    lagChanges(i-1)=lagChanges(i);
    lagChanges(i-2)=lagChanges(i);
end

qChOne   = lagChanges < thresholds(1);
qChTwo   = lagChanges >=  thresholds(1) & ...
           lagChanges < thresholds(2);
qChThree = lagChanges >= thresholds(2);


% Get the average returns 
res1 = nanols(pret(qChOne,   end), 100*const(qChOne));
res2 = nanols(pret(qChTwo,   end), 100*const(qChTwo));
res3 = nanols(pret(qChThree, end), 100*const(qChThree));


fprintf('Panel B: different unexpected oil price changes:\n')
a  = 100 * [res1.beta res2.beta res3.beta];
tA =       [res1.tstat res2.tstat res3.tstat];
mat2Tex(a, tA, {'$r^e$'});

% Panel C: Get expected oil changes
tempData = readtable('ExpectedOilKilianData.xlsx', 'Sheet','Data');
tempData = tempData(:,{'Date','QUnexpectedOilChange'});
tempData.Properties.VariableNames = {'dates','uOilChange'};

monthIndex = dates-100*floor(dates/100);
quarterEndIndex = (monthIndex==3) | (monthIndex==6) | (monthIndex==9) | (monthIndex==12);

uOilMonthly=nan(size(dates));
uOilQuarterly=nan(size(dates));

[~,ia,ib] = intersect(dates,tempData.dates);
uOilQuarterly(ia) = tempData.uOilChange(ib);
uOilQuarterly(~quarterEndIndex) = nan;

% Calculate the thresholds for the oil price changes
quarterlyChanges = (oilChanges) - (uOilQuarterly);
quarterlyChanges(dates<197501) = nan;
quarterlyChanges(isnan(quarterlyChanges)) = [];
thresholds = prctile(quarterlyChanges, [30 70]);

% Get the three samples indicators
lagChanges=lag((oilChanges),3,nan);
lagChanges(dates<197501)=nan;
for i=find(ismember(dates-100*floor(dates/100),[3 6 9 12]), 1, 'first'):3:length(lagChanges)
    lagChanges(i-1)=lagChanges(i);
    lagChanges(i-2)=lagChanges(i);
end

qChOne   = lagChanges < thresholds(1);
qChTwo   = lagChanges >=  thresholds(1) & ...
           lagChanges < thresholds(2);
qChThree = lagChanges >= thresholds(2);


% Get the average returns 
res1 = nanols(pret(qChOne,   end), 100*const(qChOne));
res2 = nanols(pret(qChTwo,   end), 100*const(qChTwo));
res3 = nanols(pret(qChThree, end), 100*const(qChThree));


fprintf('Panel C: different expected oil price changes:\n')
a  = 100 * [res1.beta res2.beta res3.beta];
tA =       [res1.tstat res2.tstat res3.tstat];
mat2Tex(a, tA, {'$r^e$'});


% Get the GDP data
fredStruct = getFredData('USREC','1970-01-02',[], 'lin', [],'eop');
fredData = array2table(fredStruct.Data, 'VariableNames',[{'date'},{'recInd'}]);
fredData.date = datetime(fredData.date, 'ConvertFrom', 'datenum');
data = fredData;
data.dates=100*year(data.date)+month(data.date);
[~,ia,ib]=intersect(dates,data.dates);
recInd = nan(size(dates));
recInd(ia) = data.recInd(ib);

strat_ret = pret(:,end);
strat_ret(dates<197501) = nan;

ind = (recInd == 1);
res1 = (nanols(strat_ret(ind),const(ind)));
res2 = (nanols(strat_ret(~ind),const(~ind)));

fprintf('Panel D: recessions vs expansions:\n')
a  = [res1.beta res2.beta];
tA = [res1.tstat res2.tstat];
mat2Tex(a, tA, {'$r^e$'});

% Bear/bull
cumMkt = makePastPerformance(mkt, 23, 0)-1;
bear = 1* cumMkt<0;

strat_ret = pret(:,end);
strat_ret(dates<197501) = nan;

ind = (bear == 1);
res1 = (nanols(strat_ret(ind),const(ind)));
res2 = (nanols(strat_ret(~ind),const(~ind)));

fprintf('Panel E: bear vs bull:\n')
a  = [res1.beta res2.beta];
tA = [res1.tstat res2.tstat];
mat2Tex(a, tA, {'$r^e$'});

% Calculate the thresholds for the oil price changes
quarterlyChanges = (oilChanges);
quarterlyChanges(dates<197501) = nan;
quarterlyChanges(isnan(quarterlyChanges)) = [];
thresholds = prctile(quarterlyChanges, [30 70]);

% Get the three samples indicators
lagChanges=lag((oilChanges),3,nan);
lagChanges(dates<197501)=nan;
for i=find(ismember(dates-100*floor(dates/100),[3 6 9 12]), 1, 'first'):3:length(lagChanges)
    lagChanges(i-1)=lagChanges(i);
    lagChanges(i-2)=lagChanges(i);
end

qChOne   = lagChanges < thresholds(1);
qChTwo   = lagChanges >=  thresholds(1) & ...
           lagChanges < thresholds(2);
qChThree = lagChanges >= thresholds(2);



% Get the GDP data
fredStruct = getFredData('USREC','1970-01-02',[], 'lin', [],'eop');
fredData = array2table(fredStruct.Data, 'VariableNames',[{'date'},{'recInd'}]);
fredData.date = datetime(fredData.date, 'ConvertFrom', 'datenum');
data = fredData;
data.dates=100*year(data.date)+month(data.date);
[~,ia,ib]=intersect(dates,data.dates);
recInd = nan(size(dates));
recInd(ia) = data.recInd(ib);

strat_ret = pret(:,end);
strat_ret(dates<197501) = nan;

ind = qChOne & (recInd == 0);
res1 = (nanols(strat_ret(ind),const(ind)));
ind = qChTwo & (recInd == 0);
res2 = (nanols(strat_ret(ind),const(ind)));
ind = qChThree & (recInd == 0);
res3 = (nanols(strat_ret(ind),const(ind)));

ind = qChOne & (recInd == 1);
res4 = (nanols(strat_ret(ind),const(ind)));
ind = qChTwo & (recInd == 1);
res5 = (nanols(strat_ret(ind),const(ind)));
ind = qChThree & (recInd == 1);
res6 = (nanols(strat_ret(ind),const(ind)));

fprintf('Panel F: double sort (returns):\n')
a  = [res1.beta res2.beta res3.beta; 
      res4.beta res5.beta res6.beta];
tA  = [res1.tstat res2.tstat res3.tstat; 
      res4.tstat res5.tstat res6.tstat];
nobs = [res1.nobs res2.nobs res3.nobs; 
      res4.nobs res5.nobs res6.nobs];

h = {'Expansions','Recessions'};

mat2Tex(a,tA, h, 2);

fprintf('Panel G: double sort (# stocks):\n')
mat2Tex(nobs,nobs,h,0);



%% Table 6 - robustness to industry classification

fprintf('\n\n\nTable 6 output:\n\n\n');

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

fprintf('\n\nPanel A:\n\n');

h = cell(1,1);
a = [];
tA = [];

for r = 1:length(varStruct)
    % Find the FF49 results
%     r = find(strcmp([varStruct.label],{'FF49'}));
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
    for j = 1:(nPtf+1)
        res = nanols(pret(:, j), const);
        a(r, j) = res.beta;
        tA(r, j) = res.tstat;
    end
    h(r,1) = cellstr(['$r_{',char(varStruct(r).label),'}^e$']);    
end

mat2Tex(a, tA, h, 2);

fprintf('\n\nPanel B:\n\n');

a = [];
tA = [];
h = cell(1,1);

for r = 1:length(varStruct)
    % Find the FF49 results
%     r = find(strcmp([varStruct.label],{'FF49'}));
    oilChanges = varStruct(r).wtiMonthlyBetas.quarterlyChanges;
    oilResponseForecast = varStruct(r).predictedCAR.Returns;


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

    for j = 1:(nPtf+1)
        res = nanols(pret(:, j), const);
        a(r, j) = res.beta;
        tA(r, j) = res.tstat;
    end
    h(r,1) = cellstr(['$r_{',char(varStruct(r).label),'}^e$']);
 
end
mat2Tex(a, tA, h, 2);



%% Table 7 - announcing vs non-announcing firms
 
fprintf('\n\n\nTable 7 output:\n\n\n');

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
nPtf = 4;


% Find the FF49 results
r = find(strcmp([varStruct.label],{'FF49'}));
oilChanges = varStruct(r).wtiMonthlyBetas.quarterlyChanges;
oilResponseForecast = varStruct(r).predictedCAR.Changes;


% Get the expected announcers
load IBQ

amonth = IBQ ~= lag(IBQ,1,nan);
amonth = isfinite(IBQ).*amonth;
indAMO = isfinite(IBQ) + lag(amonth,11,nan); % firms that announced eleven months ago are 
                                         % likely to announce in the comming month
% res4 = runUnivSort(ret, indAMO, dates ,me, 'timePeriod', 197412, 'plotFigure', 0);
     


% Do the non-announcers first
temp = FillMonths(oilResponseForecast);
temp(indAMO == 2) = nan;
ind = makeUnivSortInd(temp, nPtf);
res1 = runUnivSort(ret, ind, dates,me, 'plotFigure', 0, ...
                                       'printResults', 0, ...
                                       'factorModel', 1);

% Account for the zero-oil-price-change quarters
pret = res1.pret;

% Account for portfolio-months with fewer than 20 stocks
for i = 1:nPtf-1
    ind = res1.ptfNumStocks(:,i)<20;
    pret(ind,i) = 0;
end
ind = res1.ptfNumStocks(:,nPtf)<20;
pret(ind, nPtf) = rf(ind);
pret(:,nPtf+1) = pret(:,nPtf) - pret(:,1);

indZeroOilQtrChange = find(oilChanges == 0 & ...
                          dates>197412);
pret(indZeroOilQtrChange+1,end-1:end)=repmat(rf(indZeroOilQtrChange+1),1,2);
pret(indZeroOilQtrChange+2,end-1:end)=repmat(rf(indZeroOilQtrChange+2),1,2);
pret(indZeroOilQtrChange+3,end-1:end)=repmat(rf(indZeroOilQtrChange+3),1,2);
pret(indZeroOilQtrChange+1,1:end-2)=0;
pret(indZeroOilQtrChange+2,1:end-2)=0;
pret(indZeroOilQtrChange+3,1:end-2)=0;

fprintf('\n\nPanel A: Non-Announcers.\n\n');
printPtfResults(pret, s)


% Do the announcers next
temp = FillMonths(oilResponseForecast);
temp(indAMO == 1) = nan;
ind = makeUnivSortInd(temp, nPtf);
res2 = runUnivSort(ret, ind, dates,me, 'plotFigure', 0, ...
                                       'printResults', 0, ...
                                       'factorModel', 1);

pret = res2.pret;

% Account for portfolio-months with fewer than 20 stocks
for i = 1:nPtf-1
    ind = res2.ptfNumStocks(:,i)<20;
    pret(ind,i) = 0;
end
ind = res2.ptfNumStocks(:,nPtf)<20;
pret(ind, nPtf) = rf(ind);
pret(:,nPtf+1) = pret(:,nPtf) - pret(:,1);

% Account for the zero-oil-price-change quarters
indZeroOilQtrChange = find(oilChanges == 0 & ...
                          dates>197412);
pret(indZeroOilQtrChange+1,end-1:end)=repmat(rf(indZeroOilQtrChange+1),1,2);
pret(indZeroOilQtrChange+2,end-1:end)=repmat(rf(indZeroOilQtrChange+2),1,2);
pret(indZeroOilQtrChange+3,end-1:end)=repmat(rf(indZeroOilQtrChange+3),1,2);
pret(indZeroOilQtrChange+1,1:end-2)=0;
pret(indZeroOilQtrChange+2,1:end-2)=0;
pret(indZeroOilQtrChange+3,1:end-2)=0;

fprintf('\n\nPanel B: Announcers.\n\n');
printPtfResults(pret, s)




%% Table 8 - placebo predictability
 
fprintf('\n\n\nTable 8 output:\n\n\n');

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
oilChanges = varStruct(r).wtiMonthlyBetas.Changes;


% Run the time-series regressions
ind = makeUnivSortInd(oilChanges, nPtf, NYSE);
res = runUnivSort(ret, ind, dates, me, 'plotFigure', 0, ...
                                       'printResults', 0, ...
                                       'factorModel', 1);

% Account for the zero-oil-price-change quarters
pret = res.pret;

printPtfResults(pret, s, 1)

%% Table 9 - commodities/macro
 
fprintf('\n\n\nTable 9 output:\n\n\n');

warning('off','all')

clear
clc

load ret
load me
load dates
load nyse
load varStruct
load ff


% Find the FF49 results
r = find(strcmp([varStruct.label],{'FF49'}));
wCAR3 = varStruct(r).wCAR3;

fredKeys = {'WPU101','WPU01830131','WPU01220205','WPU10210501','WPUSI019011','WPU102402', ...
            'GDP','IPMAN','UNRATE','PCUOMFGOMFG','CPALTT01USM657N', ...
            'GS10', 'GS2', 'VIXCLS','LIOR3M'};
        
fredStruct = getFredData(char(fredKeys(1)),'1970-01-02',[], 'lin', [],'eop');
fredData = array2table(fredStruct.Data, 'VariableNames',[{'date'},fredKeys(1)]);
fredData.date = datetime(fredData.date, 'ConvertFrom', 'datenum');
fredData.dates = 100*year(fredData.date)+month(fredData.date);
fredData = fredData(:, {'date','dates',char(fredKeys(1))});

for i = 2:length(fredKeys)
    try
        if ismember(fredKeys(i), {'GDP','LIOR3M'})            
            fredStruct = getFredData(char(fredKeys(i)),'1970-01-02',[],'lin','q','eop');
            tempData = array2table(fredStruct.Data, 'VariableNames',[{'date'},fredKeys(i)]);
            tempData.date = datetime(tempData.date, 'ConvertFrom', 'datenum');
            tempData.date = tempData.date + calmonths(2);
        else
            fredStruct = getFredData(char(fredKeys(i)),'1970-01-02',[],'lin','m','eop');
            tempData = array2table(fredStruct.Data, 'VariableNames',[{'date'},fredKeys(i)]);
            tempData.date = datetime(tempData.date, 'ConvertFrom', 'datenum');
        end
        fredData = outerjoin(fredData, tempData, 'MergeKeys', 1);
    catch
        fprintf('Error with %s.\n', char(fredKeys(i)))
    end
end
data = fredData;

% Add slopes
data.slope = data.GS10-data.GS2;
data.GS10 = [];

% Add unexpected SPY
load uSPY
tempData = array2table([dates uSPY], 'VariableNames', {'dates','uSPY'});
% tempData = readtable('SPYData.xlsx', 'Sheet','UnexpectedSPY');
% tempData.Properties.VariableNames = {'dates','uSPY'};
data = outerjoin(data, tempData, 'MergeKeys', 1, 'Type', 'Left');


% Make the structure
monthIndex=dates-100*floor(dates/100);
quarterEndIndex=(monthIndex==3) | (monthIndex==6) | (monthIndex==9) | (monthIndex==12);

[~,ia,ib] = intersect(dates,data.dates);

commodityStruct=struct;
commodityStruct(1).monthly=nan(size(dates));
commodityStruct(1).quarterly=nan(size(dates));
commodityStruct(1).quarterlyChanges=nan(size(dates));
commodityStruct(1).quarterlyReturns=nan(size(dates));

varNames=data.Properties.VariableNames';
for i=2:length(varNames)-1
    eval(['commodityStruct(i-1,1).monthly(ia,1)=data.',char(varNames(i+1)),'(ib);']);
    commodityStruct(i-1,1).quarterly = commodityStruct(i-1).monthly;
    commodityStruct(i-1,1).quarterly(~quarterEndIndex) = nan;
    commodityStruct(i-1,1).quarterlyChanges = commodityStruct(i-1).quarterly-lag(commodityStruct(i-1).quarterly,3,nan);
    commodityStruct(i-1,1).quarterlyReturns = commodityStruct(i-1).quarterly./lag(commodityStruct(i-1).quarterly,3,nan)-1;    
    commodityStruct(i-1,1).label=(varNames(i+1));
end

% Make the variables we don't want changes in unexpected market 
indLevels = find(ismember([commodityStruct.label], {'uSPY'}));
for i = indLevels
    commodityStruct(i).quarterlyChanges = commodityStruct(i).quarterly;
    commodityStruct(i).quarterlyReturns = commodityStruct(i).quarterly;
end

horizon=12 * 3; % 20 quarters enough?
min_quarters=12; % 12 quarters minimum enough?

for i=1:length(commodityStruct)
    commodityBetasOLS(i,1)=makeCommodityBetas(wCAR3,commodityStruct(i),horizon,min_quarters);
end

% predictedCAR is a structure with all 4 betas - changes, returns
for i=1:length(commodityBetasOLS)
    predictedCAR(i,1)=makePredictedCAR(commodityBetasOLS(i)); 
end


for i=1:length(predictedCAR)
    try
        ind = makeUnivSortInd(predictedCAR(i).Changes, 5, NYSE);
        res(i,1) = runUnivSort(ret, ind, dates, me, 'factorModel', 4, ...
                                                    'plotFigure', 0, ...
                                                    'printResults', 0);        
    catch
        res(i, :) = res(i-1, :);
    end
end


mapping = [{'WPU101','WPU01830131','WPU01220205','WPU10210501','WPUSI019011','WPU102402', ...
            'GDP','IPMAN','UNRATE','PCUOMFGOMFG','CPALTT01USM657N', ...
            'GS2', 'VIXCLS','LIOR3M','uSPY','slope'};
            {'IronSteel','Soybeans','Corn','Gold','Copper','Aluminum', ...
            'GDP','IndProd','UnRate','PPI','CPI', ...
            'TreasuryLevel','VIX','LIBOR','UnexpectedSPX','TreasurySlope'}]';        

for i=1:size(mapping,1)
    r = find(ismember([commodityStruct.label], mapping(i,1)));
    if isfinite(r)
        commodityStruct(r).label = mapping(i,2);         
    end
end

clc

for i=1:length(commodityStruct)
    pret = res(i).pret;
    indZeroOilQtrChange = find(commodityStruct(i).quarterlyChanges == 0 & ...
                              dates>197412 & dates<dates(end-2));
    pret(indZeroOilQtrChange+1,end-1:end)=repmat(rf(indZeroOilQtrChange+1),1,2);
    pret(indZeroOilQtrChange+2,end-1:end)=repmat(rf(indZeroOilQtrChange+2),1,2);
    pret(indZeroOilQtrChange+3,end-1:end)=repmat(rf(indZeroOilQtrChange+3),1,2);
    pret(indZeroOilQtrChange+1,1:end-2)=0;
    pret(indZeroOilQtrChange+2,1:end-2)=0;
    pret(indZeroOilQtrChange+3,1:end-2)=0;
    
    pret(dates<197501,:)=nan;
    for j = 1:size(pret,2)
        if j == size(pret,2)
            tempRes = nanols(pret(:,j), const);
        else
            tempRes = nanols(pret(:,j) - rf, const);
        end
        res(i).xret(j) = tempRes.beta;
        res(i).txret(j) = tempRes.tstat;
    end
    h(i,1)={(['$r^{\text{',char(commodityStruct(i).label),'}}$'])};
    a(i,:)=res(i).xret;
    tA(i,:)=res(i).txret;
end

fprintf('\n\nPanel A:');
mat2Tex(a(1:6, :),tA(1:6, :),h(1:6, :),2);
fprintf('\n\nPanel B:');
mat2Tex(a(7:end, :),tA(7:end, :),h(7:end, :),2);

%% Table 10 - share turnover & volume table

fprintf('\n\n\nTable 10 output:\n\n\n');

fprintf('\n\nPanel A:');

clear
clc

load ret
load me
load dates
load nyse
load varStruct
load vol
load shrout

% Determine the start date & number of portfolios
s = find(dates==197501);
nPtf = 5;

% Find the FF49 results
r = find(strcmp([varStruct.label],{'FF49'}));
oilChanges = varStruct(r).wtiMonthlyBetas.quarterlyChanges;
oilResponseForecast = varStruct(r).predictedCAR.Changes;

% Run the time-series regressions
ind = makeUnivSortInd(oilResponseForecast, nPtf, NYSE);

indTO=zeros(size(ret));
for i=s-1:3:length(dates)-3
    indTO(i+1:i+3,:)=repmat(ind(i,:),3,1);
end

vol(vol<=0)=nan;
shrout(shrout<=0)=nan;

vol = vol*100;
shrout = shrout*1000;
turnover = log(1+vol./shrout);
turnover(1:s-1,:)=nan;
turnover(abs(turnover)==Inf)=nan;

ts_TO = nan(length(dates),5);
for i=1:5
    tempTO = turnover;
    tempME = lag(me,1,nan);
    index = (indTO==i) & isfinite(tempME);
    tempTO(~index) = nan;
    tempME(~index) = nan;
    ts_TO(:,i)=nansum(tempTO.*tempME./(repmat(nansum(tempME,2),1,size(ret,2))),2);    
end

indZeroOilQtrChange = find(oilChanges == 0 & ...
                          dates>197412);

ind1=find(dates>197412 & ismember(dates-100*floor(dates/100), [1 4 7 10]));
ind2=find(dates>197412 & ismember(dates-100*floor(dates/100), [2 5 8 11]));
ind3=find(dates>197412 & ismember(dates-100*floor(dates/100), [3 6 9 12]));
ind1(ismember(ind1, indZeroOilQtrChange+1)) = [];
ind2(ismember(ind2, indZeroOilQtrChange+2)) = [];
ind3(ismember(ind3, indZeroOilQtrChange+3)) = [];

a=100*[nanmean(ts_TO(ind1,:),1); ...
    nanmean(ts_TO(ind2,:),1); ...
    nanmean(ts_TO(ind3,:),1)];

h={'Month 1','Month 2','Month 3'};
mat2Tex(a,a,h,2);

fprintf('\n\nPanel B:');

%%
clear

load ret
load me
load dates
load nyse
load weightedTurnover
load varStruct

% Determine the start date & number of portfolios
nPtf = 5;


% Find the FF49 results
r = find(strcmp([varStruct.label],{'FF49'}));
oilResponseForecast = varStruct(r).predictedCAR.Changes;
s = find(dates==197501);

% Run the time-series regressions
ind = makeUnivSortInd(oilResponseForecast, nPtf, NYSE);

indFinite=find(sum(ind,2)>0);
index=zeros(size(ret));
for i=1:length(indFinite)
    index(indFinite(i)+3,:)=ind(indFinite(i),:);        
end

for i=1:5
    temp=weightedTurnover;
    temp(index~=i)=nan;    
    a(1,i)=nanmean(nanmean(temp(s:end,:),2));
end

h={'wTO'};
mat2Tex(100*a,100*a,h,2);

%% Timekeeping


fprintf('\n\n\n\nDone with table printing @ %s\n\n\n',char(datetime('now')));


