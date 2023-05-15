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

% %% Make subsample performances
% 
% clear
% clc
% 
% load ret
% load me
% load dates
% load nyse
% load ff
% load varStruct
% 
% % Determine the start date & number of portfolios
% s = find(dates==197501);
% nPtf = 5;
% 
% 
% % Find the FF49 results
% r = find(strcmp([varStruct.label],{'FF49'}));
% oilChanges = varStruct(r).wtiMonthlyBetas.quarterlyChanges;
% oilResponseForecast = varStruct(r).predictedCAR.Changes;
% 
% 
% % Run the time-series regressions
% ind = makeUnivSortInd(oilResponseForecast, nPtf, NYSE);
% res = runUnivSort(ret, ind, dates, me, 'plotFigure', 0, ...
%                                        'printResults', 0, ...
%                                        'factorModel', 1);
% 
% % Account for the zero-oil-price-change quarters
% pret = res.pret;
% indZeroOilQtrChange = find(oilChanges == 0 & ...
%                           dates>197412);
% pret(indZeroOilQtrChange+1,end-1:end)=repmat(rf(indZeroOilQtrChange+1),1,2);
% pret(indZeroOilQtrChange+2,end-1:end)=repmat(rf(indZeroOilQtrChange+2),1,2);
% pret(indZeroOilQtrChange+3,end-1:end)=repmat(rf(indZeroOilQtrChange+3),1,2);
% pret(indZeroOilQtrChange+1,1:end-2)=0;
% pret(indZeroOilQtrChange+2,1:end-2)=0;
% pret(indZeroOilQtrChange+3,1:end-2)=0;
% 
% % Panel A: oil price changes:
% % Calculate the thresholds for the oil price changes
% quarterlyChanges = (oilChanges);
% quarterlyChanges(dates<197501) = nan;
% quarterlyChanges(isnan(quarterlyChanges)) = [];
% thresholds = prctile(quarterlyChanges, [30 70]);
% 
% % Get the three samples indicators
% lagChanges=lag((oilChanges),3,nan);
% lagChanges(dates<197501)=nan;
% for i=find(ismember(dates-100*floor(dates/100),[3 6 9 12]), 1, 'first'):3:length(lagChanges)
%     lagChanges(i-1)=lagChanges(i);
%     lagChanges(i-2)=lagChanges(i);
% end
% 
% qChOne   = lagChanges < thresholds(1);
% qChTwo   = lagChanges >=  thresholds(1) & ...
%            lagChanges < thresholds(2);
% qChThree = lagChanges >= thresholds(2);
% 
% 
% % Get the average returns 
% res1 = nanols(pret(qChOne,   end), 100*const(qChOne));
% res2 = nanols(pret(qChTwo,   end), 100*const(qChTwo));
% res3 = nanols(pret(qChThree, end), 100*const(qChThree));
% 
% fprintf('Panel A: different total oil price changes:\n')
% a  = 100 * [res1.beta res2.beta res3.beta];
% tA =       [res1.tstat res2.tstat res3.tstat];
% mat2Tex(a, tA, {'$r^e$'});
% 
% % Get unexpected oil changes
% 
% load macroOtherCommodityStruct
% 
% % Calculate the thresholds for the oil price changes
% r = find(strcmp([commodityStruct.label],{'uOilChange'}));
% quarterlyChanges = commodityStruct(r).quarterlyChanges;
% quarterlyChanges(dates<197501) = nan;
% quarterlyChanges(isnan(quarterlyChanges)) = [];
% thresholds = prctile(quarterlyChanges, [30 70]);
% 
% % Get the three samples indicators
% lagChanges=lag((oilChanges),3,nan);
% lagChanges(dates<197501)=nan;
% for i=find(ismember(dates-100*floor(dates/100),[3 6 9 12]), 1, 'first'):3:length(lagChanges)
%     lagChanges(i-1)=lagChanges(i);
%     lagChanges(i-2)=lagChanges(i);
% end
% 
% qChOne   = lagChanges < thresholds(1);
% qChTwo   = lagChanges >=  thresholds(1) & ...
%            lagChanges < thresholds(2);
% qChThree = lagChanges >= thresholds(2);
% 
% 
% % Get the average returns 
% res1 = nanols(pret(qChOne,   end), 100*const(qChOne));
% res2 = nanols(pret(qChTwo,   end), 100*const(qChTwo));
% res3 = nanols(pret(qChThree, end), 100*const(qChThree));
% 
% 
% fprintf('Panel B: different unexpected oil price changes:\n')
% a  = 100 * [res1.beta res2.beta res3.beta];
% tA =       [res1.tstat res2.tstat res3.tstat];
% mat2Tex(a, tA, {'$r^e$'});
% 
% % Get expected oil changes
% 
% load macroOtherCommodityStruct
% 
% % Calculate the thresholds for the oil price changes
% r = find(strcmp([commodityStruct.label],{'uOilChange'}));
% quarterlyChanges = (oilChanges)-commodityStruct(r).quarterlyChanges;
% quarterlyChanges(dates<197501) = nan;
% quarterlyChanges(isnan(quarterlyChanges)) = [];
% thresholds = prctile(quarterlyChanges, [30 70]);
% 
% % Get the three samples indicators
% lagChanges=lag((oilChanges),3,nan);
% lagChanges(dates<197501)=nan;
% for i=find(ismember(dates-100*floor(dates/100),[3 6 9 12]), 1, 'first'):3:length(lagChanges)
%     lagChanges(i-1)=lagChanges(i);
%     lagChanges(i-2)=lagChanges(i);
% end
% 
% qChOne   = lagChanges < thresholds(1);
% qChTwo   = lagChanges >=  thresholds(1) & ...
%            lagChanges < thresholds(2);
% qChThree = lagChanges >= thresholds(2);
% 
% 
% % Get the average returns 
% res1 = nanols(pret(qChOne,   end), 100*const(qChOne));
% res2 = nanols(pret(qChTwo,   end), 100*const(qChTwo));
% res3 = nanols(pret(qChThree, end), 100*const(qChThree));
% 
% 
% fprintf('Panel C: different expected oil price changes:\n')
% a  = 100 * [res1.beta res2.beta res3.beta];
% tA =       [res1.tstat res2.tstat res3.tstat];
% mat2Tex(a, tA, {'$r^e$'});
% 
% 
% % Get the GDP data
% fredStruct = getFredData('USREC','1970-01-02','2021-12-31', 'lin', [],'eop');
% fredData = array2table(fredStruct.Data, 'VariableNames',[{'date'},{'recInd'}]);
% fredData.date = datetime(fredData.date, 'ConvertFrom', 'datenum');
% data = fredData;
% data.dates=100*year(data.date)+month(data.date);
% [~,ia,ib]=intersect(dates,data.dates);
% recInd = nan(size(dates));
% recInd(ia) = data.recInd(ib);
% 
% strat_ret = pret(:,end);
% strat_ret(dates<197501) = nan;
% 
% ind = (recInd == 1);
% res1 = (nanols(strat_ret(ind),const(ind)));
% res2 = (nanols(strat_ret(~ind),const(~ind)));
% 
% fprintf('Panel D: recessions vs expansions:\n')
% a  = [res1.beta res2.beta];
% tA = [res1.tstat res2.tstat];
% mat2Tex(a, tA, {'$r^e$'});
% 
% % Bear/bull
% cumMkt = makePastPerformance(mkt, 23, 0)-1;
% bear = 1* cumMkt<0;
% 
% strat_ret = pret(:,end);
% strat_ret(dates<197501) = nan;
% 
% ind = (bear == 1);
% res1 = (nanols(strat_ret(ind),const(ind)));
% res2 = (nanols(strat_ret(~ind),const(~ind)));
% 
% fprintf('Panel E: bear vs bull:\n')
% a  = [res1.beta res2.beta];
% tA = [res1.tstat res2.tstat];
% mat2Tex(a, tA, {'$r^e$'});

%% Make uSPY

clear
clc

% Input your WRDS username
Params.username = usernameUI();                                            

% Input your WRDS password
Params.pass     = passwordUI();                                            

% Call WRDS connection
WRDS = callWRDSConnection(Params.username,Params.pass);

% Download Tesla's stock price
ticker = 'SPY';
qry = ['select date, ret from CRSP.MSF left join CRSP.MSFHDR on msf.permno=msfhdr.permno where htick=''', ticker, ''''];
retTable = fetch(WRDS, qry);
retTable.date = datetime(retTable.date);
retTable.date = 100*year(retTable.date) + month(retTable.date);

load dates
load ff

[~, ia, ib] = intersect(dates, retTable.date);

spy = nan(size(dates));
spy(ia) = retTable.ret(ib);

% STR factor
unzip('https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/F-F_ST_Reversal_Factor_CSV.zip', 'Data');
opts =detectImportOptions('F-F_ST_Reversal_Factor.CSV');
data =readtable('F-F_ST_Reversal_Factor.CSV', opts);
[~, ia, ib] = intersect(dates, data.Var1);
STR = nan(size(dates));
STR(ia) = data.ST_Rev(ib);


% LTR factor
unzip('https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/F-F_LT_Reversal_Factor_CSV.zip', 'Data');
opts =detectImportOptions('F-F_LT_Reversal_Factor.CSV');
data =readtable('F-F_LT_Reversal_Factor.CSV', opts);
[~, ia, ib] = intersect(dates, data.Var1);
LTR = nan(size(dates));
LTR(ia) = data.LT_Rev(ib);

ind = isfinite(spy);
x = [const smb hml umd STR LTR];
res = (ols(spy(ind),x(ind,:)));

uSPY = nan(size(dates));
uSPY(ind) = res.resid;

save Data\uSPY uSPY

