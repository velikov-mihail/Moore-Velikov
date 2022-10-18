%% Timekeeping

fprintf('Now working on the making the figures. Run started @ %s.\n\n\n',char(datetime('now')));

%% Figure 1 - Two-firm figure

clear
clc

load permno
load dff
load dprc
load ddates

c1 = find(permno==89790);
c2 = find(permno==89014);

lightGray  = [0.66,0.66,0.66];
darkGray   = [0.33,0.33,0.3];
lineWidth  = 3;
markerSize = 10;
titleSize  = 25;
axisSize   = 25;

r = find(floor(ddates/100)==201410,1)-1;
r2 = 62;

mkt=1+dmkt(r-63:r+63);
mkt(1)=1;
mkt=cumprod(mkt);
mkt=mkt/mkt(64);        


x=(r-63:r+63);
y1=dprc(r-63:r+63,c1)/dprc(r,c1);
y2=dprc(r-63:r+63,c2)/dprc(r,c2);

figure('visible','off');
b = plot(x,100*[y1 y2 mkt]);

set(b(1),'LineStyle','-','Marker','s','color', 'k','LineWidth',lineWidth,'MarkerSize',markerSize)
set(b(2),'LineStyle','-','Marker','x','color',darkGray,'LineWidth',lineWidth,'MarkerSize',markerSize)
set(b(3),'LineStyle',':','color',lightGray,'LineWidth',lineWidth,'MarkerSize',markerSize)
leg=legend('MPET','COL','MKT','Location','south','Orientation','horizontal');
legend('boxoff');
xlim(x([1 end]));
ax=gca;
ax.XTick=find(ismember(ddates,[[20140701],[20140731],[20140829],[20140930],[20141031],[20141128],[20141231]]));
ax.XTickLabel=[{'01/Jul/14'},{'31/Jul/14'},{'29/Aug/14'},{'30/Sep/14'},{'31/Oct/14'},{'28/Nov/14'},{'31/Dec/14'}];
title(['Daily prices of MPET and COL'],'FontSize',titleSize);
set(gca,'FontSize',axisSize');
set(gcf, 'PaperPositionMode', 'auto');
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
set(gca,'LooseInset',get(gca,'TightInset'))

% Add the textboxes - Chevron and Exxon
textBoxText = 'Chevron and Exxon announce on 31/Oct/2014';
t = annotation('textbox', [0.3, 0.3, 0.1, 0.1], 'String', textBoxText);
t.FontSize = 16;
t.EdgeColor = [1 1 1];
annotation('arrow', [0.54, 0.693], [0.38, 0.65]);

% Add the textboxes - Magellan Petroleum
textBoxText = 'Magellan Petroleum announces on 13/Nov/2014';
t = annotation('textbox', [0.35, 0.10, 0.1, 0.1], 'String', textBoxText);
t.FontSize = 16;
t.EdgeColor = [1 1 1];
annotation('arrow', [0.59, 0.75], [0.2, 0.51]);

% Add the textboxes - Boeing
textBoxText = 'Boeing announces on 22/Oct/2014';
t = annotation('textbox', [0.66, 0.62, 0.1, 0.1], 'String', textBoxText);
t.FontSize = 16;
t.Color = 0.5 * [1 1 1];
t.EdgeColor = [1 1 1];
t = annotation('arrow', [0.66, 0.64], [0.71, 0.81]);
t.Color = 0.5 * [1 1 1];

% Add the textboxes - Rockwell Collins
textBoxText = 'Rockwell Collins announces on 31/Oct/2014';
t = annotation('textbox', [0.68, 0.67, 0.1, 0.1], 'String', textBoxText);
t.FontSize = 16;
t.Color = 0.5 * [1 1 1];
t.EdgeColor = [1 1 1];
t = annotation('arrow', [0.70, 0.693], [0.77, 0.90]);
t.Color = 0.5 * [1 1 1];

export_fig('Figures/figure1.pdf','-transparent');

%% Figure 4: Performance by industry

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

indTRet = nan(49, 3);
indMean = nan(49, 1);
indKurt = nan(49, 1);

for i = 1:49
    temp = oilChanges;
    temp(FF49~=i) = nan;
    indMean(i) = nanmean(nanmedian(temp(s:end,:),2));
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
    indTRet(i,:) = res(i).txret;    
end

%%

figure('visible','off');

subplot(2,1,1)
x = 100*indMean;
y = indTRet(:,3);

h1 = scatter(x, y, 'filled');
xlabel('Average Industry \beta');
ylabel({'T-statistic on strategy return','within industry'});
h2 = lsline;
res = nanols(y,[ones(size(y)) x]);
txt = sprintf('y = %.2f + %.2fx',res.beta(1),res.beta(2));
text(0.04,1,txt,...
        'HorizontalAlignment','center',...
        'VerticalAlignment','bottom',...
        'FontSize',12);
text(0.037,0.6,strcat('[',num2str(res.tstat(1),'%0.2f'),']'),...
        'HorizontalAlignment','center',...
        'VerticalAlignment','bottom',...
        'FontSize',12);
text(0.046,0.6,strcat('[',num2str(res.tstat(2),'%0.2f'),']'),...
        'HorizontalAlignment','center',...
        'VerticalAlignment','bottom',...
        'FontSize',12);


subplot(2,1,2)
x = indKurt;
y = indTRet(:,3);

h1 = scatter(x, y, 'filled');
xlabel('Kurtosis of Industry \beta');
ylabel({'T-statistic on strategy return','within industry'});

h2 = lsline;
res = nanols(y,[ones(size(y)) x]);
txt = sprintf('y = %.2f + %.2fx',res.beta(1),res.beta(2));
text(8,0,1,txt,...
        'HorizontalAlignment','center',...
        'VerticalAlignment','bottom',...
        'FontSize',12);
text(7.85,-0.4,strcat('[',num2str(res.tstat(1),'%0.2f'),']'),...
        'HorizontalAlignment','center',...
        'VerticalAlignment','bottom',...
        'FontSize',12);
text(8.5,-0.4,strcat('[',num2str(res.tstat(2),'%0.2f'),']'),...
        'HorizontalAlignment','center',...
        'VerticalAlignment','bottom',...
        'FontSize',12);

set(gcf, 'PaperPositionMode', 'auto');
set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
set(gca,'LooseInset',get(gca,'TightInset'))

export_fig('Figures/figure4.pdf','-transparent');


%% Figure 5: Strategy returns following large oil price changes

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
                                       'printResults', 1, ...
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


a  = 100 * [res1.beta res2.beta res3.beta];
tA =       [res1.tstat res2.tstat res3.tstat];

figure('visible','off');
bar(a,'FaceColor',[0 0 0.6]);
xlim([0.5 3.5]);
titleSize=25;
axisSize=25;
legendSize=20;


title('Average returns to strategy following quarterly oil price changes','FontSize',titleSize);
ylabel('xret in %/month','FontSize',axisSize);
xlabel('\Delta P_{t-1}^{Oil}','FontSize',axisSize);
set(gca,'XTickLabel',[{'Tertile 1'},{'Tertile 2'},{'Tertile 3'}]);
set(gca,'FontSize',axisSize');

for i=1:length(a)
    text(i,a(i)+0.07,num2str(a(i),'%0.2f'),...
            'HorizontalAlignment','center',...
            'VerticalAlignment','bottom',...
            'FontSize',axisSize);
    text(i,a(i)+0.03,strcat('[',num2str(tA(i),'%0.2f'),']'),...
            'HorizontalAlignment','center',...
            'VerticalAlignment','bottom',...
            'FontSize',axisSize);
end
set(gca,'Ylim',[0 max(a)+0.2]);
set(gcf, 'PaperPositionMode', 'auto');
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
set(gca,'LooseInset',get(gca,'TightInset'))


export_fig('Figures/figure5.pdf','-transparent');

%% Figure 6: Event-time quarterly returns relative to quarter start figure

clear
clc

load ret
load dates
load ddates
load me
load ptf_drets
load permno

eoqflag=eomflag==1 & ismember(floor(ddates/100)-100*floor(ddates/10000),[3 6 9 12]);
Qends=find(eoqflag==1);
qlength=(Qends-lag(Qends,1,nan));

qtrRelRetsShort=nan(nanmax(qlength),nanmax(qlength));
qtrRelRetsLong=nan(nanmax(qlength),nanmax(qlength));

for i=1:length(qlength)-1
    s=Qends(i)+1;
    e=Qends(i+1);
%     ddates([s e])'
    qtrRelRetsShort(i,1:(e-s+1))=ptf_drets(s:e,1)';
    qtrRelRetsLong(i,1:(e-s+1))=ptf_drets(s:e,5)';
end

short_ret=100*([1 cumprod(1+nanmean(qtrRelRetsShort,1))]-1);
long_ret=100*([1 cumprod(1+nanmean(qtrRelRetsLong,1))]-1);

temp=qtrRelRetsLong-qtrRelRetsShort;
temp(isnan(sum(temp,2)),:)=[];
temp2=[cumprod(1+temp')];
long_short=100*[(nanmean(temp2',1)-1)]';
b1=100*[1.645*(nanstd(temp2'-1,[],1))/sqrt(size(temp2,2))]';
b2=b1;
% long_short=100*([cumprod(1+nanmean(temp,1))]-1);
% b1=[100*([cumprod(1+nanmean(temp,1)-nanstd(temp,[],1)/sqrt(size(temp,1)))]-1)];
% b2=[100*([cumprod(1+nanmean(temp,1)+nanstd(temp,[],1)/sqrt(size(temp,1)))]-1)];

titleSize=25;
axisSize=25;
legendSize=20;
lineWidth=5;
markerSize=10;
lightGray=[0.66,0.66,0.66];
darkGray=[0.33,0.33,0.3];

% b=plot(,);
x=[0:length(long_short)-1]';
y=[long_short]';
figure('visible','off');
[hl,hp] = boundedline(x,y,b1,'-k');
axis square

set(hl,'LineStyle','-','color', 'k','LineWidth',lineWidth,'MarkerSize',markerSize)
xlim([0 60]);
hold on;
line([9 9],get(gca,'YLim'),'Color',darkGray)
% line([12 12],get(gca,'YLim'),'Color',darkGray)
line([18 18],get(gca,'YLim'),'Color',darkGray)
line([28 28],get(gca,'YLim'),'Color',darkGray)

title('Event-time quarterly returns','FontSize',titleSize);
ylabel('Cumulative Return (%)','FontSize',axisSize);
xlabel('Trading days relative to fiscal quarter start','FontSize',axisSize);
set(gca,'FontSize',axisSize');

set(gcf, 'PaperPositionMode', 'auto');
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
set(gca,'LooseInset',get(gca,'TightInset'))

textBoxText = {'1% of firms','announce','by 9th day'};
t = annotation('textbox', [0.28, 0.80, 0.1, 0.1], 'String', textBoxText);
t.FontSize = 16;
t.EdgeColor = [1 1 1];
annotation('arrow', [0.32, 0.345], [0.81, 0.75]);

textBoxText = {'25% of firms','announce','by 18th day'};
t = annotation('textbox', [0.345, 0.13, 0.1, 0.1], 'String', textBoxText);
t.FontSize = 16;
t.EdgeColor = [1 1 1];
annotation('arrow', [0.38, 0.41], [0.23, 0.3]);


textBoxText = {'75% of firms','announce','by 28th day'};
t = annotation('textbox', [0.49, 0.13, 0.1, 0.1], 'String', textBoxText);
t.FontSize = 16;
t.EdgeColor = [1 1 1];
annotation('arrow', [0.518, 0.488], [0.23, 0.3]);

export_fig('Figures/figure6.pdf');
% print -dpdf 'Figures/figure6.pdf' -bestfit

%% Figure 7 -trading costs

clear
clc

load ret
load me
load ff
load tcosts
load dates
load nyse
load tcosts
load varStruct

% Determine the start date & number of portfolios
s = find(dates==197501);
e = find(dates == 201712);
nPtf = 5;


% Find the FF49 results
r = find(strcmp([varStruct.label],{'FF49'}));
oilChanges = varStruct(r).wtiMonthlyBetas.quarterlyChanges;
oilResponseForecast = varStruct(r).predictedCAR.Changes;


% Run the time-series regressions
ind = makeUnivSortInd(oilResponseForecast, nPtf, NYSE);
res = runUnivSort(ret, ind, dates, me, 'plotFigure', 0, ...
                                       'printResults', 0, ...
                                       'factorModel', 1, ...
                                       'tcosts', tcosts);

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

orpNetRets = pret(s:e,end) - res.tcostsTS(s:e);

res = (nanols(orpNetRets,const(s:e)));

a = [{'Oil Surprise'},{res.beta},{sqrt(12)*nanmean(orpNetRets)/nanstd(orpNetRets)}];

% Read in the 23 anomalies
[anoms23, labels23] = getAnomalySignals('novyMarxVelikovAnomalies.csv', 1, 2);

labels=[{'Size'},{'Value'},{'Gross Profitability'},{'ValProf'},{'Accruals'},...
    {'Asset Growth'},{'Investment'},{'Piotroski F-score'}, {'Net Issuance (M)'},...    
    {'ROE'},{'Failure Probability'},{'ValMomProf'},{'ValMom'},...
    {'Idiosyncratic Volatility'},{'Momentum'},{'PEAD (SUE)'},{'PEAD (CAR3)'},...
    {'Industry Momentum'},...
    {'Industry Relative Reversals'},{'High-frequency Combo'},{'Short-run Reversals'},...
    {'Seasonality'},{'IRR (LowVol)'}];

% Store the number of anomalies
nAnoms = size(anoms23, 3);

for i = 1:nAnoms
    sortVar = anoms23(:, :, i);
    ind = makeUnivSortInd(sortVar, 5, NYSE);
    ind = 1 * (ind==1) + ...
          2 * (ind==5);
    res = runUnivSort(ret, ind, dates, me, 'tcosts', tcosts, ...
                                           'plotFigure', 0, ...
                                           'printResults', 0, ...
                                           'factorModel', 1); 
    res2 = nanols(res.netpret(s:end, end), const(s:end));
    
    a = [a; labels(i),{res2.beta},{sqrt(12)*nanmean(res.netpret(s:end))/nanstd(res.netpret(s:end))}];
                                       
end

temp=cellfun(@(x) x,a(:,3));
[B,I]=sort(temp,'descend');
head=a(I,1);
fHand = figure('visible','off');
aHand = axes('parent', fHand);
hold(aHand, 'on')
lightGray=[0.66,0.66,0.66];
darkGray=[0.33,0.33,0.3];

for i = 1:numel(B)
    if I(i)~=1
        bar(i, B(i), 'parent', aHand, 'facecolor',lightGray);
    else
        bar(i, B(i), 'parent', aHand, 'facecolor','k');
    end
end
axisSize=15;
titleSize=20;
title('Average net Sharpe ratios to different strategies','FontSize',titleSize);
ylabel('Net Sharpe ratio','FontSize',axisSize);
xlabel('Anomaly','FontSize',axisSize);
set(gca,'FontSize',axisSize);
set(gcf,'units','normalized','outerposition',[0 0 1 1])
set(gca,'LooseInset',get(gca,'TightInset'))
set(gcf, 'PaperPositionMode', 'auto');
set(gca,'XTick',1:24)
set(gca,'XTickLabel',head);
xtickangle(90);

export_fig('Figures/figure7.pdf','-transparent');