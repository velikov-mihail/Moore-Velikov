
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

indSharpe = nan(49, 3);
indMean = nan(49, 1);
indStd = nan(49, 1);
indSkew = nan(49, 1);
indKurt = nan(49, 1);

for i = 1:49
    temp = oilChanges;
    temp(FF49~=i) = nan;
    indMean(i) = mean(mean(temp(s:end,:), 2, 'omitnan'), 'omitnan');
    indStd(i) = mean(std(temp(s:end,:), [], 2, 'omitnan'), 'omitnan');
    indKurt(i) = mean(kurtosis(temp(s:end,:), [], 2), 'omitnan');
    indSkew(i) = mean(skewness(temp(s:end,:), [], 2), 'omitnan');

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
        res(i).sharpe(j) = sqrt(12)*mean(pret(:,j), 'omitnan')/std(pret(:,j), 'omitnan');
    end    
    indSharpe(i,:) = res(i).sharpe;
end


a = [100*indMean 100*indStd indSkew indKurt indSharpe(:, 3)];
h= FF49Names.longName;

[h, ID] = sort(h);
h = regexprep(h, '&', '\\&');
a = a(ID, :);

mat2Tex(a, a, h, 2);