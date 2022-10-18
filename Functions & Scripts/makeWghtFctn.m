function retMatrix = makeWghtFctn(inds, xret_tab)

load('ret.mat', 'ret');
load('RDQ.mat', 'RDQ');
load('dates.mat', 'dates');
load('FQTR.mat', 'FQTR');
load('permno.mat', 'permno');


% Store a few constants
nMonths = size(ret, 1);
nStocks = size(ret, 2);

% Pick the start date
s = find(dates==197003);

% Initialize the return matrix
retMatrix = nan(nMonths, nStocks);

% Make a vector with the month numbers
mths = dates - 100*floor(dates/100);

parfor i = s:nMonths
    % Need this to use parfar (looping over the calendar quarter-ends)
    if ismember(mths(i), [3 6 9 12])
        
        % Store the last three months' RDQs, industries, and quarter-ends
        tempRDQ = RDQ(i-2:i,:);
        tempInds = inds(i-2:i,:);
        tempDataDateQ = FQTR(i-2:i,:);
        tempPermno = repmat(permno', 3, 1);
        
        % Find the observations with fiscal quarter-ends that match the end
        % of last calendar quarter (i.e., dates(i-3())
        ind = isfinite(tempDataDateQ) & ...
              isfinite(tempRDQ)       & ...
              (tempInds)>0            & ...
              floor(tempDataDateQ/100) == dates(i-3);
          
        % Save the end of last quarter
        lastQuarterEnd = datetime(mode(tempDataDateQ(ind)),'ConvertFrom','YYYYMMDD');
        
        % Create a table with the indexed observations
        tempData = array2table([tempRDQ(ind) tempInds(ind) tempPermno(ind)], ...
                               'VariableNames', {'ddates','inds','permno'});
        tempData.daysSinceQtrStart = days(datetime(tempData.ddates,'ConvertFrom','YYYYMMDD') - lastQuarterEnd);
        tempData.rawWeight = 1./(tempData.daysSinceQtrStart.^2);
        
        
        tempRes = nan(1, length(permno));
        for j = 1:height(tempData)        
            r = find(permno==tempData.permno(j));        
            indData = tempData(tempData.inds==tempData.inds(j), :);
            indData.permno = repmat(tempData.permno(j), height(indData), 1);
            temp_xret_teb = xret_tab(xret_tab.permno==tempData.permno(j), :);
            indData = outerjoin(indData, temp_xret_teb, 'Type', 'Left', 'MergeKeys', 1);
            indData(isnan(indData.xret),:) = [];
            tempRes(r) = sum(indData.rawWeight .* indData.xret) / sum(indData.rawWeight);
        end 
        retMatrix(i, :) = tempRes;
    end
end
