function commodityBetas=makeCommodityBetas(wCAR3,commoditySpotPrices,horizon,min_quarters)
% Usage: commodityBetas=makeCommodityBetas(commoditySpotPrices,horizon,min_quarters,labels,method)
%        Inputs:
%              - wCAR3 - wCAR3 matrix
%              - commoditySpotPrices - structure containing prices & changes/returns 
%              - horizon - horizon for beta estimatino (in months)
%              - min_quarters - minimum quarters for betas estimation (in quarters)

load dates

s = find(dates==197012);
commodityPriceChangeBetas  = nan(size(wCAR3));
commodityPriceReturnsBetas = nan(size(wCAR3));

for i = s+(min_quarters*3)-3:3:(length(dates)-3)

%         fprintf('Working on month %d\n',dates(i));
         % This gives us only the quarter-end rows in lhs, i.e., March, June, Sept, Dec. Assumption is if they are recorded in lhs for, say, Dec, they are last quarter's numbers. 
        y_matrix = wCAR3(i-horizon+3:3:i,:);
        
        x1 = commoditySpotPrices.quarterlyChanges(i-horizon:3:i-3); % Lagged oil price changes
        x2 = commoditySpotPrices.quarterlyReturns(i-horizon:3:i-3); % Lagged oil price returns

        const = ones(size(x1));
        horizontal_index = find( sum(isfinite(y_matrix),1)>=min_quarters & ...
                                 isfinite(y_matrix(end,:)) );

        for j = 1:length(horizontal_index)
            k = horizontal_index(j);        
            res1 = nanols(y_matrix(:,k), [const x1]);
            commodityPriceChangeBetas(i,k) = res1.beta(2);
            res2 = nanols(y_matrix(:,k),[const x2]);
            commodityPriceReturnsBetas(i,k) = res2.beta(2);
        end
end
        
commodityBetas.Changes = commodityPriceChangeBetas;
commodityBetas.Returns = commodityPriceReturnsBetas;
commodityBetas.quarterlyChanges = commoditySpotPrices.quarterlyChanges;
commodityBetas.quarterlyReturns = commoditySpotPrices.quarterlyReturns;

