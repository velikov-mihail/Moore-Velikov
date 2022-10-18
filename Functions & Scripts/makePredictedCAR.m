function [predictedCAR] = makePredictedCAR(commodityBetas)

predictedCAR = struct;
nMonths = size(commodityBetas(1).Changes, 1);
nStocks = size(commodityBetas(1).Changes, 2);
nBetas = length(commodityBetas);

if nBetas==1    
    rptdChanges = repmat(commodityBetas.quarterlyChanges, 1, nStocks);
    predictedCAR.Changes = commodityBetas.Changes .* rptdChanges;

    rptdReturns = repmat(commodityBetas.quarterlyReturns, 1, nStocks);
    predictedCAR.Returns = commodityBetas.Returns .* rptdReturns;
else        
    predictedCAR.Changes    = zeros(nMonths, nStocks);
    predictedCAR.Returns    = zeros(nMonths, nStocks);
    predictedCAR.AvgChanges = zeros(nMonths, nStocks);
    predictedCAR.AvgReturns = zeros(nMonths, nStocks);

    for i=1:nBetas
        
        rptdChanges = repmat(commodityBetas(i).spotPrices.quarterlyChanges, 1, nMonths);
        predictedCAR.Changes=predictedCAR.Changes + commodityBetas(i).Changes.*rptdChanges;

        rptdReturns = repmat(commodityBetas(i).spotPrices.quarterlyReturns, 1, nMonths);
        predictedCAR.Returns=predictedCAR.Returns + commodityBetas(i).Returns .* rptdReturns;
        
        rptdAvgChanges = repmat(commodityBetas(i).spotPrices.AvgQuarterlyChanges, 1, nMonths);
        predictedCAR.AvgChanges=predictedCAR.AvgChanges + commodityBetas(i).AvgChanges .* rptdAvgChanges;
        
        rptdAvgReturns = repmat(commodityBetas(i).spotPrices.AvgQuarterlyReturns, 1, nMonths);
        predictedCAR.AvgReturns=predictedCAR.AvgReturns + commodityBetas(i).AvgReturns .* rptdAvgReturns;
    end
end
