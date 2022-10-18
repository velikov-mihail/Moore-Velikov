%% Timekeeping

fprintf('Now working on the daily betas. Run started @ %s.\n\n\n',char(datetime('now')));

%%

clear
clc

load dret
load dff
load ddates

b = find(isfinite(dmkt), 1, 'first');
hor  = 252;
const = ones(hor,1);
nDays = length(dmkt);

X = [dmkt lag(dmkt,1,nan)];

dcapmb = nan(size(dret));

for i = (b+hor-1):nDays
    
    dRange = (i-hor+1:i);
    thisX = X(dRange, :);
    thisY = dret(dRange,:) - repmat(drf(dRange), 1, size(dret, 2));
    index = find(isfinite(sum(thisY, 1)));
    for j = index
        res = ols(thisY(:, j),[const thisX]);
        dcapmb(i, j) = res.beta(2) + res.beta(3);
    end
end

save -v7.3 Data\dcapmb dcapmb