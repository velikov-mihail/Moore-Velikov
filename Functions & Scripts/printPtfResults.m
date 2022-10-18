function printPtfResults(pret, s, printBetaFlag)

load ff
nPtf = size(pret, 2);

% Estimate the average returns, alphas, and loadings
for i = 1:(nPtf)
    % Determine the y variable
    if i==(nPtf)
        y = pret(s:end, i);
    else
        y = pret(s:end, i) - rf(s:end);
    end
    
    % Average returns
    x = const(s:end);
    tempRes = ols(y, x);
    result.xret(i,1)  = tempRes.beta(1);
    result.txret(i,1) = tempRes.tstat(1);

    % CAPM
    x = [const(s:end) mkt(s:end)];
    tempRes = ols(y, x);
    result.alpha1(i,1)  = tempRes.beta(1);
    result.talpha1(i,1) = tempRes.tstat(1);

    % FF3
    x = [ff3(s:end, :)];
    tempRes = ols(y, x);
    result.alpha3(i,1)  = tempRes.beta(1);
    result.talpha3(i,1) = tempRes.tstat(1);
    
    % FF4
    x = [ff4(s:end, :)];
    tempRes = ols(y, x);
    result.alpha4(i,1)  = tempRes.beta(1);
    result.talpha4(i,1) = tempRes.tstat(1);
    
    % FF5
    x = [ff5(s:end, :)];
    tempRes = ols(y, x);
    result.alpha5(i,1)  = tempRes.beta(1);
    result.talpha5(i,1) = tempRes.tstat(1);
    
    % Betas
    result.betas(i,:)  = tempRes.beta(2:end);
    result.tbetas(i,:) = tempRes.tstat(2:end);
end



heads = [{'$r^e$'}, ...
         {'$\alpha^{\text{CAPM}}$'}, ...
         {'$\alpha^{\text{FF3}}$'}, ...
         {'$\alpha^{\text{FF3+UMD}}$'}, ...
         {'$\alpha^{\text{FF5}}$'}];

a = [result.xret result.alpha1 result.alpha3 result.alpha4 result.alpha5]';
tA = [result.txret result.talpha1 result.talpha3 result.talpha4 result.talpha5]';

mat2Tex(a, tA, heads, 2);

fprintf('\n\n');

if nargin>2 & printBetaFlag==1
    heads = [{'$\beta_{\text{MKT}}$'}, ...
             {'$\beta_{\text{SMB}}$'}, ...
             {'$\beta_{\text{HML}}$'}, ...
             {'$\beta_{\text{RMW}}$'}, ...
             {'$\beta_{\text{CMA}}$'}];

    a  = result.betas';
    tA = result.tbetas';

    mat2Tex(a, tA, heads, 2);
end