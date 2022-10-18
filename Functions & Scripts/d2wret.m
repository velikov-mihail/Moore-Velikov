function res = d2wret(dmat,eowflag)

% daily to weekly (returns)

tmat = dmat;
index = isnan(dmat);
tmat(index) = 0;
tmat = cumprod(1+tmat);
tmat(index) = nan;
tmat = tmat(eowflag,:);
res = tmat./lag(tmat,1,nan)-1;


end

