function results = tsx_mean(char,mv) ; 

%time-series average of the x-sectional mean (F-MacB)

% equal weight in the x-section 
if nargin ==1    
    results = nanmean(nanmean(char'));
  
% value weight in the x-section
elseif nargin==2
    mvcoincident=mv;
    mvcoincident(isfinite(char)==0)=nan;
    results = nanmean(nanmean((char.*mvcoincident)')./nanmean(mvcoincident'));
end ; 
