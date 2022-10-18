function results = prt_char(char,ind,mv) ; 

%time-series average of the x-sectional mean (F-MacB)

pnum = max(max(ind));

indtemp=ind;
V=zeros(1,pnum);

if nargin==2

for k=1:pnum 
    indtemp(ind~=k)=nan;
    indtemp(ind==k)=1;
    prtchar=char.*indtemp;
    prtchar=winsorize(prtchar,1);
    V(1,k)=tsx_mean(prtchar);
end;


elseif nargin==3

for k=1:pnum
    indtemp(ind~=k)=nan;
    indtemp(ind==k)=1;
    prtmv=mv.*indtemp;
    V(1,k)=tsx_mean(char,prtmv);
end;

end;

results = V;