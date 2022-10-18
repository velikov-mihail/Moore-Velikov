function ind = cond_index_from_inds(ind1,x2,n2,n)
% PURPOSE: assigns indexes on conditional double sorted portfolios
%---------------------------------------------------
% USAGE: ind=cond_index_from_inds(ind1,x2,n2)      
%        Function that spits out a matrix of which combined bin each firm falls
%        under as well as a matrix that shows how the individual sorts map into
%        the double one. Note that B/M is available only in June of each year!
%---------------------------------------------------
% Inputs
%        -ind1 - index of assigned firms to bins 
%        -x2 - variable according to which to assign within first
%        -n2 - number of bins in the second sort
% Output
%        -ind - matrix with assigned firms to n1*n2 bins


n1 = max(max(ind1));                           % Number of bins in first sort
ind  = zeros(size(ind1));                      % Output matrix

% for i = 1:size(ind,1)
%     if (isfinite(ind1(i,:)))
%     for j = 1:n1                             % For each common bin
%             firstInd=(ind1(i,:)==j);
%             if nargin==4
%                 tempInd2=aprts(x2(i,firstInd),n2,NYSE(i,firstInd));
%             else
%                 tempInd2=aprts(x2(i,firstInd),n2);
%             end
%             ind(i,firstInd)=tempInd2+n2*(j-1);
%     end
%     end
% end


for i=1:n1
    temp=x2;
    temp(ind1~=i)=nan;
    if nargin==4
        if max(max(n)) == 1
            load NYSE
            tempIndex=makeUnivSortInd(temp,n2,NYSE);
        else 
            load me
            tempIndex=makeUnivSortInd(temp,n2,me);
        end
            
    else
        tempIndex=makeUnivSortInd(temp,n2);
    end
    ind(ind1==i & tempIndex>0)=tempIndex(ind1==i & tempIndex>0)+n2*(i-1);    
end


