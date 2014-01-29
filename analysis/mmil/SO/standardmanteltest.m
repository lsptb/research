function [p,r0,r]=standardmanteltest(mata, matb, nperm)
% source: http://www-stat.stanford.edu/~susan/phylo/index/node64.html
% function [p,r0,r]=manteltest(mata, matb, nperm)
% modified last by JSS on 26-Oct-2010
%   - added normalization
        veca=asvect(mata); stda = std(veca); veca = veca - mean(veca);
        vecb=asvect(matb); stdb = std(vecb); vecb = vecb - mean(vecb);
        r = zeros(nperm, 1);
        [m,ncol]=size(mata);
        
        n  = ((m-1)*m)/2;
        r0 = (veca/stda)*(vecb/stdb)' / (n-1);
        
        for(i = (1:nperm)) 
           ni=randperm(n);                     
           r(i) =(veca/stda)*(vecb(ni)/stdb)'/(n-1); 
        end

        p = (sum(r >= r0))/(nperm + 1);



function vect=asvect(mat)
       [m,ncol]=size(mat);
        n=((m-1)*m)/2;
       k=1;
       vect=zeros(1,n);
       for (i = (2:m))
           for ( j = (1:(i-1)))
               vect(k)=mat(i,j);
               k=k+1;
           end
       end