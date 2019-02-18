function [ c,p ] = nancorr( x,y,N)
%NANCORR correlation coefficient ignoring nans
%   c- Pearson's correlation coefficient  
%   p- p value (random permuations)
%   N- number of interations
if(nargin<3)
    N=1e3;
end
ind=find(~isnan(x) & ~isnan(y));
x=x(ind); y=y(ind);
c=corr(x,y);
for n=1:N
    ctl(n)=corr(x,y(randperm(numel(y))));
end
p=sum(abs(ctl)>=abs(c))/N;

end

