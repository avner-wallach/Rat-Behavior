function [pvalue,D,dd]=boottest(x,y,N,func,mode)
% using bootstrap random permutations to show
% that difference in means of x and y is segnificant
% N- number of permutations
% mode- 'mean','median'
replacement=1;
if(nargin<4)
    func=@mean;
    mode='single';
elseif(nargin==4 & isstr(func))
    eval(['func=@',func,';']);
    mode='single';
end

%remove NaNs:
x=reshape(x(~isnan(x)),1,[]);
y=reshape(y(~isnan(y)),1,[]);

Nx=numel(x);
Ny=numel(y);
Mx=func(x);
My=func(y);
% eval(['Mx=',mode,'(x);']);
% eval(['My=',mode,'(y);']);
D=Mx-My;
V=[x y];
for i=1:N
    if(replacement)
        ind=randi(numel(V),numel(V),1);
    else
        ind=randperm(numel(V));
    end
    VV=V(ind);

    xx=VV(1:Nx);
    yy=VV((Nx+1):end);

    m_xx=func(xx);
    m_yy=func(yy);
%     eval(['m_xx=',mode,'(xx);']);
%     eval(['m_yy=',mode,'(yy);']);
    dd(i)=m_xx-m_yy;
end

% edges=linspace(min(dd),max(dd),50);
% h=histc(dd,edges);
% bar(edges,h);
% hold on;
% vline(D);
if(strcmp(mode,'single'))
    if(D>0)
        pvalue=sum(dd>=D)/N;
    else
        pvalue=sum(dd<=D)/N;
    end
else %double-ended
    pvalue=sum(abs(dd)>=abs(D))/N;
end
end