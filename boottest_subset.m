function pvalue=boottest_subset(x,y,N,func,mode)
% using bootstrap random permutations to show
% that difference in means of x and y is segnificant. x is a subset of y
% N- number of permutations
% mode- 'single','double'
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
D=Mx-My;
V=[y];
for i=1:N
    if(replacement)
        ind=randi(numel(V),numel(V),1);
    else
        ind=randperm(numel(V));
    end
    VV=V(ind);

    xx=VV(1:Nx);

    m_xx=func(xx);
    dd(i)=m_xx-My;
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