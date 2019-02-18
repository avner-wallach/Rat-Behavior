function [X_out,XX_out,lag]=dilute_series(X,N,ind)
if nargin<3
    ind=1:numel(X);
end
K=50;
if(numel(ind)<2)
    X_out=X(ind);
    XX_out=zeros(1,K+1);
    lag=0;
    return;
end
rec_ind=setdiff(1:numel(X),ind);
X(rec_ind)=mean(X(ind)); %replace with mean=>neutral for correlations
XX=xcov(X(:),K,'coeff')';
% [m,I]=max(XX);
XX_out=abs(XX(ceil(end/2):end)); %take only half (no ABS)
XX=abs(XX(ceil(end/2):end)); %take only half
sig_val=get_sigval(X,N,K);
maxlag=min(K,numel(XX));
index=find(XX(1:maxlag)<sig_val(1:maxlag));
try
    lag=index(1)-1;
catch
    figure; plot(XX); hold on; plot(sig_val);
end
X(rec_ind)=NaN;
X=smooth(X,lag);
ind_dilute=[ceil(lag/2):lag:numel(X)];
ind_dilute=intersect(ind_dilute,ind);
X_out=X(ind_dilute);
%X_out=decimate(X,lag,lag,'fir'); %lowpass & downsample
X_out=X_out(:);
end

function s=get_sigval(X,M,K)

for i=1:M
    x1=X(randperm(numel(X)));
    tmp=xcov(x1,K,'coeff')';
    xx(i,:)=abs(tmp(ceil(end/2):end));    
end
s=quantile(xx,0.95);
end