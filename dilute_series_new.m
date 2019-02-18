function [X_out,XX_out,lag]=dilute_series_new(X,N,ind,Nmax)
if nargin<4
    Nmax=30; %max segment length (split in half if above)
end
if nargin<3
    ind=1:numel(X);
end
K=50;
if(numel(ind)<2)
    X_out=[];
    XX_out=[];
    lag=[];
    return;
end

if(numel(X)>Nmax)
    Y{1}=X(1:floor(end/2));
    Ind{1}=ind(ind<=numel(Y{1}));
    Y{2}=X((floor(end/2)+1):end);
    Ind{2}=ind(ind>numel(Y{1}))-numel(Y{1});
else
    Y{1}=X;
    Ind{1}=ind;
end
XX_out=[];
X_out=[];
lag=[];
for i=1:numel(Y)
    rec_ind=setdiff(1:numel(Y{i}),Ind{i});
    Y{i}(rec_ind)=mean(Y{i}(Ind{i})); %replace with mean=>neutral for correlations
    XX=xcov(Y{i}(:),K,'coeff')';
    XX=abs(XX(ceil(end/2):end)); %take only half
    sig_val=get_sigval(Y{i},N,K);
    maxlag=min(K,numel(XX));
    index=find(XX(1:maxlag)<sig_val(1:maxlag));
    if(numel(index)==0)
        continue;
    end
%     try
        lag=[lag;index(1)-1];
%     catch
%         figure; plot(XX); hold on; plot(sig_val);
%     end
    Y{i}(rec_ind)=NaN;
    Y{i}=smooth(Y{i},lag(end));
    ind_dilute=[ceil(lag(end)/2):lag(end):numel(Y{i})];
    ind_dilute=intersect(ind_dilute,Ind{i});
    X_out=[X_out;Y{i}(ind_dilute)];
    XX_out=[XX_out;XX]; %take only half 

%     X_out=X_out(:);
end
end

function s=get_sigval(X,M,K)

for i=1:M
    x1=X(randperm(numel(X)));
    tmp=xcov(x1,K,'coeff')';
    xx(i,:)=abs(tmp(ceil(end/2):end));    
end
s=quantile(xx,0.95);
end