function [pval,H]=my_boxplot(cell_data,cell_names,M,func,mode)
if(nargin<5)
    mode='single'
end
if(nargin<4)
    func=@mean;
end
if(nargin<3)
    M=1e3;
end
N=numel(cell_data);
%M=5e3;
M=1e3;
X=[];
G=cell(0);
for i=1:N
    x=numel(cell_data{i});
    X=[X cell_data{i}];
    Gx=cell(1,x);
    Gx(:)={cell_names{i}};
    G=[G Gx];
end
pval=nan(N-1,N);
for i=1:N
    for j=(i+1):N
        pval(i,j)=boottest(cell_data{i},cell_data{j},M,func,mode);
    end
end

H=boxplot(X,G);
