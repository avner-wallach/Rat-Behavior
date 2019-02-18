function [lturn_ind,rturn_ind,amps,ints]=find_head_turns(beta,minamp,maxint);
% find indices of head turns 
% alpha- head angle in deg
% minamp- minimum amplitude of turn; 
% maxint-maximal interval of turn

% for i=1:numel(hdata)
%     alpha=hdata{i}(:,1);
    [pks,plcs]=findpeaks(beta);
    [trghs,tlcs]=findpeaks(-beta);
    A=sortrows([plcs pks ; tlcs -trghs]);
    amps=diff(A(:,2));
    ints=diff(A(:,1));
    inds=find(amps>=minamp & ints<=maxint);
    lturn_ind = A(inds,1);
    inds=find(amps<=-minamp & ints<=maxint);
    rturn_ind = A(inds,1);
% end
end
