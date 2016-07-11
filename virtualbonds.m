function [ bondList ] = virtualbonds(conf1, conf2)
%VIRTUALBONDS (conf1, conf2) gets coordinates of new virtual bonds
%   VIRTUALBONDS(conf1, conf2) finds the most stable bonds between chains. 
%
% MCHAIN-PROMPT Toolbox for MATLAB

% By Elizabeth Lutsenko, 2016.
  

nChains = length(conf1);

if nChains == 1
    bondList = [];
    return
end

weights = zeros(nChains, nChains);
ind = zeros(nChains, nChains, 2);

for i=1:nChains
    for j=i+1:nChains
        [weights(i,j), ind(i,j,:)] = mindist(conf1{i}, conf1{j}, conf2{i}, conf2{j});
        weights(j,i) = weights(i,j);
        ind(j,i,1) = ind(i,j,2);
        ind(j,i,2) = ind(i,j,1);
    end
end

weights = sparse(weights);

[spantree,~] = graphminspantree(weights);
spantree = full(spantree);
spantree = tril(spantree)+tril(spantree,-1)';

bondList = zeros(nChains-1,4);

bondList(:,1:2) = dfes(1, [], spantree);

for i=1:nChains-1
    bondList(i,3:4) = ind(bondList(i,1),bondList(i,2),:);
end

end

function [weight, ind]  = mindist (ch1_beg, ch2_beg, ch1_end, ch2_end) %ch = chain

ch1_len = length(ch1_beg);
ch2_len = length(ch2_beg);

weight = 1000000;
ind = [1,1];

for i=1:ch1_len
    for j=1:ch2_len
        dist_beg = sqrt(sum((ch1_beg(i,:) - ch2_beg(j,:)).^2));
        
        dist_end = sqrt(sum((ch1_end(i,:) - ch2_end(j,:)).^2));
        
        if abs(dist_beg - dist_end) < weight
            weight = abs(dist_beg - dist_end);
            ind = [i,j];
        end
    end   
end

if weight == 0, weight = realmin; end
end
