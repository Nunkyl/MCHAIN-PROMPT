function [ min_ind ] = newbegginingpoints(conf1, conf2)
%NEWBEGGININGPOINTS Finds atoms for creating virtual bonds
%   NEWBEGGININGPOINTS (conf1, conf2) counts the sum of distances from each 
%   atom to all the other atoms for the first and last conformations. 
%   Returns the IDs of the most stable atoms in the protein (one in each 
%   chain) 
%
% MCHAIN-PROMPT Toolbox for MATLAB

% By Elizabeth Lutsenko, 2016.

nChains = length(conf1);
chainLength = zeros (nChains);
for i=1:nChains
    chainLength(i) = length(conf1{i});
end

dist1 = distance(conf1);
dist2 = distance(conf2);

dist = abs(dist1 - dist2);

min_ind = zeros(nChains,1);

 for i=1:nChains
     d = dist(sum(chainLength(1:i-1))+1:sum(chainLength(1:i-1)) + ...
                            chainLength(i)); 
     [~,min_ind(i)] = min(d);
 end

end

function [dist]  = distance (conf)

new_conf = vertcat(conf{:});

N = length(new_conf);

dist = zeros(N, N);

for i=1:N
    for j=i+1:N
        dist(i,j) = sqrt(sum((new_conf(i,:) - new_conf(j,:)).^2));
        dist(j,i) = dist(i,j);
    end
end

dist = sum(dist);

end