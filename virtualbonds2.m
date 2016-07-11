function [ bondList ] = virtualbonds2(conf1, conf2)
%VIRTUALBONDS2 Returns the virtual bonds between chains
%   VIRTUALBONDS2(conf1, conf2) finds the most stable atoms in each chain. 
%   Based on these atoms bondList is formed. One row represents one bond.
%
% MCHAIN-PROMPT Toolbox for MATLAB

% By Elizabeth Lutsenko, 2016.

nChains = length(conf1);

if nChains == 1
    bondList = [];
    return
end

nChains = length(conf1);
bondList = ones(nChains-1,4);
bondList(:,2) = 2:nChains;
newPoints = newbegginingpoints(conf1, conf2);
bondList(:,3) = newPoints(1)*ones(nChains-1,1);
bondList(:,4) = newPoints(2:end);

end

