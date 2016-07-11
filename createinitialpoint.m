function [ initial_point ] = createinitialpoint(trmodel,  angleIndices)
% CREATEINITIALPOINT (trmodel,  angleIndices) creates the initial 
%   approximation for fmincon. Takes all the data about the virtual bonds 
%   (bond lengths, planar and torsion angles)for the intermediate 
%   conformations and all the torsion angles specified in angleIndices and 
%   combines them into one vector.
%
% MCHAIN-PROMPT Toolbox for MATLAB

% By Elizabeth Lutsenko, 2016.


if ~iscell(trmodel) && ~isstruct(trmodel)
    trmodel = {trmodel};
end

if ~iscell(angleIndices)
    angleIndices = {angleIndices};
end

nChains = length(trmodel);
nConfs = size(trmodel(1).r, 2) - 2;

if nChains == 1
    initial_point = (trmodel(1).psi(angleIndices{1},2:end-1));
    initial_point = reduceangles(initial_point);
    initial_point = initial_point(:);
    return
end

len = 0;
for i=1:nChains-1
    len = len + length(angleIndices{i}) + size(trmodel(i).bond_r,1) + ...
        size(trmodel(i).bond_alpha,1) + size(trmodel(i).bond_psi,1);
end
len = len + length(angleIndices{nChains});

chainData = zeros(len, nConfs);

curInd = 1;
% include virtual r
for i=1:nChains-1
    chainData(curInd,:) = trmodel(i).bond_r(2:end-1);
    curInd = curInd + size(trmodel(i).bond_r,1);
end

% include virtual alpha
for i=1:nChains-1
    n = curInd + size(trmodel(i).bond_alpha,1) - 1;
    chainData(curInd:n,:) = trmodel(i).bond_alpha(:,2:end-1);
    curInd = n + 1;
end

% include virtual psi
for i=1:nChains-1
    n = curInd + size(trmodel(i).bond_psi,1) - 1;
    chainData(curInd:n,:) = trmodel(i).bond_psi(:,2:end-1);
    curInd = n + 1;
end

% include psi
for i=1:nChains
    n = curInd + length(angleIndices{i}) - 1;
    chainData(curInd:n,:) = trmodel(i).psi(angleIndices{i},2:end-1);
    curInd = n + 1;
end

% reduce the angles
n = size(trmodel(1).bond_r,1) * (nChains-1) + 1;
chainData(n:end,:) = reduceangles(chainData(n:end,:));

% stretch everything into one column
chainData = chainData';
initial_point = chainData(:);

end



