function [ trmodel ] = modelopt( trmodel, angleIndices, vars )
% MODELOPT inserts all the variables into their places in trmodel
%
% MCHAIN-PROMPT Toolbox for MATLAB

% By Elizabeth Lutsenko, 2016.

if ~iscell(trmodel) && ~isstruct(trmodel)
    trmodel = {trmodel};
end

if ~iscell(angleIndices)
    angleIndices = {angleIndices};
end

nConfs = size(trmodel(1).r, 2) - 2;
nChains = length(trmodel);

if nChains == 1
    trmodel(1).psi(angleIndices{1}, 2:end-1) = reshape(vars, ...
        length(angleIndices{1}), nConfs);
    return
end

curInd = 1;
% Restore virtual r
for i=1:nChains-1
    trmodel(i).bond_r(2:end-1) = vars(curInd:curInd+nConfs-1);
    curInd = curInd + nConfs;
end

% Restore virtual alpha
for i=1:nChains-1
    trmodel(i).bond_alpha(:,2:end-1) = reshape(vars(curInd:curInd+...
        nConfs*2-1), nConfs, 2)'; 
    curInd = curInd + nConfs*2;
end

% Restore virtual psi
for i=1:nChains-1
    trmodel(i).bond_psi(:,2:end-1) = reshape(vars(curInd:curInd+...
        nConfs*3-1), nConfs, 3)';
    curInd = curInd + nConfs*3;
end

% Restore psi
for i=1:nChains
    n = length(angleIndices{i});
    trmodel(i).psi(angleIndices{i},2:end-1) = reshape(vars(curInd:...
        curInd+nConfs*n-1),...
         nConfs, n)';
    curInd = curInd + nConfs*n;
end


end

