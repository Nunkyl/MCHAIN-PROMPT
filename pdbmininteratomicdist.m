function minDist = pdbmininteratomicdist(PDBStruct)
%PDBMININTERATOMICDIST Returns minimal interatomic distances
%   PDBMININTERATOMICDIST(PDBStruct) returns the vector of the minimal
%   interatomic distances between atoms within transformation 
%   configurations. The transformation is specified by the PDB structure
%   PDB struct. Each vector element corresponds to the certain PDB model.
%
%   See also trmmininteratomicdist
%
% MCHAIN-PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2014.

coords = pdbextractcoords(PDBStruct);
nModels = length(coords);
minDist = nan(nModels, 1);

for iModel = 1:nModels
	minDist(iModel) = min(pdist(coords{iModel}));
end

end

