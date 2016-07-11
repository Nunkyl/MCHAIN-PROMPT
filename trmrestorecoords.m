function coords = trmrestorecoords(trmodel)
%TRMRESTORECOORDS Restore Cartesian coordinates of transformation atoms
%   TRMRESTORECOORDS(trmodel) returns a cell array of matrices containing
%   Cartesian coordinates of the atoms that constitute configurations of
%   the transformation trmodel.
%
%   See also restorecoords
%
% MCHAIN-PROMPT Toolbox for MATLAB

% By Gaik Tamazian (2014), Elizabeth Lutsenko (2016).

if ~iscell(trmodel) && ~isstruct(trmodel)
    trmodel = {trmodel};
end

nConf = size(trmodel(1).psi, 2);
nChains = length(trmodel);
nElements = zeros(1, nChains);
for i=1:nChains
    nElements(i) = size(trmodel(i).m, 1);
end
nAtoms = sum(nElements);

coords = cell(1, nConf);
coords{1} = vertcat(trmodel(:).StartCoords);
firstConfTranslation = repmat(mean(coords{1}, 1), nAtoms, 1);

for i = 2:nConf
    coords{i} = restorecoords(trmodel, i);
    coords{i} = coords{i}*trmodel(1).U{i};
    
    % apply the translation
    currTranslation = repmat(mean(coords{i}, 1), nAtoms, 1);
    coords{i} = coords{i} - currTranslation + firstConfTranslation;       
end

end

