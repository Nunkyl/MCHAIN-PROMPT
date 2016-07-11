function coords = pdbextractcoords(PDBStruct)
%PDBEXTRACTCOORDS Extract Cartesian atom coordinates from a PDB structure
%   PDBEXTRACTCOORDS(PDBStruct) extracts coordinates of the atoms that
%   constitute a protein from a PDB structure. The function returns a cell
%   array of matrices. Each matrix corresponds to a separate model from the
%   PDB structure.
%
%   See also restorecoords trmrestorecoords
%
% MCHAIN-PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2014.

nModels = length(PDBStruct.Model);
coords = cell(1, nModels);

for iModel = 1:nModels
    coords{iModel} = [[PDBStruct.Model(iModel).Atom.X]' ...
        [PDBStruct.Model(iModel).Atom.Y]' ...
        [PDBStruct.Model(iModel).Atom.Z]'];
end

end

