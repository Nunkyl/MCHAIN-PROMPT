function backbonePDBStruct = pdbbackbone(chains)
%PDBBACKBONE Return a PDB structure only with backbone atoms
% PDBBACKBONE(PDBStruct) returns a PDB structure that contains only
% backbone atoms (that is, atoms, which names are N, C, and CA).
%
% See also pdbextractcoords
%
% MCHAIN-PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2014.

backbonePDBStruct = chains;
nChains = length(backbonePDBStruct);

if isfield(backbonePDBStruct{1}, 'Model')
    for i = 1:nChains
        backbonePDBStruct{i} = backbonePDBStruct{i}.Model.Atom;
    end
end

for i = 1:nChains
    backbonePDBStruct{i} = ...
        backbonePDBStruct{i}(ismember(...
        {backbonePDBStruct{i}.AtomName}, {'N' 'C' 'CA'}));
end


end

