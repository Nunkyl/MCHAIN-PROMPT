function trmodel = trmcreate(PDBStruct1, PDBStruct2, nConf)
%TRMCREATE Create a transformation model
%   TRMCREATE(PDBStruct1, PDBStruct2, nModels) creates a transformation
%   model from two PDB structures which contain a single model and
%   correspond to the same protein. Note that amino acid content of both 
%   structures must be the same. The parameter nConf specifies the number 
%   of intermediate configurations in the transformation model.
%
%   See also createmodel pdb2trm
%
% MCHAIN-PROMPT Toolbox for MATLAB

% By Gaik Tamazian (2014), Elizabeth Lutsenko (2016).

if ~(iscell(PDBStruct1) && iscell(PDBStruct2))
    PDBStruct1 = {PDBStruct1};
    PDBStruct2 = {PDBStruct2};
end

% Remove non-backbone atoms from the specified PDB structures
PDBStruct1 = pdbbackbone(PDBStruct1);
PDBStruct2 = pdbbackbone(PDBStruct2);

nChains = length(PDBStruct1);

% Get atomic masses of backbone atoms.
trmodel = struct('m',           atomicmass(getSymbols(PDBStruct1)), ...
                 'StartCoords', atomiccoords(PDBStruct1));
 
% Get details about the bonds between the chains.
bondInd = virtualbonds2(atomiccoords(PDBStruct1), atomiccoords(PDBStruct2));
nBonds = size(bondInd,1);             
trmodel(1).bondInd = bondInd;

for i=1:nChains
    % Process side chains atoms - get a vector of their masses and add 
    % them to atomic masses of alpha carbons. Also add a mass of one 
    % hydrogen atom.
    
    alphaCarbonAtoms = PDBStruct1{i}( ...
        ismember({PDBStruct1{i}.AtomName}, {'CA'}));

    trmodel(i).m(2:3:end) = trmodel(i).m(2:3:end) + ...
        sidechainmass({alphaCarbonAtoms.resName}) + ...
        cell2mat(atomicmass({{'H'}}));

    % Add masses of hydrogen atoms to nitrogen atoms of the backbone.
    trmodel(i).m(1:3:end) = trmodel(i).m(1:3:end) + ...
        cell2mat(atomicmass({{'H'}}));

    % Add masses of oxygen atoms to carbon atoms of the backbone.
    trmodel(i).m(3:3:end) = trmodel(i).m(3:3:end) + ...
        cell2mat(atomicmass({{'O'}}));

    % Add a mass of a hydrogen atom to N-end of the protein.
    trmodel(i).m(1) = trmodel(i).m(1) + cell2mat(atomicmass({{'H'}}));

    % Add a mass of a hydroxil group to C-end of the protein.
    trmodel(i).m(end) = trmodel(i).m(end) + ...
        sum(cell2mat(atomicmass({{'O', 'H'}})));
end


for i=1:nChains
    trmodel(i).r     = zeros(length(trmodel(i).m) - 1, nConf+2);
    trmodel(i).alpha = zeros(length(trmodel(i).m) - 2, nConf+2);
    trmodel(i).psi   = zeros(length(trmodel(i).m) - 3, nConf+2);
end

for i=1:nBonds
    trmodel(i).bond_r     = zeros(1, nConf+2);
    trmodel(i).bond_alpha = zeros(2, nConf+2);
    trmodel(i).bond_psi   = zeros(3, nConf+2);
end

% Calculate bond lengths, planar and torsion angles of the first and last
% configurations.
for i=1:nChains
    [trmodel(i).r(:,1), trmodel(i).alpha(:,1), trmodel(i).psi(:,1)] = ...
        createmodel(PDBStruct1, i);
    [trmodel(i).r(:,end), trmodel(i).alpha(:,end), trmodel(i).psi(:,end)] = ...
        createmodel(PDBStruct2, i);
end

for i=1:nBonds
    [trmodel(i).bond_r(:,1), trmodel(i).bond_alpha(:,1), trmodel(i).bond_psi(:,1)] = ...
        createmodel(PDBStruct1, i, bondInd);
    [trmodel(i).bond_r(:,end), trmodel(i).bond_alpha(:,end), trmodel(i).bond_psi(:,end)] = ...
        createmodel(PDBStruct2, i, bondInd);
end

for i=1:nChains
    trmodel(i).r(:,2:end-1) = ...
        interpolate(trmodel(i).r(:,1), trmodel(i).r(:,end), nConf);
    trmodel(i).alpha = ...
        circinterp(trmodel(i).alpha(:,1), trmodel(i).alpha(:,end), nConf);
    trmodel(i).psi = ...
        circinterp(trmodel(i).psi(:,1), trmodel(i).psi(:,end), nConf);
end;

for i=1:nBonds
    trmodel(i).bond_r(:,2:end-1) = ...
        interpolate(trmodel(i).bond_r(:,1), trmodel(i).bond_r(:,end), nConf);
    trmodel(i).bond_alpha = ...
        circinterp(trmodel(i).bond_alpha(:,1), trmodel(i).bond_alpha(:,end), nConf);
    trmodel(i).bond_psi = ...
        circinterp(trmodel(i).bond_psi(:,1), trmodel(i).bond_psi(:,end), nConf);
end;

% Calculate the rotation matrix for superposition of the second
% conformation to the first one.
trmodel(1).U = cell(nConf+2, 1);
trmodel(1).U{1} = eye(3);
prevCoords = vertcat(trmodel(:).StartCoords);
for j = 2:nConf+2
    currCoords = restorecoords(trmodel, j);
    [~,prevCoords,t] = procrustes(prevCoords, currCoords, ...
        'scaling', false, 'reflection', false);
    trmodel(1).U{j} = t.T;
end

end

function result = interpolate(start, finish, M)
    result = ...
        repmat((M:-1:1)/(M+1), length(start), 1).*repmat(start, 1, M) + ...
        repmat((1:M)/(M+1), length(finish), 1).*repmat(finish, 1, M);
end

function coords = atomiccoords(PDBStruct)
    coords =  cell(1, length(PDBStruct));
    for i=1:length(PDBStruct)
        coords{i} = [[PDBStruct{i}.X]' [PDBStruct{i}.Y]' [PDBStruct{i}.Z]']; %each one is a column
    end
end

function symbols = getSymbols(PDBStruct)
    symbols = cell(1, length(PDBStruct));
    for i=1:length(PDBStruct)
        symbols{i} = {PDBStruct{i}.element};
    end
end

