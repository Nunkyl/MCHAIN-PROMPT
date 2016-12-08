function [ pdbStructure, firstChains, lastChains ] = getdata( name )
% GETDATA reads the initial pdb file 
%   Returns the initial pdb structure and two cell arrays of chains for the
%   first and last conformations (each element of the cell array contains 
%   data for one chain).
%
% MCHAIN-PROMPT toolbox for MATLAB

% By Elizabeth Lutsenko, 2016.

toolboxPath = fileparts(which('demo'));
pdbStructure = pdbread(fullfile(toolboxPath, name));

firstConf = pdbStructure; firstConf.Model = firstConf.Model(1);
lastConf = pdbStructure; lastConf.Model = lastConf.Model(end);

lastLetterCode = double(firstConf.Model.Atom(end).chainID);
nChains = lastLetterCode - 65 + 1;

if nChains <= 0 
    nChains = 1;
end

firstChains = cell(1, nChains);
lastChains = cell(1, nChains);
letter = 'A';

if nChains == 1
    firstChains{1} = firstConf.Model.Atom;
    lastChains{1} = lastConf.Model.Atom;
else
    for i=1:nChains
        ind = [firstConf.Model.Atom.chainID] == letter;
        firstChains{i} = firstConf.Model.Atom(ind);
        lastChains{i} = lastConf.Model.Atom(ind);
        letter = char(double(letter) + 1);
    end
end

end

