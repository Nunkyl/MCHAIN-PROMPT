function [h, rmsdValues] = pdbplottrfcmp(pdbStruct1, pdbStruct2, name)
%PDBPLOTTRFCMP Plot RMSDs between configurations of the PDB structures.
%   PDBPLOTTRFCMP(pdbStruct1, pdbStruct2) plots RMSDs between
%   configurations of the specifed PDB structures representing a protein
%   transformation. It also returns the RMSD values as an optional output
%   argument.
%
%   See also pdbplotadjrmsd pdbplotfixedrmsd
%
% MCHAIN-PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2015.

nModels = length(pdbStruct1.Model);
rmsdValues = zeros(nModels, 1);
coords1 = pdbextractcoords(pdbStruct1);
coords2 = pdbextractcoords(pdbStruct2);
for i = 1:nModels
    % superpose points of the second model to points of the first one
    [~, coords2{i}] = procrustes(coords1{i}, coords2{i}, ...
        'scaling', false, 'reflection', false);
    rmsdValues(i) = mean(sqrt(sum((coords1{i} - coords2{i}).^2, 2)));
end

h = plot(transpose(rmsdValues), '-o');
xlabel('Configuration Number');
ylabel('RMSD in AA');

% If you fancy saving this plot provide a name for it.
if exist('name', 'var')
    name = strcat(name, '_cmp');
    print(name, '-dpng');
end
end

