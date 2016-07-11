function h = pdbplotfixedrmsd(pdbStruct, modelNo, name)
%PDBPLOTFIXEDRMSD Plot RMSDs between models and the specified one
%   PDBPLOTFIXEDRMSD(pdbStruct, modelNo) plots RMSDs between all 
%   models of transformations specified in the cell array pdbStruct 
%   and the model with the specified number modelNo.
%
%   See also pdbplotadjrmsd trmplotadjrmsd trmplotfixedrmsd
%
% MCHAIN-PROMPT Toolbox for MATLAB

% By Gaik Tamazian (2014), Elizabeth Lutsenko (2016).

if ~iscell(pdbStruct)
    pdbStruct = {pdbStruct};
end

nTrans = length(pdbStruct);
nModels = length(pdbStruct{1}.Model);
rmsdValues = zeros(nTrans, nModels);

for i = 1:nTrans
    coords = pdbextractcoords(pdbStruct{i});
    for j = 1:nModels
        % superpose the current model to the specified one
        [~, coords{j}] = procrustes(coords{modelNo}, coords{j}, ...
            'scaling', false, 'reflection', false);
        rmsdValues(i,j) = mean(sqrt(sum((coords{j} - ...
            coords{modelNo}).^2,2)));
    end
end

h = plot(transpose(rmsdValues),'-o');
xlabel('Configuration Number');
ylabel('RMSD in A');

% plot will be saved if a name is provided
if exist('name', 'var')
    name = strcat(name, '_fixed');
    print(name, '-dpng');
end

end

