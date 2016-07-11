function [h, rmsdValues] = pdbplotadjrmsd(pdbStruct, name)
%PDBPLOTADJRMSD Plot RMSDs between adjacent models
%   PDBPLOTADJRMSD(pdbStruct, fitModels) plots RMSDs between adjacent 
%   models of PDB structures specified in the cell array pdbStruct.
%
%   See also pdbplotfixedrmsd trmplotadjrmsd trmplotfixedrmsd
%
% MCHAIN-PROMPT Toolbox for MATLAB

% By Gaik Tamazian (2014), Elizabeth Lutsenko (2016).

if ~iscell(pdbStruct)
    pdbStruct = {pdbStruct};
end

nTrans = length(pdbStruct);
nModels = length(pdbStruct{1}.Model);
rmsdValues = zeros(nTrans, nModels-1);

for i = 1:nTrans 
    coords = pdbextractcoords(pdbStruct{i});
    for j = 2:nModels 
        % superpose the current model to the previous one
        [~, coords{j}] = procrustes(coords{j-1}, coords{j}, ...
            'scaling', false, 'reflection', false);
        rmsdValues(i,j-1) = mean(sqrt(sum((coords{j} - ...
            coords{j-1}).^2,2)));
    end
end

h = plot(transpose(rmsdValues),'-o');
xlabel('Configuration Pair');
ylabel('RMSD in A');

% modify x axis tick labels
ax = gca;
xticks = get(ax,'XTickLabel');
for j = 1:nModels-1
    xticks{j} = [int2str(j),'-',int2str(j+1)];
end
set(ax,'XTickLabel',xticks);

% plot will be saved if a name is provided for it
if exist('name', 'var')
    name = strcat(name, '_adj');
    print(name, '-dpng');
end

end

