function [] = trmplotadjrmsd(trmodels)
%TRMPLOTADJRMSD Plot RMSDs between adjacent configurations
%   TRMPLOTADJRMSD(trmodels) plots RMSDs between adjacent configurations of
%   transformation models specified in a cell array trmodels.
%
%   See also trmplotfixedrmsd trmplottranglediff
%
% MCHAIN-PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2014.

if ~iscell(trmodels) && ~isstruct(trmodels)
    trmodels = {trmodels};
end

nTrans = length(trmodels);
nConf = size(trmodels{1}(1).r,2);
rmsdValues = zeros(nTrans,nConf-1);

for i = 1:nTrans
    coords = trmrestorecoords(trmodels{i});
    for j = 2:nConf
        % superpose the current configuration to the previous one
        [~, coords{j}] = procrustes(coords{j-1}, coords{j}, ...
            'scaling', false, 'reflection', false);
        rmsdValues(i,j-1) = mean(sqrt(sum((coords{j} - ...
            coords{j-1}).^2,2)));
    end
end

plot(transpose(rmsdValues),'-o');
xlabel('Configuration Pair');
ylabel('RMSD in AA');

% modify x axis tick labels
ax = gca;
xticks = get(ax,'XTickLabel');
for j = 1:nConf-1
    xticks{j} = [int2str(j),'-',int2str(j+1)];
end
set(ax,'XTickLabel',xticks);

end

