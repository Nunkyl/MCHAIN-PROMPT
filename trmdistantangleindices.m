function angleIndices = trmdistantangleindices(trmodel, nAngles)
%TRMDISTANTANGLEINDICES Get indices of the most distant torsion angles
%   TRMDISTANTANGLEINDICES(trmodel, nAngles) returns indices of nAngles
%   torsion angles which values in the first and last model configurations
%   differ at most.
%
%   See also trmcostangles
%
% MCHAIN-PROMPT Toolbox for MATLAB

% By Gaik Tamazian (2014), Elizabeth Lutsenko (2016).

if ~iscell(trmodel) && ~isstruct(trmodel)
    trmodel = {trmodel};
end

nChains = length(trmodel);
if (~exist('nAngles', 'var'))
    nAngles = zeros(1, nChains);
    for i=1:nChains
        nAngles(i) = floor(2/3*(length(trmodel(i).psi)));
    end
end

angleIndices = cell(1, nChains);

for i=1:nChains
    [~,angleIndices{i}] = sort(abs(circdist(trmodel(i).psi(:,1), ...
        trmodel(i).psi(:,end))), 'descend');
end

for i=1:nChains
    angleIndices{i} = angleIndices{i}(1:nAngles(i));
end

if length(angleIndices) == 1
    angleIndices = angleIndices{1};
end

end


