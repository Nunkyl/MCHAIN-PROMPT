function transformations = trmvariateinterpolations(trmodel, nAngles)
%TRMVARIATEINTERPOLATIONS Vary interpolation modes for torsion angles
%   TRMVARIATEINTERPOLATIONS(trmodel, nAngles) returns a cell array of 
%   2^nAngles transformations which differ in the direction the most 
%   differing angles in the specified transformation trmodel were 
%   interpolated.
%
%   See also trmcreate
%
% MCHAIN-PROMPT Toolbox for MATLAB

% By Gaik Tamazian (2014), Elizabeth Lutsenko (2016).

if ~iscell(trmodel) && ~isstruct(trmodel)
    trmodel = {trmodel};
end

nChains = size(trmodel, 2);

if (~exist('nAngles', 'var'))
    nAngles = 2*ones(1, nChains);
end

nTransformations = sum(nAngles);

transformations = cell(2^nTransformations,1);
interpmask = 1 - getbinarycode([],nTransformations);
nConf = size(trmodel(1).r, 2);

% get indices of N torsion angles which differ the most in the first and
% last configurations of the specified transformation

I = trmdistantangleindices(trmodel, nAngles);

for i = 1:2^nTransformations
    newModel = trmodel;
    for j=1:nChains
        newModel(j).psi(I{j},:) = circinterp(newModel(j).psi(I{j},1), ...
            newModel(j).psi(I{j},end), nConf-2, ...
                transpose(interpmask(i, sum(nAngles(1:j-1)) + ...
                    1:sum(nAngles(1:j-1))+nAngles(j))));
        newModel(j).psi = reduceangles(newModel(j).psi);
    end
    transformations{i} = newModel;
end


end

