function [ lowerBound, upperBound ] = getbounds(trmodel, angleIndices)
% GETBOUNDS forms the constraints for fmincon.
%   Forms the constraints for fmincon. Bond lengths can change within 20% 
%   of the initial length. All angles must be between -180 and 180 degrees. 
%
% MCHAIN-PROMPT tooolbox for MATLAB

% By Elizabeth Lutsenko, 2016.

if ~iscell(trmodel) && ~isstruct(trmodel)
    trmodel = {trmodel};
end

if ~iscell(angleIndices)
    angleIndices = {angleIndices};
end

nChains = length(trmodel);
nConfs = size(trmodel(1).r, 2) - 2;

if nChains == 1
    n = length(angleIndices{1}) * nConfs;
    lowerBound = -pi*ones(n,1);
    upperBound = pi*ones(n,1);
    return
end

minRadius = zeros(nChains-1);
maxRadius = zeros(nChains-1);

for i=1:nChains-1
    if trmodel(i).bond_r(1,1) > trmodel(i).bond_r(1,end)
        minRadius(i) = trmodel(i).bond_r(1,end);
        maxRadius(i) = trmodel(i).bond_r(1,1);
    else
        minRadius(i) = trmodel(i).bond_r(1,1);
        maxRadius(i) = trmodel(i).bond_r(1,end);
    end
end

lowerRadius = minRadius - minRadius*0.2;
upperRadius = maxRadius + maxRadius*0.2;

len = 0;
for i=1:nChains-1
    len = len + length(angleIndices{i}) + size(trmodel(i).bond_r,1) + ...
        size(trmodel(i).bond_alpha,1) + size(trmodel(i).bond_psi,1);
end
len = len + length(angleIndices{nChains});

lowerBounds = ones(len, nConfs);
upperBounds = ones(len, nConfs);

%include r
for i=1:nChains-1
    lowerBounds(i,:) = lowerBounds(i,:) *  lowerRadius(i);
    upperBounds(i,:) = upperBounds(i,:) *  upperRadius(i);
end

%include angles
lowerBounds(nChains:end,:) = lowerBounds(nChains:end,:) * (-pi);
upperBounds(nChains:end,:) = upperBounds(nChains:end,:) * pi;

%stretch everything into one column
lowerBounds = lowerBounds';
upperBounds = upperBounds';
lowerBound = lowerBounds(:);
upperBound = upperBounds(:);

end
   
