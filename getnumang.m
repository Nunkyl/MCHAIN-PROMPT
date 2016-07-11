function [nAngles] = getnumang(model, angles)
% GETNUMANG returns the number of angles that have to be interpolated in
% different ways (long or short arc) for each chain.
%   *model* is transformation model
%   *angles* is a boarder value of psi; angles that are greater than it  
%   should be interpolated (not nec)  
%
% MCHAIN-PROMPT toolbox for MATLAB

% By Elizabeth Lutsenko, 2016.

if ~iscell(model) && ~isstruct(model)
    model = {model};
end

nChains = size(model, 2);

if (~exist('angles', 'var'))
    angles = 140*ones(1, nChains);
end

nAngles = zeros(1, nChains);

for i=1:nChains
    nAngles(i) = sum(abs(circdist(model(i).psi(:,1), ...
        model(i).psi(:,end)))*180/pi >= angles(i));
end

end

