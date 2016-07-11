function [r, alpha, psi] = createmodel(PDBStruct, curNum, bondInd)
%CREATEMODEL Get transformation bond lengths, planar and torsion angles
%   [r, alpha, psi] = CREATEMODEL(PDBStruct, curNum, bondInd) calculates 
%   bond lengths, planar and torsion angles of the transformation specified 
%   by the structure PDBStruct or of the bonds between chains specified in 
%   bondInd for the conformation curNum. The bond lengths are returned in 
%   the vector r, the planar angles in the vector alpha and the torsion 
%   angles in the vector psi. 
%
%   See also restorecoords trmcreate
%
% MCHAIN-PROMPT Toolbox for MATLAB

% By Gaik Tamazian (2014), Elizabeth Lutsenko (2016).

if ~exist('bondInd', 'var')
    X = [PDBStruct{curNum}.X];
    Y = [PDBStruct{curNum}.Y];
    Z = [PDBStruct{curNum}.Z];
    [r, alpha, psi] = ralphapsi(X, Y, Z);
else
    chainInd1 = bondInd(curNum,1);
    chainInd2 = bondInd(curNum,2);
    point1 = bondInd(curNum,3);
    point2 = bondInd(curNum,4);
    
    chain1 = zeros(length([PDBStruct{chainInd1}.X]),3);
    chain2 = zeros(length([PDBStruct{chainInd2}.X]),3);
    
    chain1(:,1) = [PDBStruct{chainInd1}.X];
    chain1(:,2) = [PDBStruct{chainInd1}.Y];
    chain1(:,3) = [PDBStruct{chainInd1}.Z];
    
    chain2(:,1) = [PDBStruct{chainInd2}.X];
    chain2(:,2) = [PDBStruct{chainInd2}.Y];
    chain2(:,3) = [PDBStruct{chainInd2}.Z];
    
    bond = zeros(6,3);

    if point1 >= 3
        bond(1:3,:) = [chain1(point1-2,:); ...
            chain1(point1-1,:);chain1(point1,:)];
    else
        bond(1:3,:) = [chain1(point1+2,:); ...
            chain1(point1+1,:);chain1(point1,:)];
    end
    
    if point2 >= 3
        bond(4:6,:) = [chain2(point2,:); ...
            chain2(point2-1,:);chain2(point2-2,:)];
    else
        bond(4:6,:) = [chain2(point2,:); ...
            chain2(point2+1,:);chain2(point2+2,:)];
    end
    
    [r, alpha, psi] = ralphapsi(bond(:,1)', bond(:,2)', bond(:,3)');
    
    r = r(3);
    alpha = alpha(2:3);
    psi = psi(1:3);
end

end

function [r, alpha, psi] = ralphapsi(X, Y, Z)
% Calculate bond lengths and planar and torsion angles for each
% configuration.
% In matrices r, alpha and psi, the ith column corresponds to the ith
% configuration.

deltax = diff([X' Y' Z']);

r = sqrt(sum(deltax.^2, 2));

alpha = acos(dot(deltax(1:end-1,:), deltax(2:end,:), 2) ./ ...
    (r(1:end-1) .* r(2:end)));

N = cross(deltax(1:end-1,:), deltax(2:end,:), 2);

psi = atan2(r(2:end-1) .* ...
    dot(deltax(1:end-2,:), N(2:end,:), 2), ...
    dot(N(1:end-1,:), N(2:end,:), 2));
end