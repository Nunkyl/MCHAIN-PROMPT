function coordsUnited = restorecoords(model, confNum)
%RESTORECOORDS Restore Cartesian coordinates of atoms
%   restorecoords(r, alpha, psi) restores Cartesian coordinates of the 
%   atoms from the specified bond length, planar and torsion angle vectors.
%
%   The first point is zero, the second point lies on OX axis and the third
%   one lies on XOY coordinate plane. The Natural Extension Reference frame
%   method is used for placing atoms [1].
%
%   References:
%
%       [1] Parsons, J., Holmes, J. B., Rojas, J. M., Tsai, J., & Strauss, 
%       C. E. M. (2005). Practical conversion from torsion space to 
%       Cartesian space for in silico protein synthesis. Journal of %
%       computational chemistry, 26(10), 1063 8. doi:10.1002/jcc.20237
%
%   See also createmodel trmrestorecoords
%
% MCHAIN-PROMPT Toolbox for MATLAB

% By Gaik Tamazian (2014), Elizabeth Lutsenko (2016).


if ~iscell(model) && ~isstruct(model)
    model = {model};
end

nChains = length(model);

coords = cell(1, nChains);
for i=1:nChains
    nAtoms = length(model(i).m);
    coords{i} = zeros(nAtoms, 3);
end

% working with chain A
beg(1,:) = [0 0 0];
beg(2,:) = [model(1).r(1, confNum) 0 0];
beg(3,:) = beg(2,:) + ...
    [model(1).r(2, confNum)*cos(model(1).alpha(1, confNum)) ...
    model(1).r(2, confNum)*sin(model(1).alpha(1, confNum)) 0 ];

r = model(1).r(3:end, confNum);
alpha = model(1).alpha(2:end, confNum);
psi = model(1).psi(:, confNum);

coords{1} = restore (beg, r, alpha, psi);

% working with all the other chains
for curBond=1:nChains-1
    
    chainInd1 = model(1).bondInd(curBond,1);
    chainInd2 = model(1).bondInd(curBond,2);
    point1 = model(1).bondInd(curBond,3);
    point2 = model(1).bondInd(curBond,4);
                 
    if point1 >= 3
        beg = [coords{chainInd1}(point1-2,:); ...
            coords{chainInd1}(point1-1,:); ...
            coords{chainInd1}(point1,:)];
    else
        beg = [coords{chainInd1}(point1+2,:); ...
            coords{chainInd1}(point1+1,:); ...
            coords{chainInd1}(point1,:)];
    end


    if point2 >= 3
        r1 = model(chainInd2).r(1:point2-1,confNum);
        alpha1 = model(chainInd2).alpha(1:point2-2,confNum);
        psi1 = model(chainInd2).psi(1:point2-3,confNum);

        r2 = model(chainInd2).r(point2:end,confNum);
        alpha2 = model(chainInd2).alpha(point2-1:end,confNum);
        psi2 = model(chainInd2).psi(point2-2:end,confNum);

        r1 = [model(curBond).bond_r(:,confNum); flipud(r1)];
        alpha1 = [model(curBond).bond_alpha(:,confNum); flipud(alpha1)];
        psi1 = [model(curBond).bond_psi(:,confNum); flipud(psi1)];

        % restore the first half of the coordinates
        coords1 = restore (beg, r1, alpha1,psi1);
        
        % cut the extra 3 points out
        coords1 = coords1(4:end,:);
        
        % flip the array back
        coords1 = flipud(coords1);

        % restore the second half of the coordinates
        beg = [coords1(point2-2,:);coords1(point2-1,:);coords1(point2,:)];
        
        coords2 = restore (beg, r2, alpha2, psi2);

        % cut the extra 3 points out
        coords2 = coords2(4:end,:);

        % glue together two parts of the chains
        coords{chainInd2} = [coords1; coords2];   

    elseif point2 == 1
        
        % add the bond
        r = [model(curBond).bond_r(:,confNum); ...
            model(chainInd2).r(:,confNum)];
        alpha = [model(curBond).bond_alpha(:,confNum); ...
            model(chainInd2).alpha(:,confNum)];
        psi = [model(curBond).bond_psi(:,confNum); ...
            model(chainInd2).psi(:,confNum)];

        % restore the coordinates
        coords{chainInd2} = restore (beg, r, alpha, psi);

        % cut the extra 3 points out
        coords{chainInd2} = coords{chainInd2}(4:end,:);
        
    else

        % add the bonds between chains
        r = [model(curBond).bond_r(:,confNum); ...
            model(chainInd2).r(2:end,confNum)];
        alpha = [model(curBond).bond_alpha(:,confNum); ...
            model(chainInd2).alpha(2:end,confNum)];
        psi = [model(curBond).bond_psi(:,confNum); ...
            model(chainInd2).psi(2:end,confNum)];

        % restore the coordinates
        coor = restore (beg, r, alpha, psi);
        
        % cut the extra 2 points out
        coor = coor(3:end,:);

        % add the remaining first point
        coor(1,1) = model(chainInd2).r(1, confNum).* ...
            cos(model(chainInd2).alpha(1, confNum));
        coor(1,2) = model(chainInd2).r(1, confNum).* ...
            sin(model(chainInd2).alpha(1, confNum)).* ...
            cos(model(chainInd2).psi(1, confNum));
        coor(1,3) = model(chainInd2).r(1, confNum).* ...
            sin(model(chainInd2).alpha(1, confNum)).* ...
            sin(model(chainInd2).psi(1, confNum));

        bc = coor(2,:) - coor(3,:); bc = bc./norm(bc);
        n = cross3d(coor(3,:) - coor(4,:), bc); n = n./norm(n);
        M = [bc', cross3d(n, bc)', n']; 
        coor(1,:) = (M*coor(1,:)')' + coor(2,:);

        coords{chainInd2} = coor;

    end       
end

coordsUnited = vertcat(coords{:});

end

function [coor] = restore (beg, r, alpha, psi)

        coor = zeros(length(r)+3,3);

        coor(1:3,:) = beg;

        coor(4:end,1) = r(:).* cos(alpha(:));
        coor(4:end,2) = r(:).* sin(alpha(:)).*cos(psi(:));
        coor(4:end,3) = r(:).* sin(alpha(:)).*sin(psi(:));

        for j = 4:length(r)+3
            bc = coor(j-1,:) - coor(j-2,:); bc = bc./norm(bc);
            n = cross3d(coor(j-2,:) - coor(j-3,:), bc); 
            n = n./norm(n);
            M = [bc', cross3d(n, bc)', n']; 
            coor(j,:) = (M*coor(j,:)')' + coor(j-1,:);
        end
end

