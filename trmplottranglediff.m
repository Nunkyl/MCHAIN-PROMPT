function [] = trmplottranglediff(trmodel, sortAngles, angleNum, markerType,name)
%TRMPLOTTRANGLEDIFF Plot differences between transformation torsion angles
%   TRMPLOTTRANGLEDIFF(trmodel,isSorted,angleNum,markerType) plots absolute 
%   values of circular differences between torsion angles of the first and 
%   last configurations of the specified transformation model trmodel. If 
%   sortAngles is true, then the plotted values are sorted in the 
%   descending order. By default, sortAngles is set to false. The parameter
%   angleNum specifies the number of angles to be plotted for each chain; 
%   by default, all angles are shown. The parameter markerType specifies 
%   the plot marker style; its default value is '-'.
%
%   See also trmplotadjrmsd trmplotfixedrmsd
%
% MCHAIN-PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2014.

if ~iscell(trmodel) && ~isstruct(trmodel)
    trmodel = {trmodel};
end

if nargin < 2                
    sortAngles = false;
end

if nargin < 3
    for i=1:length(trmodel)
        angleNum(i) = size(trmodel(i).psi,1);
    end
end

if nargin < 4
    markerType = '-';
end

torsAngleDiff = cell(1,length(trmodel));
for i=1:length(trmodel)
    torsAngleDiff{i} = abs(circdist(trmodel(i).psi(:,1), ...
        trmodel(i).psi(:,end)));
    if sortAngles
        torsAngleDiff{i} = sort(torsAngleDiff{i},'descend');
    end
    torsAngleDiff{i} = torsAngleDiff{i}(1:angleNum(i));
end

for i=1:length(trmodel)
    figure;
    plot(torsAngleDiff{i}*180/pi, markerType);
    
    axis([0 length(torsAngleDiff{i}) 0 180]);
    ylabel('Absolute Circular Distance in Degrees');
    if sortAngles
        xlabel('Torsion Angle Rank');
        title(['50 Torsion Angles Differing the Most between the First and ', ...
            'Last Configurations', ' Chain ', int2str(i)]);
    else
        xlabel('Torsion Angle Number');
        title(['Difference Between Torsion Angles of the First and ', ...
            'Last Configurations', ' Chain ', int2str(i)]);
    end
   grid;
   
if exist('name', 'var')
    name = strcat(name, '_ang_dif');
    print(name, '-dpng');
end

end

end
