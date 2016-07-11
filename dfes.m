function [ edges ] = dfes(beg, edges, graph)
%DFES Depth first edge search
%   Forms a list of bonds between chains. Later these bonds will be created
%   and reconstructed in this list order.
%
% MCHAIN-PROMPT toolbox for MATLAB

% By Elizabeth Lutsenko, 2016.

for i=1:size(graph,2)
    if graph(beg,i) ~= 0 && ~belongs(edges,beg,i)
        edges = [edges; [beg,i]];
        edges = dfes( i, edges, graph );
    end
end

end

function [answ] = belongs (edges, beg,i)
% Determines whether this particular edge [beg,i] has already been found.
    answ = false;
    if sum(sum((edges == beg) + (edges == i),2) == 2) > 0
        answ = true;
    end
    
end