% This function takes in a filename where the GUI network text file is
% located, reads it in, and creates the Adjacency matrix, shortest paths
% (spaths) matrix, and x,y coords for the location of every node. It also
% returns n, the number of nodes in the network.
function [loc, spaths, A, n] = create_GUI_network(filename, edgesPresent)
    data = importdata(filename);
    n = data(1);
    data(1) = [];
    loc = data(1:2*n);
    locx = loc(1:2:end);
    locy = loc(2:2:end);
    loc = [locx locy];
    %loc = unique(reshape(loc, 2, n)', 'rows')
    if edgesPresent
        edges = data(2*n + 1:size(data,1));
        n_edges = length(edges)/2;
        edges = unique(reshape(edges, 2, n_edges)', 'rows');
    
        A = sparse([edges(:,1); edges(:,2)], [edges(:,2); edges(:,1)], [ones(size(edges,1),1); ones(size(edges,1),1)]);
        D = squareform(pdist(loc));
        A = A.*D;
        sparseA = sparse(A);
        spaths = calc_shortest_paths(n, sparseA);
        A = full(A);
    else
        A = zeros(n);
        spaths = zeros(n);
    end
 

end
