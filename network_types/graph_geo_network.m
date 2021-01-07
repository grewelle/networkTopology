    % Constructing a graph from a given degree sequence: deterministic

function [loc, spaths, A, n] = graph_geo_network(filename)
    data = importdata(filename);
    n = data(1);
    data(1) = [];
    loc = data(1:2*n);
    locx = loc(1:2:end);
    locy = loc(2:2:end);
    loc = [locx locy];
    edges = data(2*n + 1:size(data,1));
    n_edges = length(edges)/2;
    edges = unique(reshape(edges, 2, n_edges)' + 1, 'rows');
    A = sparse([edges(:,1); edges(:,2)], [edges(:,2); edges(:,1)], [ones(size(edges,1),1); ones(size(edges,1),1)]);
    D = squareform(pdist(loc));
    A = A.*D;
    sparseA = sparse(A);
    spaths = calc_shortest_paths(n, sparseA);
    A = full(A);

end