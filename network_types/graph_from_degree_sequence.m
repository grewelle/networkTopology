% Constructing a graph from a given degree sequence: deterministic

function [loc, spaths, A, n, centrality] = graph_from_degree_sequence(filename, containsCentrality, loc)
    data = importdata(filename);
    n = data(1);
    data(1) = [];
    if containsCentrality
        centrality_py = data(1);
        data(1) = [];
    end
    n_edges = length(data)/2;
    edges = unique(reshape(data, 2, n_edges)' + 1, 'rows');
    [loc, A, spaths, centrality] = rand_coord_network(n, edges, loc);

%end