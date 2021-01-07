%% Generate random x, y coordinates for all nodes from a uniform distribution,
% Calc distances between the coordinates, create Adjacency matrix with
% distances, and calculate shortest paths along the network

function [loc, A, spaths, centrality] = rand_coord_network(n, edges, loc)
  A = sparse([edges(:,1); edges(:,2)], [edges(:,2); edges(:,1)], [ones(size(edges,1),1); ones(size(edges,1),1)]);
  if size(loc,1) ==0
    x1 = 0;
    x2 = 100;
    x = x1 + (x2-x1).*rand(n,1);
    y = x1 + (x2-x1).*rand(n,1);     
    loc = [x y];
  end
  fullA = full(A);
  centrality = calc_centralization(fullA,n);
  D = squareform(pdist(loc));
  A = A.*D;
  sparseA = sparse(A);
  spaths = calc_shortest_paths(n, sparseA);
end