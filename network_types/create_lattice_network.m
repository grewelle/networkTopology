% This function returns loc, which contains the x,y coords for all of the
% nodes; spaths, which contains the shortest distances between nodes along
% the network paths; and A, the adjacency matrix, which describes the
% connectivity and distances between nodes in the network
function [loc, spaths, A] = create_lattice_network(nodes, dist)
   n = sqrt(nodes);
   m = sqrt(nodes);
   x_coords = linspace(0,dist*(n-1),n);
   y_coords = linspace(0,dist*(m-1),m);
   [x,y] = meshgrid(x_coords,y_coords);
   x = reshape(x,n*m, 1);
   y = reshape(y,n*m, 1);
   loc = [x y];
   D = squareform(pdist(loc));
   A = D == dist;
   A = A.*dist;
   sparseA = sparse(A);
   spaths = calc_shortest_paths(nodes, sparseA);

end
