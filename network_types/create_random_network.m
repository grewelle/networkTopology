
% This function returns the x,y coordinates of a random connected network
% (i.e. one connected component)
% INPUTS: n = number of nodes, n_add_edges = the number of additional
% connections in the network beyond the random spanning tree - this
% controls sparsity of the graph
% OUTPUTS: loc = matrix where each row corresponds to the coords of a node.
% 1st col = x coord, 2nd col = y coord; spaths = distance matrix with the
% shortest distances along the network between each node; A = adjacency
% matrix for network

% To generate the random network:
% 1. Generate random x, y coordinates for all nodes from a uniform
% distribution
% 2. Generate a (uniformly chosen) random spanning tree with N nodes and N - 1 edges.
% 3. Until the requested number of edges has been reached, add an edge between any two random nodes.
% adapted from: https://gist.github.com/bwbaugh/4602818
% http://stackoverflow.com/questions/2041517/random-simple-connected-graph-generation-with-given-sparseness

% Author: Emily Alsentzer, Summer 2015

function [loc, spaths, A] = create_random_network(n,x1, x2, n_add_edges, loc)
 %% 1. Generate random x, y coordinates for all nodes from a uniform distribution

 %% 2. Generate a (uniformly chosen) random spanning tree with n nodes

  node_list =[1:n];
  S = [1:n];
  T = [];
  curr_node_i = randi(length(S));
  curr_node = S(curr_node_i);
  T = [T curr_node];
  S(curr_node_i) = [];
  edges = [0 0];
  while ~isempty(S)
        node_s_i = randi(length(S));
        node_s = S(node_s_i);
        node_t_i = randi(length(T));
        node_t = T(node_t_i);
        if ~ismember([node_t node_s], edges, 'rows')
            edges = [edges; node_s node_t];
            T = [T node_s];
            S(node_s_i) = [];
        end
  end
  edges(1,:)=[]; % get rid of the [0 0] edge that was added as a placeholder 
  
  %% 3. Until the requested number of edges has been reached, add an edge between any two random nodes.

  add_edges = [0 0];
  while size(add_edges,1) < n_add_edges
      randedge_1 = randi(n);
      randedge_2 = randi(n);
      if ~ismember([randedge_1 randedge_2], edges, 'rows') && ~ismember([randedge_1 randedge_2], edges, 'rows') && randedge_1 ~= randedge_2 && ~ismember([randedge_1 randedge_2], add_edges, 'rows') && ~ismember([randedge_2 randedge_1], add_edges, 'rows');
            add_edges = [add_edges; randedge_1 randedge_2];
      end
  end
  add_edges(1,:)=[]; % get rid of the [0 0] edge that was added as a placeholder 

  edges = [edges; add_edges];
  edges = unique(edges, 'rows'); % matrix of edges where each row is an edge, [origin dest]
  
  %% calculate distances between each edge based on their locations
  dist = zeros(size(edges,1),1);
  for i = 1:size(edges,1)
      node1 = edges(i,1);
      node2 = edges(i,2);
      dist(i) = pdist([loc(node1,:); loc(node2,:)]);
  end
  
sparseA = sparse([edges(:,1); edges(:,2)], [edges(:,2); edges(:,1)], [dist; dist])
A = full(sparseA);
spaths = calc_shortest_paths(n, sparseA);
  
end
