% Given n, the number of nodes, and sparseA, a sparse representation of the
% adjacency matrix and distances between each node, this function will
% calculate the shortest paths between each node in the graph. 

% The sparse adjacency matrix will look something like this:
%   (15,1)      42.8072
%   (14,2)      41.6724
%   (38,3)      34.2849
%    (6,4)      36.5288
%    (8,4)      60.8937
%   (21,4)     102.4177
%   (15,5)      21.5090
%   (32,5)      31.7474

% Author: Emily Alsentzer, Summer 2015


function spaths = calc_shortest_paths(n, sparseA)
   shortest_paths = zeros(n,n);
   for i = 1:n
       [s_dist,path,pred] = graphshortestpath(sparseA,i,'directed',false);
       shortest_paths(i,:) = s_dist;
   end
   spaths = shortest_paths;
end