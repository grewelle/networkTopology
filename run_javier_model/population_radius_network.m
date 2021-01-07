% This method finds s, which refers to matrix where sij is the total population living within a
% radius dij around the origin, excluding the origin and destination
% populations, and dij being the distance between nodes i and j
% loc refers to x,y coords of each node - rows are nodes, col 1 is
% x coord, col 2 is y coord
% n_node refers to number of nodes
% H is a vector of populations at each node
% spaths is 
function S = population_radius_network(n_node, loc, H, spaths)
    %D = squareform(pdist(loc)); use this if you want standard radiation model
    D = spaths; %instead of using euclidean dist, incorporate network into calc distances between each node
    for i = 1:n_node
        dest_nodes = D(i,:);
        for j = 1:n_node
            if (i==j)
                S(i,j) = 0;
            else
                dist = D(i,j);
                nearby_node_i = find(dest_nodes <= dist);
                nearby_populations = H(nearby_node_i);
                % need to exclude population of origin and destination 
                S(i,j) = sum(nearby_populations)-H(i)-H(j);
            end
        end
    end    
end
   