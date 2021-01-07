% This method finds s, which refers to matrix where sij is the total population living within a
% radius dij around the origin, excluding the origin and destination
% populations, and dij being the distance between nodes i and j
% loc refers to x,y coords of each node - rows are nodes, col 1 is
% x coord, col 2 is y coord
% n_node refers to number of nodes
function S = population_radius(n_node, loc, H)
    D = squareform(pdist(loc));
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
   