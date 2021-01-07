%format of the file that is read in:

% n
% x1 y1 x2 y2    each line corresponds to the start and end points of an edge
% x3 y3 x4 y4
% ...

function [loc, spaths, A, n] = graph_nav_small_world(filename)
    filename = '/Users/emilyalsentzer/Documents/DeLeo/small_world_network.txt'
    data = importdata(filename, ' ') %reads all of the data into a single vector
    n = data(1)^2; %first line of textfile is the square root of the number of nodes
    data(1) = [];
    data = data + 1;
    n_edges = length(data)/4;
    locx = data(1:2:end); %every other element in data is an x coord
    locy = data(2:2:end); %every other element in data is a y coord
    loc = unique([locx locy], 'rows');
    loc = sortrows(loc);
    edges = reshape(data, 4, n_edges)';
    edge_list = [];
    for i = 1:size(edges,1)
         edge1 = sqrt(n)*(edges(i,1)-1) + edges(i,2); %calculates the edge number given the x and y coords
         edge2 = sqrt(n)*(edges(i,3)-1) + edges(i,4);
         edge_list = [edge_list; edge1 edge2];
    end
    edge_list = unique(edge_list, 'rows');
    %A = sparse([edge_list(:,1); edge_list(:,2)], [edge_list(:,2); edge_list(:,1)], [ones(size(edge_list,1),1); ones(size(edge_list,1),1)]);
    A = sparse([edge_list(:,1)], [edge_list(:,2)], [ones(size(edge_list,1),1)]);
    D = squareform(pdist(loc));
    A = A.*D;
    sparseA = sparse(A);
    spaths = calc_shortest_paths(n, sparseA);
end

