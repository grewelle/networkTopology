% This function returns the Q matrix, where Qij is the probability that a
% person travelling out of node i reaches node j
% radiation model - Javier (PLOS paper)
% Qij = Hi*Hj / (Hi + sij)(Hi + Hj + sij)
% n_node refers to number of nodes
% H refers to a vector of populations at each node
% s refers to matrix where sij is the total population living within a
% radius dij around the origin, excluding the origin and destination
% populations, and dij being the distance between nodes i and j
% Authors: Nick White and Emily Alsentzer, Summer 2015
function Q = radiation_model(n_node, H, s)
    Q = zeros(n_node, n_node);
    for i = 1:n_node
        for j = 1:n_node
            if(i==j) 
                Q(i,j) = 0;
            else
                Q(i,j) = (H(i)*H(j)) / ((H(i)+s(i,j)) * (H(i) + H(j) + s(i,j)));
            end
        end
        Q(i,:) = Q(i,:)./sum(Q(i,:));
    end
end
