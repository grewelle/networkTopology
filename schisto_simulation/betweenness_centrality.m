% Calculates the betweenness centrality of each node and ranks nodes by BC
% in descending order.
% Note: assumes shortest paths between any two nodes are unique.

function bc = betweenness_centrality(G, n)
    bc = zeros(n,2);
    for iter = 1:n
        sp = cell(n,n);
        for i = 1:n
            for j = i+1:n
                sp{i,j} = shortestpath(G,i,j);
            end
        end
        bc_temp = zeros(n,2);
        for node = 1:n
            if ismember(node, bc(:,2))
                bc_temp(node,1) = -inf;
                bc_temp(node,2) = node;
                continue
            end
            for i = 1:n
                if i == node
                    continue
                end
                for j = i+1:n 
                    if j == node
                        continue
                    end
                    if ismember(node, sp{i,j})
                        bc_temp(node,1) = bc_temp(node,1) + 1;
                    end
                end
            end
            bc_temp(node,2) = node;
        end
        bc_temp = sortrows(bc_temp, [-1 2]);
        bc(iter,:) = bc_temp(1,:);
    end
end