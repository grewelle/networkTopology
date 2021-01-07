% Adaptive betweenness centrality:
% Ranks nodes by iteratively calculating betweenness centrality and
% then removing the top-ranked node from the graph

function bc = betweenness_centrality_final(G, n)
    bc = zeros(n, 1);
    for i = 1:n
        bc_ranks = centrality(G, 'betweenness');
        [~, indices] = sort(bc_ranks, 'descend');
        bc(i) = indices(1);
        G = rmedge(G, bc(i), linspace(1,n,n));
    end
end