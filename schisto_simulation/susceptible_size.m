% Calculates the susceptible size of a graph for a given sequence of nodes
% to treat.

function r = susceptible_size(G, ss, n)
    cc = conncomp(G);
    r = sum(cc==mode(cc));
    for q = 1:n
        G = rmedge(G, ss(q), linspace(1,n,n));
        cc = conncomp(G);
        r = r + sum(cc==mode(cc));
    end
    r = r/(n*(n+1));
end