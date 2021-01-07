%  Input:
%    pos - N x dim matrix of position of nodes, OR, if a number N is given
%          instead, then the N nodes are given positions in a ring
%    p - connection probability OR a vector of length N denoting the
%        probabilities of each number inputs
%    q - the rewiring probability [0,1] (note that for large q asymmetry
%        may arise near the diagonal between the above-diagonal and
%        below-diagonal elements)

function [loc, spaths, A] =create_smworld_network(n, p, q, scale)
    [A, loc] = wattsstrogatz(n,p,q, scale);
    dist = squareform(pdist(loc));
    A = A.*dist;
    sparseA = sparse(A);
    spaths = calc_shortest_paths(n, sparseA)

end