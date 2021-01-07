% This function calculates the populations at each node according to a zipf
% distribution. A and K are parameters of the zipf distribution, and scale
% tells you the max population of the largest city

function H = zipf_population(n, a, k)
    x = [1:n];
    H = a.*x.^(-k);
end