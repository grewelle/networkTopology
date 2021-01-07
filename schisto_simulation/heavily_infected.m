function [ p ] = heavily_infected( k,mu,thresh )
% This returns the portion of a population that is heavily infected,
% defined by having a certain number of worms determined by the parameter
% thresh.

% The distribution is negative binomial
ps = binom(k,mu);

% Sum all probabilities greater than the threshold
p = sum(ps(thresh+1:end));


end

