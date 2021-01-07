function [ p ] = binom( k,mu )
% This function spits out the binomial distribution as a function of mean
% and the clumping parameter.

% value of p(0)
p = zeros(51,1);
p(1) = (1+mu/k)^(-k);

% iterate
for x=1:50
    p(x+1) = p(x)*((k+x-1)/x)*(mu/(mu+k));
end


end

