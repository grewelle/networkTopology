function [n] = n_heavily_infected( H,mwb,k,thresh )
% This function takes in an array of population sizes and mean worm burdens
% and returns the number of heavily infected individuals

% k is clumping parameter
% thresh is the havily infected threshold

p = zeros(size(mwb));
for i=1:length(mwb)
    p(i) = heavily_infected(k,mwb(i),thresh);
end

n = p.*H;

end

