function n = n_heavily_infected_revised( H,mwb,k,thresh )
% This function takes in an array of population sizes and mean worm burdens
% and returns the number of heavily infected individuals

% k is clumping parameter
% thresh is the havily infected threshold

p = zeros(size(mwb,1), size(mwb,2));
for i=1:size(mwb,1)
    for j=1:size(mwb,2)
        p(i,j) = heavily_infected(k,mwb(i,j),thresh);
    end
end

n = sum(p.*repmat(H,1,size(mwb,2)))';

end

