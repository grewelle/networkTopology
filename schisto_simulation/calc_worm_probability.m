% Given the mean worm burden and populations at each node, this function
% calculates the number of people who have x number of worms form 0 to
% max_n_worms as well as the number of people with x or more worms at each
% node

% max_n_worms = max number of worms for which to generate probability distribution values 
% mu = mean worm burden
% k = parameter of the negative binomial function
% H = populations at each node
function [worm_prob, cum_worm_prob, indiv_worms, indiv_cum_worms] = calc_worm_probability(max_n_worms, mu, k, H)
    worm_prob = zeros(length(mu),max_n_worms); %keep track of probability of having x number of worms
    cum_worm_prob = zeros(length(mu),max_n_worms); %keep track of cumulative probability

    worm_prob(:,1) =(1+mu/k).^-k;    %probability of having 0 worms
    cum_worm_prob(:,1) = worm_prob(1);
    for i = 2:length(worm_prob(1))
        worm_prob(:,i) = worm_prob(:,i-1).*((k+i-1)/i).*(mu/(mu+k));
        cum_worm_prob(:,i) = cum_worm_prob(:,i-1) + worm_prob(:,i);
    end
    
    size(H)
    size(worm_prob(1,:))
    for j = 1:length(mu)
        indiv_worms(j,:) =  worm_prob(j,:) .* H(j); %number of people with x number of worms
        indiv_cum_worms(j,:) =  (1-cum_worm_prob(j,:)).*H(j); % number of people with x worms or more
    end
end