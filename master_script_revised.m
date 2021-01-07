% This master script is used to run the javier model, schisto simulation,
% optimization, and analysis of G0

clear all; close all;
%%
network_type = 'random_geo'; % change the network topology here
population_opt = 'javier'; % change the network population distribution here
isNetwork = true; %true if we care about road networks, false if want fully connected network


% load constants
javierConstants
    
n = 100; % number of nodes (may be overwritten, depending on network topology

% parameters used if you want to specifiy different m for each node
cities = [];
m_max = 1;

switch network_type  
    case 'lattice'
    % lattice network, requires that n is a perfect square
        n = 100;
        dist = 10;
        [loc, spaths, A] = create_lattice_network(n, dist);
    case 'random'
    % generates a network with random locations and random connections. Two
    % nodes next to each other in space are not necessarily connected
        n_add_edges = 150;
        x1 = 0; 
        x2 = 100;
        x = x1 + (x2-x1).*rand(n,1);
        y = x1 + (x2-x1).*rand(n,1);     
        loc = [x y];
        [loc, spaths, A] = create_random_network(n, x1, x2, n_add_edges, loc);
    case 'small_world'
    % circular small world phenomenon
        p = 0.1;
        q =0.1;
        dist = 50;
        [loc, spaths, A] =create_smworld_network(n, p, q, dist);
    case 'nav_small_world'
    % lattice small world
        filename =  '../../data/final_networks/small_world_network.txt';
        [loc, spaths, A, n] = graph_nav_small_world(filename);
    case 'havel_hakimi' 
    % generate network from degree sequence using havel hakimi algorithm
        degree_seq = 3 * ones(1,n);
        filename = '../../data/final_networks/havel_hakimi.txt';
        [loc, spaths, A, n, centrality] = graph_from_degree_sequence(filename, 'false', []);
    case 'geo_threshold'
        filename = '../../data/final_networks/geo_threshold_network.txt';
        [loc, spaths, A, n] = graph_from_degree_sequence(filename, 'false', []);
    case 'random_geo'
        %filename = '../data/final_networks/random_geo_network.txt';
        filename = '../../data/networks_jonas/random_geo_network_javier.txt'
        [loc, spaths, A, n] = graph_geo_network(filename);
    case 'GUI'
    % Generates a network that is created using the Java GUI tool for
    % network building
        GUIfilename = '../../data/networks_jonas/simplenet_urban.txt';
        %other options for files:
        %GUIfilename = '../data/final_networks/effect_ruralconnectivity/pptedgelist.txt';
        %GUIfilename = '../data/java_networks/edgelist10.26.txt';
        %GUIfilename = '../data/java_networks/edgelist.txt';

        [loc, spaths, A, n] = create_GUI_network(GUIfilename, true); 
    case 'GUI_no_edge'
        GUIfilename = '../../data/java_networks/edgelist.txt';
        [loc, spaths, A, n] = create_GUI_network(GUIfilename, false); 
    case 'random_degree_seq'
        rdsfilename = '../../data/final_networks/random_degree_seq.txt';
        [loc, spaths, A, n, centrality] = graph_from_degree_sequence(rdsfilename, 'true', []);
    case 'expected_degree'
        n = 10;
        x1 = 0;
        x2 = 100;
        x = x1 + (x2-x1).*rand(n,1);
        y = x1 + (x2-x1).*rand(n,1);     
        loc = [x y];
        loc(1,:) = [50 50];
        loc(2,:) = [40 55];
        loc(3,:) = [45 35];
        filename = '../../data/final_networks/expected_degree/expected_degree_graph_1.txt';


        [loc, spaths, A, n, centrality] = graph_from_degree_sequence(filename, 'true', loc);

    case 'sparseness'
        n = 100;
        sparseness = (n*n - n)/2.2;
        x1 = 0;
        x2 = 100;
        x = x1 + (x2-x1).*rand(n,1);
        y = x1 + (x2-x1).*rand(n,1);     
        loc = [x y];
 
        [loc, spaths, A, n] = vary_sparseness(n, sparseness, loc);
    otherwise % default
        n_add_edges = 200;
        x = x1 + (x2-x1).*rand(n,1);
        y = x1 + (x2-x1).*rand(n,1);     
        loc = [x y];
        [loc, spaths, A] = create_random_network(n, 0, 100, n_add_edges, loc);
end

N = 5000 * ones(1,n); % snail populations default

%% Options for population distributions 
% Calculating populations at each node
switch population_opt
    case 'zipf' 
    % Zipf Distribution
        % Rank locations by distance from center, Burkina Faso
        center = [12.36566, -1.53388];
        dists = sqrt(sum((loc - center).^2,2));
        [~, ind] = sort(dists);
        a = 1e6;  % 2e5;
        k = 2;
        zipf = zipf_population(n, a, k);
        H(ind) = zipf;  % sort(zipf, 'descend');
        %plot([1:n],sort(H)) %plot resulting populations
        n_cities = sum(sum(H > 50000));
        n_towns = sum(sum(H > 10000)) - n_cities;
        n_villages = sum(sum(H < 2000));
        N = 6000 * ones(1, n);
        % N = [300 * ones(1,n_cities) 1000 * ones(1,n_towns) 2500 * ones(1,n-n_cities-n_towns-n_villages) 6000 * ones(1,n_villages)]; % snail populations
        % N = [max(normrnd(500, 100, 1, n_cities), 0) max(normrnd(1000, 200, 1, n_towns), 0) max(normrnd(2500, 500, 1, n-n_cities-n_towns-n_villages), 0) max(normrnd(6000, 1000, 1, n_villages), 0)];
        % rand_i = randperm(length(H));
        % H=H(rand_i);
        % N=N(rand_i);
    case 'rural' 
    % Rural setting where population is drawn randomly 
    % from uniform distribution between x1 and x2
        seed = 1;
        x1 = 750;
        x2 = 2000;
        y1 = 5000;
        y2 = 7000;
        s = rng;
        H = x1 + (x2 - x1) .* rand(1,n);
        rng(s);
        % N = y1 + (y2 - y1) .* rand(1,n); %6000 * ones(1,n); % snail populations
        N = max(normrnd((y1+y2)/2, 1000, 1, n), 0);

    case 'identical'
        % Every node has the same population
        H = 1000 * ones(1,n);
    case 'urban centers'
        %H = [75000 normrnd(5000,100,1,2) 75000 normrnd(5000,100,1,11)];
        %N = [1300 2400*ones(1,2) 1300 2400*ones(1,11)];
        H = [75000 normrnd(5000,100,1,2) 75000 normrnd(5000,100,1,11)];
        N = [500 2500*ones(1,2) 500 2500*ones(1,11)];
    case 'high urban centers'
        %H = [75000 normrnd(5000,100,1,2) 75000 normrnd(5000,100,1,11)];
        %N = [1300 2400*ones(1,2) 1300 2400*ones(1,11)];
        H = [75000 normrnd(5000,100,1,2) 75000 normrnd(5000,100,1,11)];
        N = [500 4000*ones(1,2) 500 4000*ones(1,11)];
    case 'city'
    % Generates a range of population distributions by keeping cities = 5%
    % of total number of nodes, but varying % makeup of towns
        allH = [];
        allN  =[];
        city_percent = [.05; .05; .05; .05; .05];
        %city_percent = [.01; .1; .2; .4];

        town_percent = [.05; .1; 0.3; 0.5; 0.7];
        %town_percent = [.2; .2; 0.2; 0.2];
        rand_i = randperm(n);

        for w = 1:length(city_percent)
            n_cities = round(n.*city_percent(w));
            n_towns = round(n.*town_percent(w));
            n_villages = n - n_cities - n_towns;
%             H = [200000*ones(1,n_cities) 50000*ones(1,n_towns) 2000*ones(1,n_villages)];
%             N = [400*ones(1,n_cities) 800*ones(1,n_towns) 3000*ones(1,n_villages)];  
            H = [normrnd(200000,10000,1,n_cities) normrnd(50000,5000,1,n_towns) normrnd(2000,100, 1,n_villages)];
            N =  [normrnd(75,20, 1,n_cities) normrnd(300,50, 1,n_towns) normrnd(1000,100, 1,n_villages)];
            H=H(rand_i);
            N=N(rand_i);
            allH = [allH; H];
            allN = [allN; N];
        end
        varyPop = town_percent;
            
    case 'normal dist'
    % Generates human and snail populations according to a normal
    % distribution. The number of large cities, cities, and towns can be
    % specified to make up a certain percentage of the population.
        n_lg_cities = 2;
        n_cities = round(n*.05);
        n_towns = round(n*.15);
        n_villages = n - n_cities - n_towns- n_lg_cities;
        H = [500000*ones(1,n_lg_cities) normrnd(200000,10000,1,n_cities) normrnd(50000,5000,1,n_towns) normrnd(2000,100, 1,n_villages)];
        N = [50*ones(1,n_lg_cities) normrnd(75,20, 1,n_cities) normrnd(300,50, 1,n_towns) normrnd(1000,100, 1,n_villages)];
        rand_i = randperm(length(H));
        H=H(rand_i);
        N=N(rand_i);
        plot(sort(H))
    case 'custom'
        % Rapidly changes depending on the test, used for short term
        % hypotheses
        seed = 1;
        x1 = 750;
        x2 = 2000;
        y1 = 5000;
        y2 = 7000;
        s = rng;
        p = unifrnd(0.5, 1.5) * ones(1,n);
        H = 1500 * ones(1,n);
        rng(s);
        N = 6000 * ones(1,n) .* p; % [5 5 5 5 5 6 6 7 7 7] * 1000;
    case 'gauss'
    % Vary population distribution by using a gaussian distribution with a
    % range of variances to determine the effect on G0. Use this with
    % vary_populations function at bottom
        allH = [];
        allN  =[];
        s = [300, 500, 700, 1000, 1500, 1700, 2000];
        for t = 1:length(s)
            H = abs(normrnd(3000,s(t),1,n));
            allH = [allH; H];
            allN = [allN; 800 * ones(1,n)];
        end
        varyPop = s';
    case 'populations'
    % Vary population distribution by using a zipf distribution with a
    % range of k values to determine the effect on G0. Use this with
    % vary_populations function at bottom
        allH = [];
        allN = [];
        a = 2e6;
        k = [ 0.7; 0.8; 0.9; 1.0; 1.2; 1.3; 1.7; 1.9; 2.1];% 1.7; 2.0]
        rand_i = randperm(length(n));
        for t = 1:length(k)
            H = zipf_population(n, a, k(t));
            n_cities = sum(sum(H > 100000));
            N = [100 * ones(1,n_cities) 1200 * ones(1,n-n_cities)]; % snail populations
            H=H(rand_i);
            N=N(rand_i);
            allH = [allH; H];
            allN = [allN; N];

        end
        varyPop = k;
    case 'senegal' %n = 697, 557 without zeros
    % Reads in the Senegal population data that Sanna gave me, which
    % contains population info for all villages within 10km of water
        [num, txt, raw] = xlsread('../../data/networks_jonas/Senegal_populations_filtered.xls', 1);
        H = num(2:end,3)';
        n = length(H);
        plot(sort(H))
        N = [5000 * ones(1, n)]; % snail populations
    case 'javier'
        raw = csvread('../../data/networks_jonas/burkinafaso_javier.csv', 1, 0);
        H = raw(:,3)';
        n = length(H);
        plot(sort(H))
        snail_densities = raw(:,4)';
        N = 6000 * ones(size(H));
        % N = 10000 * snail_densities; % snail populations

end

V = ones(1,n); %water volume, set to 1 because it isn't considered in this model
theta = exposure_rate(H, H_trans, alpha, theta_urb, theta_rur); % calculate exposure rate as a function of population
thetap = theta; % we assume contamination rate = exposure rate, which isn't necessarily true

%%
%  Calculate Mean Worm Burden at Each Node 
m = linspace(0,1,5);
total_wb = [];
total_ind_hburden = [];
max_n_worms = 20;
k = .25;
thresh = 20;
%for mob = 1:length(m)
   mob=length(m);
   m(mob) = 0.3;
    mobility = m(mob) * ones(n,1); % mobility
    W0 = 5*ones(1,n);
    Y0 = 0.05*ones(1,n);
    C0 = 80*ones(1,n);
    M0 = 100*ones(1,n);
    [Y_init, t_init, R0s_init] = schisto_model_revised('none','none','none', n, N', H', loc, spaths, A, mobility, isNetwork, W0, Y0, C0, M0, k, thresh, 0);    
    worms = reshape(Y_init(:,1,:),length(t_init),n);
    snails = reshape(Y_init(:,2,:),length(t_init),n);
    cercariae = reshape(Y_init(:,3,:),length(t_init),n);
    miracidia = reshape(Y_init(:,4,:),length(t_init),n);
    W0 = worms(end,:);
    Y0 = snails(end,:);
    C0 = cercariae(end,:);
    M0 = miracidia(end,:);
    [Y1, t1, treated_stats1, nmwb1] = schisto_model_revised('const','none','none', n, N', H', loc, spaths, A, mobility, isNetwork, W0, Y0, C0, M0, k, thresh, 1);        
    [Y2, t2, treated_stats2, nmwb2] = schisto_model_revised('prop','none','none', n, N', H', loc, spaths, A, mobility, isNetwork, W0, Y0, C0, M0, k, thresh, 1);
    [Y3, t3, treated_stats3, nmwb3] = schisto_model_revised('bc','none','none', n, N', H', loc, spaths, A, mobility, isNetwork, W0, Y0, C0, M0, k, thresh, 1);
    [Y4, t4, treated_stats4, nmwb4] = schisto_model_revised('none','none','none', n, N', H', loc, spaths, A, mobility, isNetwork, W0, Y0, C0, M0, k, thresh, 1);
    
    nmwb_perc_rec100 = 100*[(nmwb1(1)-nmwb1(101))/nmwb1(1) (nmwb2(1)-nmwb2(101))/nmwb2(1) (nmwb3(1)-nmwb3(101))/nmwb3(1) (nmwb4(1)-nmwb4(101))/nmwb4(1)];
    nmwb_perc_rec20 = 100*[(nmwb1(1)-nmwb1(21))/nmwb1(1) (nmwb2(1)-nmwb2(21))/nmwb2(1) (nmwb3(1)-nmwb3(21))/nmwb3(1) (nmwb4(1)-nmwb4(21))/nmwb4(1)];
    nmwb_perc_rec10 = 100*[(nmwb1(1)-nmwb1(11))/nmwb1(1) (nmwb2(1)-nmwb2(11))/nmwb2(1) (nmwb3(1)-nmwb3(11))/nmwb3(1) (nmwb4(1)-nmwb4(11))/nmwb4(1)];
    nmwb_perc_rec05 = 100*[(nmwb1(1)-nmwb1(6))/nmwb1(1) (nmwb2(1)-nmwb2(6))/nmwb2(1) (nmwb3(1)-nmwb3(6))/nmwb3(1) (nmwb4(1)-nmwb4(6))/nmwb4(1)];
    labels = {'Const';'Prop';'BC';'None'};
    
    figure; bar(nmwb_perc_rec10, 'k')
    set(gca,'XTickLabel',labels,'XTick',1:4)
    xlabel('Treatment strategy')
    ylabel('Percent reduction in mean worm burden')
    title('Reduction in network mean worm burden, 10 years')
    
    [ R0, G0s, Q2 ] = runjavier_network( n, N, H, loc, spaths, A, cities, m_max); 
    
%     figure; bar(nmwb_perc_rec05, 'k')
%     set(gca,'XTickLabel',labels,'XTick',1:4)
%     xlabel('Treatment strategy')
%     ylabel('Percent reduction in mean worm burden')
%     title('Reduction in network mean worm burden, 5 years')

%     mean_worm_burden = reshape(Y(:,1,:),length(t),n);
%     mwb = mean_worm_burden(end,:);
%     ind_high_burden = n_heavily_infected(H,mwb,k,thresh);
%     ind_med_infection = n_heavily_infected(H,mwb,k,5); %thresh = 1
%% calculate cumulative TWB over time for all nodes
%   [twb_1yr, twb_5yr, twb_10yr] = calculate_total_twb(n, t, Y, H, loc, A, m(mob));

%%
%     %ind_high_burden(ind_high_burden < 1) = 0;
%     figure; bar([1:n],ind_high_burden)
%     title1 = sprintf('Number of Individuals with > 20 worms at Each Node at m = %d', m(mob));
%     title (title1)
%     xlabel('Node')
%     ylabel('Number of Individuals')
%     
%      %ind_high_burden(ind_high_burden < 1) = 0;
%     figure; bar([1:n],ind_med_infection)
%     title4 = sprintf('Number of Individuals with > 5 worms at Each Node at m = %d', m(mob));
%     title (title4)
%     xlabel('Node')
%     ylabel('Number of Individuals')
%     
%     total_worm_burden = mwb.*H;
%     %total_worm_burden(total_worm_burden <1) = 0;
%     figure; bar([1:n], total_worm_burden); 
%     title2 = sprintf('Equilibrium Total Worm Burden at Each Node at m = %d',  m(mob));
%     title(title2)
%     xlabel('Node')
%     ylabel('Total Worm Burden')
    
%     figure; bar([1:n], mwb); 
%     title3 = sprintf('Equilibrium Mean Worm Burden at Each Node at m = %d',  m(mob));
%     title(title3)
%     xlabel('Node')
%     ylabel('Mean Worm Burden')
%     
%     total_wb = [total_wb; sum(mwb.*H)];
%     total_ind_hburden = [total_ind_hburden; sum(ind_high_burden)];
% 
%     graph_mwb(loc, H, mwb, n, A, m(mob)); % map of nodes in space, color coded according to MWB
%end
% 
%% CONTROL
 %controlSimulation(n, mwb, H, N, loc, spaths, A, mobility, isNetwork)
%%
% plot mobility vs TWB 
% figure;
% plot(m, total_wb);
% xlabel('Mobility')
% ylabel('Total Worm Burden')
% title('Mobility''s Effect on Total Worm Burden Across the Network')
% 
% % plot mobility vs TWB above a certain threshold, thresh
% figure;
% plot(m, total_ind_hburden);
% xlabel('Mobility')
% ylabel('Number of Individuals with High Worm Burden')
% title('Mobility''s Effect on # Individuals with Extensive Worm Burden')

%% Generate mobility vs G0 graphs  

% calculate mobility vs G0 graphs for fully connected networks
%[ R0, G0s, Q1 ] = runjavier( n, N, H, loc, cities, m_max ) 

% calculate mobility vs G0 graphs for networks with explicit road topologies
[ R0, G0s, Q2 ] = runjavier_network( n, N, H, loc, spaths, A, cities, m_max); 

%% Determine the effect of sparseness, hubness, and population distribution on mobility vs G0

% determine how sparseness (density) of a network affects mobility vs G0
%figure;
%test_sparseness(n, N, H, cities, m_max) 

% determine how hubness (centralization) of a network affects mobility vs G0
%test_hubness(n, N, H, cities, m_max) 

% determine how population distribution of a network affects mobility vs G0
%test_populations(n,allH,allN, cities, m_max, loc, spaths, A, varyPop) 

%% optimal control
% determine optimal control strategies
% G0cvxsingular( n, N, H, loc, spaths, A, isNetwork) 

