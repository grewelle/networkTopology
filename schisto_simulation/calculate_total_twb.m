% This function calculates the total twb for all nodes in a network
% over a defined time period. This essentially measures the burden of
% disease over 1, 5, and 10 year periods

function [total_twb_1,total_twb_2, total_twb_3, total_twb_5, total_twb_10, total_infected_1,total_infected_2, total_infected_3, total_infected_5, total_infected_10] = calculate_total_twb(n, t, Y, H, loc, A, m)
    mean_worm_burden = reshape(Y(:,1,:),length(t),n);
%     graph_mwb(loc, H, mean_worm_burden(1,:), n, A, m); % map of nodes in space, color coded according to MWB
%     figure; bar([1:n], mean_worm_burden(end,:)); 
%     title3 = sprintf('Equilibrium Mean Worm Burden at Each Node at m = %d',  m);
%     title(title3)
%     xlabel('Node')
%     ylabel('Mean Worm Burden')
    
    % # infected people at each time step for each node
    thresh = 2;
    k = .25;
    stepwise_nInfected = zeros(95,n); 
    for time = 1:length(t)
        stepwise_nInfected(time,:) = n_heavily_infected(H,mean_worm_burden(time,:),k,thresh);
    end

    % total worm burden at each time step for each node
    twb = mean_worm_burden;
    for row = 1:n
        twb(row,:) = mean_worm_burden(row,:).* H;
    end
    
%     figure;
%     plot(t,twb);
%     hold on; line([365*1,365*1],[0,4e5],'Color','red');
%     hold on; line([365*5,365*5],[0,4e5],'Color','blue');
%     hold on; line([365*10,365*10],[0,4e5],'Color','green');
%     
    index_1yr = find(t > 365);
    index_1yr = index_1yr(1);
    index_2yr = find(t > 365*2);
    index_2yr = index_2yr(1);
    index_3yr = find(t > 365*3);
    index_3yr = index_3yr(1);
    index_5yr = find(t > 365*5);
    index_5yr = index_5yr(1);
    index_10yr = find(t > 365*10);
    index_10yr = index_10yr(1);

    total_twb_1 = 0;
    total_infected_1 = 0;
    for node = 1:n
        node_twb = trapz(t(1:index_1yr),twb(1:index_1yr,node));
        node_nInfected = trapz(t(1:index_1yr), stepwise_nInfected(1:index_1yr,node));
        total_twb_1 = total_twb_1 + node_twb;
        total_infected_1 = total_infected_1 + node_nInfected;
    end
    
    total_twb_2 = 0;
    total_infected_2 = 0;
    for node = 1:n
        node_twb2 = trapz(t(1:index_2yr),twb(1:index_2yr,node));
        node_nInfected2 = trapz(t(1:index_2yr), stepwise_nInfected(1:index_2yr,node));
        total_twb_2 = total_twb_2 + node_twb2;
        total_infected_2 = total_infected_2 + node_nInfected2;
    end
    
    
    total_twb_3 = 0;
    total_infected_3 = 0;
    for node = 1:n
        node_twb3 = trapz(t(1:index_3yr),twb(1:index_3yr,node));
        node_nInfected3 = trapz(t(1:index_3yr), stepwise_nInfected(1:index_3yr,node));
        total_twb_3 = total_twb_3 + node_twb3;
        total_infected_3 = total_infected_3 + node_nInfected3;
    end
 
    total_twb_5 = 0;
    total_infected_5 = 0;
    for node = 1:n
        node_twb5 = trapz(t(1:index_5yr),twb(1:index_5yr,node));
        node_nInfected5 = trapz(t(1:index_5yr), stepwise_nInfected(1:index_5yr,node));
        total_twb_5 = total_twb_5 + node_twb5;
        total_infected_5 = total_infected_5 + node_nInfected5;
    end
 
    total_twb_10 = 0;
    total_infected_10 = 0;
    for node = 1:n
        node_twb10 = trapz(t(1:index_10yr),twb(1:index_10yr,node));
        node_nInfected10 = trapz(t(1:index_10yr), stepwise_nInfected(1:index_10yr,node));
        total_twb_10 = total_twb_10 + node_twb10;
        total_infected_10 = total_infected_10 + node_nInfected10;
    end
end