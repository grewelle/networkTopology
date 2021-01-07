function  controlSimulation(n, mwb, H, N, loc, spaths, A, mobility, isNetwork)

    newstart_mwb = mwb;
    newstart_mwb(newstart_mwb < 0) = 0; %testing!
    newTWB = newstart_mwb.*H;
    k = .25;
    thresh = 0;
    ind_high_burden = n_heavily_infected(H,mwb,k,thresh);
    
    node1 = 1; %13
    node2 = 4; %14
    node3 = 9; %15
    reducePercent = 1;%0.85;
%%

    %% reducing by same TWB in each node
    reducedTWBamt = (newTWB(node1) + newTWB(node2) + newTWB(node3))*reducePercent;
    reduceAllAmt = 1- (reducedTWBamt/sum(newTWB));
    reduce148Amt = 1- (reducedTWBamt/(newTWB(1) + newTWB(4) + newTWB(8)));
    reduce131415Amt = 1- (reducedTWBamt/(newTWB(node1) + newTWB(node2) + newTWB(node3)));
    if(reduce148Amt <0) 
        disp('error')
        reducedTWBamt = (newTWB(1) + newTWB(4) + newTWB(8))*reducePercent;
        reduceAllAmt = 1- (reducedTWBamt/sum(newTWB));
        reduce131415Amt = 1- (reducedTWBamt/(newTWB(node1) + newTWB(node2) + newTWB(node3)));
        reduce148Amt = 1- (reducedTWBamt/(newTWB(1) + newTWB(4) + newTWB(8)));

    end
%%
    reducedTWB_1 =  reduceAllAmt*newTWB; 
    reducedMWB_1 = reducedTWB_1./H;
    [Y, t, R0s] = schisto_model( n, N', H', loc, spaths, A, mobility, isNetwork, reducedMWB_1);
    [twb_1yr_4,twb_2yr_4,  twb_3yr_4, twb_5yr_4, twb_10yr_4, infected_1yr_4,infected_2yr_4, infected_3yr_4, infected_5yr_4, infected_10yr_4] = calculate_total_twb(n, t, Y, H, loc, A, mobility);
    mean_worm_burden = reshape(Y(:,1,:),length(t),n);
    final_twb4 = sum(mean_worm_burden(end,:).*H);

    
    reducedTWB_2 = newTWB;
    reducedTWB_2(node1) = reduce131415Amt* reducedTWB_2(node1); 
    reducedTWB_2(node2) = reduce131415Amt* reducedTWB_2(node2); 
    reducedTWB_2(node3) = reduce131415Amt* reducedTWB_2(node3); 
    reducedMWB_2 = reducedTWB_2./H;
    [Y, t, R0s] = schisto_model( n, N', H', loc, spaths, A, mobility, isNetwork, reducedMWB_2);
    [twb_1yr_5, twb_2yr_5, twb_3yr_5, twb_5yr_5, twb_10yr_5, infected_1yr_5,infected_2yr_5, infected_3yr_5, infected_5yr_5, infected_10yr_5] = calculate_total_twb(n, t, Y, H, loc, A, mobility);
    mean_worm_burden = reshape(Y(:,1,:),length(t),n);
    final_twb5 = sum(mean_worm_burden(end,:).*H);
    
    
    reducedTWB_3 = newTWB;
    reducedTWB_3(1) = reduce148Amt* reducedTWB_3(1); 
    reducedTWB_3(4) = reduce148Amt* reducedTWB_3(4); 
    reducedTWB_3(8) = reduce148Amt* reducedTWB_3(8); 
    reducedMWB_3 = reducedTWB_3./H;
    [Y, t, R0s] = schisto_model( n, N', H', loc, spaths, A, mobility, isNetwork, reducedMWB_3);
    [twb_1yr_6, twb_2yr_6, twb_3yr_6, twb_5yr_6, twb_10yr_6, infected_1yr_6,infected_2yr_6, infected_3yr_6, infected_5yr_6, infected_10yr_6] = calculate_total_twb(n, t, Y, H, loc, A, mobility);
    mean_worm_burden = reshape(Y(:,1,:),length(t),n);
    final_twb6 = sum(mean_worm_burden(end,:).*H);
    
    twb_matrix = [twb_1yr_4, twb_1yr_5, twb_1yr_6;  twb_3yr_4,  twb_3yr_5, twb_3yr_6; twb_5yr_4,  twb_5yr_4, twb_5yr_6];
    n_inf_matrix = [infected_1yr_4, infected_1yr_5, infected_1yr_6; infected_3yr_4,infected_3yr_5,infected_3yr_6; infected_5yr_4,infected_5yr_5,infected_5yr_6];
    figure; bar(twb_matrix);
    set(gca,'XTickLabel',{'1 Year', '3 Years', '5 Years'})
    ylabel('Cumulative Total Worm Burden')
    legend('Treat all', 'Treat Nodes 13,14, & 15', 'Treat Nodes 1,4, & 8');
%     figure; bar(n_inf_matrix);
%     set(gca,'XTickLabel',{'1 Year', '3 Years', '5 Years'})
%     ylabel('Cumulative Number of Infected Individuals')
%     legend('Treat all', 'Treat Nodes 13,14, & 15', 'Treat Nodes 1,4, & 8');
% 


%%
%     %% reducing by same percentage 
%     %reduce all nodes by 18%
%     reducedTWB_1 =  .82*newTWB; 
%     reducedMWB_1 = reducedTWB_1./H;
%     [Y, t, R0s] = schisto_model( n, N', H', loc, spaths, A, mobility, isNetwork, reducedMWB_1);
%     [twb_1yr_1, twb_5yr_1, twb_10yr_1] = calculate_total_twb(n, t, Y, H);
% 
%     %reduce nodes 13,14,15 by 90%
%     reducedTWB_2 = newTWB;
%     reducedTWB_2(13) = 0.1* newTWB(13); 
%     reducedTWB_2(14) = 0.1* newTWB(14); 
%     reducedTWB_2(15) = 0.1* newTWB(15); 
%     reducedMWB_2 = reducedTWB_2./H;
%     [Y, t, R0s] = schisto_model( n, N', H', loc, spaths, A, mobility, isNetwork, reducedMWB_2);
%     [twb_1yr_2, twb_5yr_2, twb_10yr_2] = calculate_total_twb(n, t, Y, H);
% 
%     %reduce nodes 1,4,8 by 90%
%     reducedTWB_3 = newTWB;
%     reducedTWB_3(1) = 0.1* newTWB(1); 
%     reducedTWB_3(4) = 0.1* newTWB(4); 
%     reducedTWB_3(8) = 0.1* newTWB(8); 
%     reducedMWB_3 = reducedTWB_3./H;
%     [Y, t, R0s] = schisto_model( n, N', H', loc, spaths, A, mobility, isNetwork, reducedMWB_3);
%     [twb_1yr_3, twb_5yr_3, twb_10yr_3] = calculate_total_twb(n, t, Y, H);

%%
%     % reduce center nodes 1 and 4 by 90%
%     reducedTWB_4 = newTWB;
%     reducedTWB_4(1) = 0.1* newTWB(1); 
%     reducedTWB_4(4) = 0.1* newTWB(4); 
%     reducedMWB_4 = reducedTWB_4./H;
%     [Y, t, R0s] = schisto_model( n, N', H', loc, spaths, A, mobility, isNetwork, reducedMWB_4);
%     [twb_1yr_4,twb_2yr_4,  twb_3yr_4, twb_5yr_4, twb_10yr_4, infected_1yr_4,infected_2yr_4, infected_3yr_4, infected_5yr_4, infected_10yr_4] = calculate_total_twb(n, t, Y, H, loc, A, mobility)
%     mean_worm_burden = reshape(Y(:,1,:),length(t),n);
%     finalmwb = mean_worm_burden(end,:);
% 
%     % reduce all nodes by an equivalent drop in TWB
%     reducedTWBamt = .9*(newTWB(1) + newTWB(4));
%     reduceAllAmt = 1- (reducedTWBamt/sum(newTWB));
%     reducedTWB_5 =  reduceAllAmt*newTWB; 
%     reducedMWB_5 = reducedTWB_5./H;
%     [Y, t, R0s] = schisto_model( n, N', H', loc, spaths, A, mobility, isNetwork, reducedMWB_5);
%     [twb_1yr_5, twb_2yr_5, twb_3yr_5, twb_5yr_5, twb_10yr_5, infected_1yr_5,infected_2yr_5, infected_3yr_5, infected_5yr_5, infected_10yr_5] = calculate_total_twb(n, t, Y, H, loc, A, mobility)
%     mean_worm_burden2 = reshape(Y(:,1,:),length(t),n);
%     finalmwb2 = mean_worm_burden2(end,:);
%     
    

    %% reduce number of snails
    
%     % reduce snails in nodes 1,  4, and 8
%     newN_1 = N;
%     newN_1(1) = newN_1(1) - 6000;
%     newN_1(4) = newN_1(4) - 6000;
%     newN_1(8) = newN_1(8) - 6000;
% 
%     [Y, t, R0s] = schisto_model( n, newN_1', H', loc, spaths, A, mobility, isNetwork, mwb);
%     [twb_1yr_6, twb_2yr_6, twb_3yr_6, twb_5yr_6, twb_10yr_6, infected_1yr_6,infected_2yr_6, infected_3yr_6 infected_5yr_6, infected_10yr_6] = calculate_total_twb(n, t, Y, H, loc, A, mobility)
%     mean_worm_burden = reshape(Y(:,1,:),length(t),n);
%     finalmwb = mean_worm_burden(end,:);
% 
%     % reduce snails  in nodes 13, 14,  and 15
%     newN_3 = N;
%     newN_3(13) = newN_3(13) - 6000;
%     newN_3(14) = newN_3(14) - 6000;
%     newN_3(15) = newN_3(15) - 6000;
%     [Y, t, R0s] = schisto_model( n, newN_3', H', loc, spaths, A, mobility, isNetwork, mwb);
%     [twb_1yr_8, twb_2yr_8, twb_3yr_8, twb_5yr_8, twb_10yr_8, infected_1yr_8,infected_2yr_8, infected_3yr_8, infected_5yr_8, infected_10yr_8] = calculate_total_twb(n, t, Y, H, loc, A, mobility)
%    
% 
%     % reduce  snails in each node
%     newN_2 = N-1200;
%     [Y, t, R0s] = schisto_model( n, newN_2', H', loc, spaths, A, mobility, isNetwork, mwb);
%     [twb_1yr_7,twb_2yr_7,twb_3yr_7, twb_5yr_7, twb_10yr_7, infected_1yr_7,infected_2yr_7,infected_3yr_7, infected_5yr_7, infected_10yr_7] = calculate_total_twb(n, t, Y, H, loc, A, mobility)
%     mean_worm_burden2 = reshape(Y(:,1,:),length(t),n);
%     finalmwb2 = mean_worm_burden2(end,:);
%     
%     twb_matrix = [twb_1yr_7, twb_1yr_8, twb_1yr_6;  twb_3yr_7,  twb_3yr_8, twb_3yr_6; twb_5yr_7,  twb_5yr_8, twb_5yr_6];
%     n_inf_matrix = [infected_1yr_7, infected_1yr_8, infected_1yr_6; infected_3yr_7,infected_3yr_8,infected_3yr_6; infected_5yr_7,infected_5yr_8,infected_5yr_6];
%     figure; bar(twb_matrix);
%     set(gca,'XTickLabel',{'1 Year', '3 Years', '5 Years'})
%     ylabel('Cumulative Total Worm Burden')
%     legend('Treat all', 'Treat Nodes 13,14, & 15', 'Treat Nodes 1,4, & 8');
%     figure; bar(n_inf_matrix);
%     set(gca,'XTickLabel',{'1 Year', '3 Years', '5 Years'})
%     ylabel('Cumulative Number of Infected Individuals')
%     legend('Treat all', 'Treat Nodes 13,14, & 15', 'Treat Nodes 1, 4, & 8');
% 


    %% reduce contamination theta_p
    thetaReduction = 0.1;
    redIndices = [1,4,8];
    [Y, t, R0s] = schisto_model_control( n, N', H', loc, spaths, A, mobility, isNetwork, mwb, thetaReduction, redIndices);
    [twb_1yr_8, twb_2yr_8, twb_3yr_8, twb_5yr_8, twb_10yr_8, infected_1yr_8,infected_2yr_8,infected_3yr_8, infected_5yr_8, infected_10yr_8] = calculate_total_twb(n, t, Y, H, loc, A, mobility)
    
    thetaReduction = 0.1;
    redIndices = [13,14,15];
    [Y, t, R0s] = schisto_model_control( n, N', H', loc, spaths, A, mobility, isNetwork, mwb, thetaReduction, redIndices);
    [twb_1yr_9, twb_2yr_9, twb_3yr_9, twb_5yr_9, twb_10yr_9, infected_1yr_9,infected_2yr_9,infected_3yr_9, infected_5yr_9, infected_10yr_9] = calculate_total_twb(n, t, Y, H, loc, A, mobility)
    
    thetaReduction = 0.02;
    redIndices = [1:15];
    [Y, t, R0s] = schisto_model_control( n, N', H', loc, spaths, A, mobility, isNetwork, mwb, thetaReduction, redIndices);
    [twb_1yr_10, twb_2yr_10, twb_3yr_10, twb_5yr_10, twb_10yr_10, infected_1yr_10,infected_2yr_10,infected_3yr_10, infected_5yr_10, infected_10yr_10] = calculate_total_twb(n, t, Y, H, loc, A, mobility)
    
    twb_matrix = [twb_1yr_10, twb_1yr_9, twb_1yr_8;  twb_3yr_10,  twb_3yr_9, twb_3yr_8; twb_5yr_10,  twb_5yr_9, twb_5yr_8];
    n_inf_matrix = [infected_1yr_10, infected_1yr_9, infected_1yr_8; infected_3yr_10,infected_3yr_9,infected_3yr_8; infected_5yr_10,infected_5yr_9,infected_5yr_8];
    figure; bar(twb_matrix);
    set(gca,'XTickLabel',{'1 Year', '3 Years', '5 Years'})
    ylabel('Cumulative Total Worm Burden')
    legend('Treat all', 'Treat Nodes 13,14, & 15', 'Treat Nodes 1,4, & 8');
%     figure; bar(n_inf_matrix);
%     set(gca,'XTickLabel',{'1 Year', '3 Years', '5 Years'})
%     ylabel('Cumulative Number of Infected Individuals')
%     legend('Treat all', 'Treat Nodes 13,14, & 15', 'Treat Nodes 1, 4, & 8');


    
end