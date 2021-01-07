%% Modeling spread of schistosomiasis
function [y_shaped, t, treated_stats, network_mwb] = schisto_model_revised(treat_human, treat_snail, treat_san, n, N, H, loc, spaths, A, m, isNetwork, W0, Y0, C0, M0, k, thresh, show_plots)
    % Load the constants according to Javier
    javierConstants;
    
    % Choose parameters
    V = ones(n,1);              % water volume
    alpha = 2.5;
    gamma = 1/(3 * 365);
    
    % Generate thetas
    theta = exposure_rate(H, H_trans, alpha, theta_urb, theta_rur);
    thetap = theta;
    
    % Choose treatment levels
    p = 0.05;                   % Default: 0.1
    r = 0.05;
    s = 0.5;
    plim = 0.7;
    rlim = 0.7;
    slim = 0.7;
    tnpt = p*sum(H);            % Total people treated
    tnrt = r*sum(N);            % Total snails removed
    tnst = s*sum(thetap);       % Total "contamination" removed
    
    [Pt, P_cent] = proportion_treated(treat_human, n, W0, H, p, plim, tnpt, A);
    [Rt, R_cent] = proportion_treated(treat_snail, n, W0, H, r, rlim, tnrt, A);
    [St, S_cent] = proportion_treated(treat_san, n, W0, H, s, slim, tnst, A);

    % Generate Q
    if isNetwork
        S = population_radius_network(n, loc, H, spaths);
     else
         S = population_radius(n, loc, H);
    end
    Q = radiation_model(n, H, S);

    R0s = calc_R0( a, b, theta, thetap, piC, piM, H, N, gamma, nu, muC, muM, V );

%%

    tspan = 365*linspace(0,100,101)';               
     y0_all_nodes = [W0; Y0; C0; M0];
     y0_all_nodes =reshape(y0_all_nodes, 4*n, 1);

    %%
    [t,y] = ode15s(@(t,y) schisto_eq_revised(t, y, a,b, m, gamma, theta, thetap, nu, muC, muM, Q, piC, piM, V, H, N, n, Pt, Rt, St), tspan, y0_all_nodes);
    y_shaped = reshape(y, length(t), 4, n); % (t,y,n)
    
    [~, P_cent_ind] = sort(P_cent);
    [~, R_cent_ind] = sort(R_cent);
    treated_stats = [H R0s (Pt > 0) P_cent_ind (Rt > 0) R_cent_ind];

    %% plot results
if(show_plots)
%     figure;
%     subplot(2,2,1)
%     plot(t,reshape(y_shaped(:,1,:),length(t),n));
%     hold on; line([365*1,365*1],[0,20],'Color','red');
%     hold on; line([365*5,365*5],[0,20],'Color','blue');
%     hold on; line([365*10,365*10],[0,20],'Color','green');
% 
%     title4 = sprintf('m = %d', m);
%     title(title4);
%     xlabel('Time');
%     ylabel('Mean Worm Burden');
%     %legend('Node 1', 'Node 2', 'Node 3', 'Node 4', 'Node 5');
% 
%     subplot(2,2,2)
%     plot(t,reshape(y_shaped(:,2,:),length(t),n));
%     title('Model of a Schistosomiasis Outbreak');
%     xlabel('Time');
%     ylabel('Prevalence of Infected Snails');
%     % legend('Node 1', 'Node 2', 'Node 3', 'Node 4', 'Node 5');
% 
%     subplot(2,2,3)
%     plot(t,reshape(y_shaped(:,3,:),length(t),n));
%     title('Model of a Schistosomiasis Outbreak');
%     xlabel('Time');
%     ylabel('Cercarial Density');
%     % legend('Node 1', 'Node 2', 'Node 3', 'Node 4', 'Node 5');
% 
%     subplot(2,2,4)
%     plot(t,reshape(y_shaped(:,4,:),length(t),n));
%     title('Model of a Schistosomiasis Outbreak');
%     xlabel('Time');
%     ylabel('Miricidial Density');
%     %legend('Node 1', 'Node 2', 'Node 3', 'Node 4', 'Node 5');
    
    %%
    mwb = reshape(y_shaped(:,1,:),length(t),n)';
    popmat = repmat(H.*(1-Pt),1,length(t));
    network_mwb = (sum(mwb.*popmat)/sum(H));
    nhinf = n_heavily_infected_revised(H,mwb,k,thresh);
    
%     scrsz = get(groot, 'ScreenSize');
%     figure('position', [(scrsz(3)-800)/2 (scrsz(4)-250)/2 800 250]);
%     subplot(1,3,1)
%     plot(t/365,network_mwb);
%     xlabel('Time (years)');
%     ylabel('Network mean worm burden');
%     axis([0 80 0 4]);
%     
%     subplot(1,3,2)
%     plot(t/365,max(mwb)');
%     xlabel('Time (years)');
%     ylabel('Maximum local mean worm burden');
%     axis([0 80 0 16]);
%     
%     subplot(1,3,3)
%     plot(t/365,nhinf);
%     xlabel('Time (years)');
%     ylabel('Total heavily infected');
%     axis([0 80 0 5000]);
    
end
    
end
% legend('Mean Worm Burden', 'Prevalence of Infected Snails', 'Cercarial Density1', 'Miricidial Density1', 'Mean Worm Burden 2', 'Prevalence of Infected Snails 2', 'Cercarial Density 2', 'Miricidial Density 2');