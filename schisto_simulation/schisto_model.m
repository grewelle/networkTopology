%% Modeling spread of schistosomiasis
function [y_shaped, t, R0s] = schisto_model( n, N, H, loc, spaths, A, m, isNetwork, W0)
    % Load the constants according to Javier
    javierConstants;

    % Choose parameters
    V = ones(n,1);              % water volume
    alpha = 2.5;
    gamma = 1/(3 * 365);


    % Generate thetas
    theta = exposure_rate(H, H_trans, alpha, theta_urb, theta_rur);
    thetap = theta;

    % Generate Q
    if isNetwork
        S = population_radius_network(n, loc, H, spaths);
     else
         S = population_radius(n, loc, H);
    end
    Q = radiation_model(n, H, S);

    R0s = calc_R0( a, b, theta, thetap, piC, piM, H, N, gamma, nu, muC, muM, V );

%%

    tspan = [0 20000]';               
     %W0 = 5*ones(1,n);%1./H';%.*ones(1,n)%H'%-30/200000.*H' + 30%randr( 0, 10, 1, n);
     Y0 = 0.05*ones(1,n);%.03 randr( 0.01, .05, 1, n); %-0.05/200000.*H' + 0.05;
     C0 = 80*ones(1,n);% 80 randr( 5, 60, 1, n);
     M0 = 100*ones(1,n);% 100 randr( 25, 100, 1, n);
     %W0 = 10*ones(1,n);%1./H';%.*ones(1,n)%H'%-30/200000.*H' + 30%randr( 0, 10, 1, n);
     %Y0 = 0.05*ones(1,n);%.03 randr( 0.01, .05, 1, n); %-0.05/200000.*H' + 0.05;
     %C0 = 80*ones(1,n);% 80 randr( 5, 60, 1, n);
     %M0 = 100*ones(1,n);% 100 randr( 25, 100, 1, n);
     y0_all_nodes = [W0; Y0; C0; M0];
     y0_all_nodes =reshape(y0_all_nodes, 4*n, 1);

    %%
    [t,y] = ode15s(@(t,y) schisto_eq(t, y, a,b, m, gamma, theta, thetap, nu, muC, muM, Q, piC, piM, V, H, N, n), tspan, y0_all_nodes);
    y_shaped = reshape(y, length(t), 4, n); % (t,y,n)

    %% plot results
    figure;
    subplot(2,2,1)
    plot(t,reshape(y_shaped(:,1,:),length(t),n));
    hold on; line([365*1,365*1],[0,20],'Color','red');
    hold on; line([365*5,365*5],[0,20],'Color','blue');
    hold on; line([365*10,365*10],[0,20],'Color','green');

    title4 = sprintf('m = %d', m);
    title(title4);
    xlabel('Time');
    ylabel('Mean Worm Burden');
    %legend('Node 1', 'Node 2', 'Node 3', 'Node 4', 'Node 5');

    subplot(2,2,2)
    plot(t,reshape(y_shaped(:,2,:),length(t),n));
    title('Model of a Schistosomiasis Outbreak');
    xlabel('Time');
    ylabel('Prevalence of Infected Snails');
    % legend('Node 1', 'Node 2', 'Node 3', 'Node 4', 'Node 5');

    subplot(2,2,3)
    plot(t,reshape(y_shaped(:,3,:),length(t),n));
    title('Model of a Schistosomiasis Outbreak');
    xlabel('Time');
    ylabel('Cercarial Density');
    % legend('Node 1', 'Node 2', 'Node 3', 'Node 4', 'Node 5');

    subplot(2,2,4)
    plot(t,reshape(y_shaped(:,4,:),length(t),n));
    title('Model of a Schistosomiasis Outbreak');
    xlabel('Time');
    ylabel('Miricidial Density');
    %legend('Node 1', 'Node 2', 'Node 3', 'Node 4', 'Node 5');
end
% legend('Mean Worm Burden', 'Prevalence of Infected Snails', 'Cercarial Density1', 'Miricidial Density1', 'Mean Worm Burden 2', 'Prevalence of Infected Snails 2', 'Cercarial Density 2', 'Miricidial Density 2');