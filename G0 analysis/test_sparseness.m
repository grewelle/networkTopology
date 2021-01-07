function test_sparseness(n, N, H, cities, m_max)

% CONSTANTS
javierConstants;

iter = 50;
alpha = 2.5;
gamma = 1/(3 * 365);

% Hydrological connectivity parameters, all set to zero
V = ones(1,n);
P = zeros(n);
sC = zeros(n);
sM = zeros(n);
lC = zeros(n,1);
lM = zeros(n,1);

% Generate thetas
theta = exposure_rate(H, H_trans, alpha, theta_urb, theta_rur);
thetap = theta;

% Calculate the original R0 values
R0 = calc_R0( a, b, theta, thetap, piC, piM, H, N, gamma, nu, muC, muM, V );

 x1 = 0;
 x2 = 100;
 x = x1 + (x2-x1).*rand(n,1);
 y = x1 + (x2-x1).*rand(n,1);     
 loc = [x y];
 
%sparseness = linspace(1,(n*n - n)/2 - n,4);
sparseness = [0 100 500 700 1000 1175];
%sparseness = [18 19 20 21];
length(sparseness)
%colors = ['b'; 'g'; 'c'; 'm'; 'k'; 'r'; 'b'; 'm'];%; 'r'];

for q = 1:length(sparseness)
    [loc, spaths, A, n] = vary_sparseness(n, sparseness(q),loc);

    % Generate Q
    S = population_radius_network(n, loc, H, spaths);
    Q = radiation_model(n, H, S);


    % Evaluate the stability for many different values of m
    ms = [0:1/iter:1]'*ones(1,n);
    ms(:,cities) = m_max*ones(iter + 1,length(cities));
    %ms(:,cities) = [0:m_max/iter:m_max]'*ones(1,length(cities))

    for i=1:size(ms,1)
        % Choose mobility rates
        m = ms(i,:);

        % Create G0
        G0 = make_G0(m, R0, a, b, piM, piC, gamma, nu, muC, muM, theta, ...
            thetap, P, sC, sM, V, N, H, Q, lC, lM, n);

        % Find dominant eigenvalue and add to our G0 vector
        G0s(i) = max(eig(G0));

    end

%     % plot cities and populations
%     figure; title('map of patches');
%     circles(loc(:,1),loc(:,2),log10(H)/10,'color','b');
%     for i = 1:n
%         for j = 1:n
%             if A(i,j) > 0
%                 hold on;
%                 line([loc(i,1) loc(j,1)],[loc(i,2) loc(j,2)],'color', 'b')
%             end
%         end
%     end
%     labels = cellstr( num2str([1:size(loc,1)]') );  %' # labels correspond to their order
%     text(loc(:,1), loc(:,2), labels, 'VerticalAlignment','bottom', 'HorizontalAlignment','right')
%     axis equal; axis off;
%     gplot(A,loc)

%     % plot mobility matrix
%     figure; title('Mobility from source i to destination j');
%     imagesc(Q); colorbar; caxis([0 1]);
%     xlabel('destination j'); ylabel('souce i'); title('Mobility Matrix: Network Model');

    % % plot original R0s, thetas and pop sizes
    % figure; subplot(3,1,1); bar(log10(H)); ylabel('population');
    % subplot(3,1,2); bar(theta); ylabel('transmission rate');
    % subplot(3,1,3); bar(R0); ylabel('R_{0}'); xlabel('node');
    % 
    % % draw a line at R0=1
    % hold on; line([.5,n+.5],[1,1],'Color','red');

    % plot stability stats
    %figure 1;
    hold on;
    plot(ms(:,1),G0s); %, 'Color', colors(q)
    xlabel('mobility m'); ylabel('G0'); title('G0 as a function of Mobility: Network Model');
    title('Stability with varying mobility');

    hold on;
end

max_edges = (n*(n-1))/2; %max number of edges in a graph
%sparseness refers to # edges removed
labels = strread(num2str((max_edges - sparseness)./max_edges),'%s');
legend(labels)
% draw a line at G0=1
hold on; line([0,1],[1,1],'Color','red');

end