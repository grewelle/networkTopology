function [ R0, G0s, Q ] = runjavier_network( n, N, H, loc, spaths, A, cities, m_max)
% This runs a javier model simulation and plots the resulting stability
% curves. It is adapted to take into account the network structure ( Sij in
% the calculation of the Q matrix depends on the shortest path along the
% network, not just the euclidean distance between the two points)

% n - number of nodes
% N - vector of snail populations at each node
% H - vector of human populations at each node
% loc - x and y coords of each node
% A - the adjacency matrix, which tells the connections
% spaths is a distance matrix which tells the shortest path between node i
% and j in the network

% NOT USED:
% cities - list of nodes that are considered cities
% m_max - maximum mobility for the nodes that are designated cities


 % CONSTANTS
javierConstants;

iter = 50;
alpha = 2.5;
gamma = 1/(3 * 365);

% Generate thetas
theta = exposure_rate(H, H_trans, alpha, theta_urb, theta_rur)
thetap = theta;

% Hydrological connectivity parameters, all set to zero
V = ones(1,n);
P = zeros(n);
sC = zeros(n);
sM = zeros(n);
lC = zeros(n,1);
lM = zeros(n,1);

% Generate Q
S = population_radius_network(n, loc, H, spaths);
Q = radiation_model(n, H, S);

% Calculate the original R0 values
R0 = calc_R0( a, b, theta, thetap, piC, piM, H, N, gamma, nu, muC, muM, V );

% Evaluate the stability for many different values of m
ms = [0:1/iter:1]'*ones(1,n);
%ms(:,cities) = m_max*ones(iter + 1,length(cities));

for i=1:size(ms,1)
    % Choose mobility rates
    m = ms(i,:);
    % Create G0
    G0 = make_G0(m, R0, a, b, piM, piC, gamma, nu, muC, muM, theta, ...
        thetap, P, sC, sM, V, N, H, Q, lC, lM, n);
    % Find dominant eigenvalue and add to our G0 vector
    G0s(i) = max(eig(G0));   
end

%plot cities and populations, draws a line between nodes if they are
%connected according to the adjacency matrix
figure; title('map of patches');
% circles(loc(:,1),loc(:,2),log10(H)/10,'color','b');
for i = 1:n
    for j = 1:n
        if A(i,j) > 0
        	hold on;
        	line([loc(i,1) loc(j,1)],[loc(i,2) loc(j,2)],'color', 'b')
        end
    end
end
labels = cellstr( num2str([1:size(loc,1)]') );  %' # labels correspond to their order
text(loc(:,1), loc(:,2), labels, 'VerticalAlignment','bottom', 'HorizontalAlignment','right')
axis equal; axis off;
gplot(A,loc)

% % plot mobility matrix
% figure; title('Mobility from source i to destination j');
% imagesc(Q); colorbar; caxis([0 1]);
% xlabel('destination j'); ylabel('souce i'); title('Mobility Matrix: Network Model');

% plot original R0s, thetas and pop sizes
figure; subplot(3,1,1); bar(log10(H)); ylabel('population');
subplot(3,1,2); bar(theta); ylabel('transmission rate');
subplot(3,1,3); bar(R0); ylabel('R_{0}'); xlabel('node');

% % draw a line at R0=1
hold on; line([.5,n+.5],[1,1],'Color','red');

% plot stability stats
figure; title('Stability with varying mobility');
plot(ms(:,1),G0s);
xlabel('mobility m'); ylabel('G0'); title('G0 as a function of Mobility: Network Model');

% draw a line at G0=1
hold on; line([0,1],[1,1],'Color','red');
hold on;

end