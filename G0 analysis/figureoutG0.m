clear all; close all;
%% Parameters
javierConstants;
n = 3;

cities = [];
m_max = 1;
iter = 50;

% Hydrological connectivity parameters, all set to zero
V = ones(1,n);
P = zeros(n);
sC = zeros(n);
sM = zeros(n);
lC = zeros(n,1);
lM = zeros(n,1);

m = 0.5*ones(1,n);
V = ones(1,n); %water volume, set to 1 because it isn't considered in this model

%finalN = linspace(200,2000, 10);
%finalH = linspace(100000,10000,10);
% finalN = [570 550 600 700 1050];
% finalH = [40000 20000 15000 10000 5000];
% colors = ['b'; 'g'; 'm'; 'r'; 'c'];
% for i = 1:length(finalN)
%     N =  2000*ones(1,n-1);%[2000 finalN(i) 2000]
%     N = [finalN(i) N];
%     H = 10000*ones(1,n-1);
%     H = [finalH(i) H];
%     theta = exposure_rate(H, H_trans, alpha, theta_urb, theta_rur); % calculate exposure rate as a function of population
%     thetap = theta; % we assume contamination rate = exposure rate, which isn't necessarily true
% 
%     R0 = calc_R0(a,b,theta, thetap, piC, piM, H, N, gamma, nu, muC, muM, V)
%     GUIfilename = '/Users/emilyalsentzer/Documents/DeLeo/SchistoGUI/edgelist.txt';
%     [loc, spaths, A, n] = create_GUI_network(GUIfilename, true); 
%     
%     % Generate Q
%     S = population_radius_network(n, loc, H, spaths);
%     Q = radiation_model(n, H, S);
%     
%     % Evaluate the stability for many different values of m
%     ms = [0:1/iter:1]'*ones(1,n);
%     ms(:,cities) = m_max*ones(iter + 1,length(cities));
%     %ms(:,cities) = [0:m_max/iter:m_max]'*ones(1,length(cities))
% 
%     for j=1:size(ms,1)
%         % Choose mobility rates
%         m = ms(j,:);
% 
%         % Create G0
%         G0 = make_G0(m, R0, a, b, piM, piC, gamma, nu, muC, muM, theta, ...
%             thetap, P, sC, sM, V, N, H, Q, lC, lM, n);
% 
%         % Find dominant eigenvalue and add to our G0 vector
%         G0s(j) = max(eig(G0));
% 
%     end
%     plot(ms(:,1),G0s, 'Color', colors(i));
%     xlabel('mobility m'); ylabel('G0'); title('G0 as a function of Mobility: Network Model');
%     title('Stability with varying H and N');
%     hold on;
%end
%%
N = [2000 2000 570];
H = [10000 10000 40000];
theta = exposure_rate(H, H_trans, alpha, theta_urb, theta_rur); % calculate exposure rate as a function of population
thetap = theta; % we assume contamination rate = exposure rate, which isn't necessarily true

test = calc_R0(a,b,theta, thetap, piC, piM, H, N, gamma, nu, muC, muM, V)
%% Diagonalize everything

m = 0.5*ones(1,n);

I = eye(n);
m = diag(m);
R0 = diag(R0);
H = diag(H);
N = diag(N);
V = diag(V);
theta = diag(theta);
thetap = diag(thetap);           

%% Calculate G0
A = a*b*piC*piM/(2*gamma*nu*muC*muM);
term1 = ((I-m)^2)*R0;

% Second term of the sum
term2 = N*thetap * ((I-m)*m*H*Q + Q'*H*m*(I-m) + Q'*m*m*H*Q) * theta;
% term3 = N*thetap *theta* ((I-m)*m*H*Q + Q'*H*m*(I-m) + Q'*m*m*H*Q);
% term4 = N*thetap *theta* m*((I-m)*H*Q + Q'*H*(I-m) + Q'*m*H*Q);
% term5 = N*thetap *theta* (m*(H*Q + Q'*H)+ m^2*(-H*Q -Q'*H + Q'*H*Q));
% term6 = N*thetap *theta* (m*(H*Q + Q'*H)+ m^2*(Q'*H*Q -H*Q-Q'*H));
% term7 = (I -2*m + m^2)*H + (m*(H*Q + Q'*H)+ m^2*(Q'*H*Q -H*Q-Q'*H));
term8 = m*(H*Q + Q'*H - 2*H)+ m^2*(H + Q'*H*Q - H*Q - Q'*H); %YES

% when H is all of the same: 
%term9 = 2*m*H*(Q - I)+ m^2*H*(I + Q'*Q - 2*Q);


G0 = term1 + A*term2;
% G0_2 = term1 + A*term3;
% G0_3 = term1 + A*term4;
% G0_4 = term1 + A*term5;
% G0_5 = term1 + A*term6;
% G0_6 = (A*N*thetap *theta) * term7;
% G0_7 = (A*N*thetap *theta) * (H + term8);
G0_8 = R0 + (A*N*thetap *theta)*term8 %most up to date
%%
mvals = linspace(0,1,1001);
all_Go = zeros(n,n,length(mvals));
BL_eigs = zeros(length(mvals));
lambdas = zeros(length(mvals),1);
for i = 1:length(mvals)
    m = diag(mvals(i)*ones(1,n));
    Go = R0 + (A*N*thetap *theta)* (m*(H*Q + Q'*H - 2*H)+ m^2*(H + Q'*H*Q - H*Q - Q'*H));
    BL = a*piM/2 *thetap*(I-m +Q'*m)*H*(I-m + m*Q) * theta;
    BL_eigs(i) = max(eig(BL));
    all_Go(:,:,i) = Go;
    lambdas(i) = max(eig(Go));
    plot(mvals(i), max(eig(Go)), '-');
    hold on;

end
hold on; line([0,1],[1,1],'Color','red'); hold off

% [minG0, loc] = min(lambdas);
% minM = mvals(loc);


% max(eig(G0))
eigG0 = max(eig(G0_8))

%%
%BL = a*piM/2 *thetap*(I-m +Q'*m)*H*(I-m + m*Q) * theta
%eigBL = max(eig(BL))

