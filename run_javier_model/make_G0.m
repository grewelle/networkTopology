function [ G0 ] = make_G0( m, R0, a, b, piM, piC, gamma, nu, muC, muM, ...
                    theta, thetap, P, sC, sM, V, N, H, Q, lC, lM, n)
                
% This function generates G0 given all the constants
% We will ignore all the hydrological parts

% Diagonalize everything
I = eye(n);
m = diag(m);
R0 = diag(R0);
H = diag(H);
N = diag(N);
V = diag(V);
theta = diag(theta);
thetap = diag(thetap);

% First term of the sum
term1 = ((I-m)^2)*R0;

% Second term of the sum
A = a*b*piC*piM/(2*gamma*nu*muC*muM);
term21 = N*inv(V)*thetap;
term22 = (I-m)*m*H*Q + Q'*H*m*(I-m) + Q'*m*m*H*Q;
term23 = theta*inv(V);
term2 = A*term21*term22*term23;

% Third term of the sum is all hydrologic
term3 = zeros(n);

% Add them all together
G0 = term1 + term2 + term3;

end

