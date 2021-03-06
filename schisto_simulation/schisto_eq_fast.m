%% schistosomiasis equations - Javier
function yp = schisto_eq_fast(t, y, a,b, m, gamma, theta, thetap, nu, muC, muM, Q, piC, piM, V, H, N, n)
% yp(:,1) = dW/dt; yp(:,2) = dY/dt; yp(:,3) = dC/dt; yp(:,4) = dM/dt
% y(:,1) = W; y(:,2) = Y; y(:,3) = C; y(:,4) = M
    y = reshape(y, [4,n])';
    yp_matrix = zeros(n,4); %preallocate matrix for yp
    yp_matrix(:,1) = a*((1-m).*theta.*y(:,3) + m.*(Q*(theta.*y(:,3)))) - gamma*y(:,1);
    yp_matrix(:,2) = b*y(:,4).*(1-y(:,2)) - nu*y(:,2);
    yp_matrix(:,3) = piC./V.*N.* y(:,2) - muC*y(:,3);
    yp_matrix(:,4) = piM./V.* thetap.* [(1-m).*H.*y(:,1)/2 + Q*(m.*H.*y(:,1)/2)] - muM*y(:,4);
    yp = reshape(yp_matrix', [4*n,1]); 
    fprintf('Working at t = %g.\n', t)
end