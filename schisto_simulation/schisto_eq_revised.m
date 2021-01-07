%% schistosomiasis equations - Javier
function yp = schisto_eq_revised(t, y, a,b, m, gamma, theta, thetap, nu, muC, muM, Q, piC, piM, V, H, N, n, Pt, Rt, St)
% yp(:,1) = dW/dt; yp(:,2) = dY/dt; yp(:,3) = dC/dt; yp(:,4) = dM/dt
% y(:,1) = W; y(:,2) = Y; y(:,3) = C; y(:,4) = M
    y = reshape(y, [4,n])';
    yp_matrix = zeros(n,4); %preallocate matrix for yp
    yp_matrix(:,1) = a*(1-Pt).*((1-m).*theta.*y(:,3) + m.*(Q*(theta.*y(:,3)))) - gamma*y(:,1);
    yp_matrix(:,2) = b*y(:,4).*(1-y(:,2)) - nu*y(:,2);
    yp_matrix(:,3) = (piC./V).*(1-Rt).*N.* y(:,2) - muC*y(:,3);
    yp_matrix(:,4) = (piM./V).* (1 - St) .* thetap.* [(1-m).*(1-Pt).*H.*(y(:,1)/2) + Q'*(m.*(1-Pt).*H.*y(:,1)/2)] - muM*y(:,4); %ADDED a TRANSPOSE IN FRONT OF Q _ MAY BE WRONG
    yp = reshape(yp_matrix', [4*n,1]); 
    %fprintf('Working at t = %g.\n', t)
end

% 
%     yp = zeros(4,n_patches); %preallocate matrix for yp
%     yp(:,1) = a*((1-m).*theta.*y(:,3) + m.*(Q*(theta.*y(:,3)))) - gamma*y(:,1);
%     yp(:,2) = b*y(:,4).*(1-y(:,2)) - nu*y(:,2);
%     yp(:,3) = pi_c./V.*N.* y(:,2) - mu_c*y(:,3);
%     yp(:,4) = pi_m./V.* theta_p.* [(1-m).*H.*y(:,1)/2 + Q*(m.*H.*y(:,1)/2)] - mu_m*y(:,4)

%   yp = zeros(4*n, 1); %preallocate matrix for yp
%     for i =1:n
%         yp(1*n) = a*((1-m(n)).*theta(n).*y(3*n) + m(n).*(Q*(theta(n).*y(n*3)))) - gamma*y(1*n);
%         yp(2*n) = b*y(n*4).*(1-y(n*2)) - nu*y(n*2);
%         yp(3*n) = pi_c./V(n).*N(n).* y(n*2) - mu_c*y(n*3);
%         yp(4*n) = pi_m./V(n).* theta_p(n).* [(1-m(n)).*H(n).*y(n*1)/2 + Q*(m(n).*H(n).*y(n*1)/2)] - mu_m*y(n*4)
%     end