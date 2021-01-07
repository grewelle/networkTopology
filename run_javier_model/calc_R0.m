% Calculates a vector of R0 values. The vector contains an R0 value for each
% node in the network.
function R0 = calc_R0(a,b,theta, thetap, piC, piM, H, N, gamma, nu, muC, muM, V)
    R0 = (a*b*piC*piM.*theta.*thetap.*H.*N)./(2*gamma*nu*muC*muM.*V.*V); 

end