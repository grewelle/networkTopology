function [ theta ] = exposure_rate( H, H_trans, a, theta_urb, theta_rur )
% This function returns the exposure rates based on population.
% H -- vector of node populations
% H_trans -- the transition population between urban and rural
% a -- parameter which controls the speed of transition from urban to rural
% theta_urb -- the urban transmission rate
% theta_rur -- the rural transmission rate

% log transform the populations
logH = log10(H);

% create theta
theta = theta_urb+(theta_rur-theta_urb)./(1+exp(a*(logH-log10(H_trans))));

end