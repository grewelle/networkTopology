% This function loads in all the constants for a schistosomiasis model
% according to Javier's paper.

% gamma - per capita mortality rate of schistosome host
% nu - per capita mortality rate of snails
% a - probability of successful cercarial infection in humans
% b - probability of successful miracidial infection in snails
% piC - cercarial emission rate per infected snail
% piM - miracidial emission rate per worm pair
% muC - per capita mortality rate of cercaria
% muM - per capita mortality rate of miracidia

% Javier's constants
%timescale = hours
% gamma = 1/(3 * 365*24);
% nu = 1/(.1 * 365*24); % comes from sanna's PNAS paper %1/(.1 * 12);
% a = (1e-5)/24;
% b = (4e-6)/24; % comes from Sanna's PNAS paper
% piC = 100/24;
% piM = 300/24;
% muM = (3.94*100)/24; % multiplied by 100 to be consistent with Sanna's PNAS paper
% muC = (.91*100)/24; % multiplied by 100 to be consistent with Sanna's PNAS paper

%timescale = months
% gamma = 1/(3 * 12);
% nu = 1/(.1 * 12); % comes from sanna's PNAS paper %1/(.1 * 12);
% a = 1e-5 *30;
% b = 4e-6 *30; % comes from Sanna's PNAS paper
% piC = 100*30;
% piM = 300*30;
% muM = 3.94*100*30; % multiplied by 100 to be consistent with Sanna's PNAS paper
% muC = .91*100* 30; % multiplied by 100 to be consistent with Sanna's PNAS paper


% timescale = days
gamma = 1/(3 * 365);
nu = .12; % comes from sanna's PNAS paper %1/(.1 * 365);
a = 1e-5;
b = 4e-6; % comes from Sanna's PNAS paper
piC = 100;
piM = 300;
muM = 3.94*100; % multiplied by 100 to be consistent with Sanna's PNAS paper
muC = .91*100; % multiplied by 100 to be consistent with Sanna's PNAS paper


% Set baseline contact rate parameters 

% theta_urb = minimum exposure rate for huge cities
% theta_rur = maximum exposure rate for small villages
% alpha - a parameter that controls the steepness of the inverse sigmoid
% function
% H_trans = a parameter that controls the population size at which there is
% a large change in the exposure rate

theta_urb = .0001;
theta_rur = .9;
alpha = 2.5;
H_trans = 7000;  % 30000;

