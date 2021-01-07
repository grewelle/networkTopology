% This master script is used to run the javier model, schisto simulation,
% optimization, and analysis of G0

clear all; close all;
%%

% load constants
javierConstants
    
n = 3; % number of nodes (may be overwritten, depending on network topology
N = [900 1700 1700]; % snail populations
H = 10000 * ones(1,n);
cities = [];
m_max = 1;

dist = 100;
loc = [0,0; 100,0; 200,0];
D = squareform(pdist(loc));
A = D == dist;
A = A.*dist;
sparseA = sparse(A);
spaths = calc_shortest_paths(n, sparseA);


% calculate mobility vs G0 graphs for networks with explicit road topologies
[ R0, G0s, Q2 ] = runjavier_network( n, N, H, loc, spaths, A, cities, m_max); 
R0
