m = [0:.1:1];
for i = 1:length(m)
    a = 2e5;
    k = 1.3;
    m = m(i) * ones(n,1);
    [loc, spaths, A, n] = create_GUI_network(GUIfilename, true); 
    H = zipf_population(n, a, k)';
    S = population_radius_network(n, loc, H, spaths);
    Q = radiation_model(n, H, S);
    G0term =  m*(H*Q + Q'*H - 2*H) + m^2*(H + Q'*H*Q - H*Q - Q'*H);
end