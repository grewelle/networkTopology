function [twb_1yr_1, twb_1yr_2, twb_1yr_3] = controlSimulation(n, mwb, H, N, loc, spaths, A, mobility, isNetwork)

    newstart_mwb = mwb;
    newstart_mwb(newstart_mwb < 0) = 0; %testing!
    newTWB = newstart_mwb.*H;

    %% reducing by same TWB in each node
    reducedTWBamt = newTWB(13) + newTWB(14) + newTWB(15);
    reduceAllAmt = 1- (reducedTWBamt/sum(newTWB));
    reduce148Amt = 1- (reducedTWBamt/(newTWB(1) + newTWB(4) + newTWB(8)));
    if(reduce148Amt <0) 
        disp('error')
    end

    reducedTWB_1 =  reduceAllAmt*newTWB; %reduce all nodes by 18%
    reducedMWB_1 = reducedTWB_1./H;
    [Y, t, R0s] = schisto_model( n, N', H', loc, spaths, A, mobility, isNetwork, reducedMWB_1);
    [twb_1yr_1, twb_5yr_1, twb_10yr_1] = calculate_total_twb(n, t, Y, H);

    reducedTWB_2 = newTWB;
    reducedTWB_2(13) = 0* reducedTWB_2(13); 
    reducedTWB_2(14) = 0* reducedTWB_2(14); 
    reducedTWB_2(15) = 0* reducedTWB_2(15); 
    reducedMWB_2 = reducedTWB_2./H;
    [Y, t, R0s] = schisto_model( n, N', H', loc, spaths, A, mobility, isNetwork, reducedMWB_2);
    [twb_1yr_2, twb_5yr_2, twb_10yr_2] = calculate_total_twb(n, t, Y, H);

    reducedTWB_3 = newTWB;
    reducedTWB_3(1) = reduce148Amt* reducedTWB_3(1); 
    reducedTWB_3(4) = reduce148Amt* reducedTWB_3(4); 
    reducedTWB_3(8) = reduce148Amt* reducedTWB_3(8); 
    reducedMWB_3 = reducedTWB_3./H;
    [Y, t, R0s] = schisto_model( n, N', H', loc, spaths, A, mobility, isNetwork, reducedMWB_3);
    [twb_1yr_3, twb_5yr_3, twb_10yr_3] = calculate_total_twb(n, t, Y, H);

end