function centralization = calc_centralization(A, n)
    degree_seq = sum(full(A),2);
    [max_degree, max_index] = max(degree_seq);
 	max_centrality = (n-1)*(n-2);
 	diff_centrality = 0;
	for i=1:length(degree_seq)
		diff = max_degree - degree_seq(i);
		diff_centrality = diff_centrality + diff;
    end
    centralization = diff_centrality./max_centrality;
end