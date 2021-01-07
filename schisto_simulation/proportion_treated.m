% Calculate the proportion to treat at each node based on constant,
% proportional, betweenness centrality, or susceptible size strategies.

function [Pt, centrality] = proportion_treated(treat_strat, n, W0, H, p, plim, tnpt, A)
    switch treat_strat
        case 'none'             % Baseline, no treatment
            Pt = zeros(n,1);
            centrality = zeros(n,1);
        case 'const'            % Constant proportion treated everywhere
            Pt = p*ones(n,1);
            centrality = zeros(n,1);
        case 'prop'             % Treat proportional to expected mean worm burden
            mwb_eq = W0';
            nmwb = sum(mwb_eq.*H)/sum(H);
            scale_factors = mwb_eq/nmwb;
            Pt = p*ones(n,1).*scale_factors;
            Pt(Pt>plim) = plim;
            if(sum(H.*Pt)>tnpt)
                Pt = Pt*tnpt/sum(H.*Pt);
            end
            centrality = zeros(n,1);
        case 'bc'
            if issymmetric(A)
                G = graph(A);
            else
                G = digraph(A);
            end
            bc = betweenness_centrality_final(G,n);
            Pt = zeros(n,1);
            treat = tnpt;
            for i = 1:n
                curr_node = bc(i);
                % Include below code if skipping treatment for high-pop nodes
                % if H(curr_node) > 10000
                %     continue
                % end
                Pt(curr_node) = plim*(plim*H(curr_node)<=treat) + treat/H(curr_node)*(plim*H(curr_node)>treat);
                treat = treat - Pt(curr_node)*H(curr_node);
                if treat == 0
                    break
                end
            end
            centrality = bc;
        case 'ss'
            if issymmetric(A)
                G = graph(A);
            else
                G = digraph(A);
            end
            ss = betweenness_centrality_final(G,n);
            r = susceptible_size(G,ss,n);
            for i = 1:1000
                rows = datasample(linspace(1,n,n),2,'Replace',false);
                ss([rows(1) rows(2)]) = ss([rows(2) rows(1)]);
                r_temp = susceptible_size(G,ss,n);
                if r_temp <= r
                    r = r_temp;
                else 
                    ss([rows(1) rows(2)]) = ss([rows(2) rows(1)]);
                end
            end
            Pt = zeros(n,1);
            treat = tnpt;
            for i = 1:n
                curr_node = ss(i);
                Pt(curr_node) = plim*(plim*H(curr_node)<=treat) + treat/H(curr_node)*(plim*H(curr_node)>treat);
                treat = treat - Pt(curr_node)*H(curr_node);
                if treat == 0
                    break
                end
            end
            centrality = ss;
    end
end