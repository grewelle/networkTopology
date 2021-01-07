%sparseness: fully connected = 0, not connected = (n^n - n)/2
function [loc, spaths, A, n] = vary_sparseness(n, sparseness,loc)

    A = ones(n);
    %i = randsample([1:n],sparseness); 
    %j = randsample([1:n],sparseness); 
    
%     k = 1:sparseness;
%     q = floor(sqrt(8*(k-1) + 1)/2 + 3/2);
%     p = k - (q-1).*(q-2)/2;
%     rand_i = randperm(length(q));
%     q=q(rand_i);        
%     p=p(rand_i);
%     indices = [q' p']
%     [k;p;q]';

    
    vals = ((1:n)'*ones(1, n));           % each column has the numbers 1 -> n
    idxs1 = squareform(tril(vals', -1))'; 
    idxs2 = squareform(tril(vals, -1))';   
    all_pairs = [idxs1, idxs2];        % this contains all possible pairs
    idx_to_use = randperm( size(all_pairs, 1), sparseness );  % choosing random k pairs
    pairs = all_pairs(idx_to_use, :);
    
    
    tally = 1;
    i = 1;
    while tally < sparseness
        A(pairs(i,1),pairs(i,2)) = 0;
        A(pairs(i,2),pairs(i,1)) = 0;
        if(sum(all(A==0,2))) > 0
            A(q(i),p(i)) = 1;
            A(p(i),q(i)) = 1;
        else
            tally = tally + 1;
        end
        i = i + 1;
    end
    D = squareform(pdist(loc));
    A = A.*D;
    spaths = calc_shortest_paths(n, sparse(A));

end