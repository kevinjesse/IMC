function [R, N_i, N_rel] = NDCG(N, DCG_SIZE, order)
    [a,b] = size(N);
    [~,ci] = sort(N, 2, 'descend');
    N_i = zeros(a, b);
    N_rel = zeros(a, b);
    for i = 1:a
        rel = DCG_SIZE;
        for ii = 1:DCG_SIZE
            N_i(i, ci(i, ii)) = ii;
            N_rel(i, ci(i, ii)) = rel;
            rel = rel - 1;
        end
    end

    if ~exist('order','var')
        R = dcg(N_rel,N_i);
    else
        R = dcg(N_rel,order);
    end
end