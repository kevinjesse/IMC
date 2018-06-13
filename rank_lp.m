function [aggr_rank] = rank_lp(ranks)
    [n_voters, n_candidates] = size(ranks);
    edge_weights = rank_build(ranks);
    c = -1*reshape(edge_weights', [], 1)';
    pairwise_constraints = zeros(((n_candidates * (n_candidates - 1)) / 2), n_candidates.^2);
    comb = combnk(0:n_candidates-1,2);
    permu = perms(0:n_candidates-1);
    permu = flipud(permu(:,1:3));
    %permu = combinator(n_candidates,3, 'p') - 1;
    ti = comb(:,1);
    tj = comb(:,2);
    
    %Pairwise constraints
    [g,h] = size(pairwise_constraints);
    for x=1:length(ti)
        i = ti(x);
        j = tj(x);
        a = zeros(1, h); a(ind(i,j,n_candidates)) = 1; a(ind(j,i,n_candidates)) = 1;
        pairwise_constraints(x,:) = a;
    end
    
    ti = permu(:,1);
    tj = permu(:,2);
    tk = permu(:,3);
    %triangle_constraints
    triangle_constraints = zeros(((n_candidates * (n_candidates - 1) * (n_candidates -2))), n_candidates.^2);
    for x = 1:length(ti)
        i = ti(x);
        j = tj(x);
        k = tk(x);
        %disp('h')
        a = zeros(1, h); a(ind(i,j,n_candidates)) = 1; a(ind(j,k,n_candidates)) = 1; a(ind(k,i,n_candidates)) = 1;
        
        triangle_constraints(x,:) = a;
    end
    
    constraints = vertcat(pairwise_constraints, triangle_constraints);
    
    [a,b] = size(pairwise_constraints);
    [d,e] = size(triangle_constraints);
    constraints_rhs = [ones(1,a),ones(1,d)];
    constraints_signs = [zeros(1,a),ones(1,d)];
    xint = 1:n_candidates^2;
    [obj,x,duals] = lp_solve(c, constraints, constraints_rhs, constraints_signs, [],[],xint);
    rankmat = reshape(x, n_candidates, n_candidates);
    aggr_rank = sum(rankmat)+1
end




