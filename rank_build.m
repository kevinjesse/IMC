function [edge_weights] = rank_build(ranks)
   [n_voters, n_candidates] = size(ranks);
   edge_weights = zeros(n_candidates,n_candidates);
   c = combnk(1:n_candidates,2);
   [ni, nj] = size(c);
   ti = c(:,1);
   tj = c(:,2);
   for ind= 1:ni
       i = ti(ind);
       j = tj(ind);
       pref = ranks(:,i) - ranks(:,j);
       h_ij = sum(pref < 0);
       h_ji = sum(pref > 0);
       if h_ij > h_ji
           edge_weights(i,j) = h_ij-h_ji;
       elseif h_ij < h_ji
           edge_weights(j,i) = h_ji-h_ij;
       end
   end
end

