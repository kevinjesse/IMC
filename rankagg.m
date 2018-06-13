function [best_ranks] = rankagg(varargin)
%     celldisp(varargin)
%     pause
    [a,b] = size(varargin{1});
    best_ranks = zeros(a,b);
    narg = nargin;
    varg = zeros(length(varargin), a, b);
    for i = 1:length(varargin)
        varg(i,:,:) = varargin{i};
    end
   
    for i = 1:a
        i
        agg = zeros(narg, b);
        for r =1:narg
            agg(r, :) = vertcat(varg(r,i, :));
        end
        size(agg)
        best_rank = zeros(1,b);
        p = 1;
        
        %this can be b (1188) or a DCG depending on length of rank used for
        %performance comparison
        while p < 11
             winner = election(agg, 'Schulze');
%             j

%              size(winner)
            l = length(winner);
%             winner(randperm(l))
            best_rank(p:p+l-1) = winner(randperm(l));

            for v = 1:l
                agg(agg == winner(v)) = [];
            end
            p=p+l;
            %p
        end
        best_ranks(i, :) = best_rank;
    end

end
