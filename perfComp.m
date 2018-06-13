function [r] = perfComp(ideal, other)
    
    [n, m] = size(ideal);
%     [~,ci] = sort(ideal, 2, 'descend');
    Obs_i = zeros(n, m);
    Obs_rel = zeros(n, m);
    DCG_SIZE = 400;
    for i = 1:n
        rel = DCG_SIZE;
        for ii = 1:DCG_SIZE
            Obs_i(i, ideal(i, ii)) = ii;
            Obs_rel(i, ideal(i, ii)) = rel;
            rel = rel - 1;
        end
    end

    ObsDCG = dcg(Obs_rel,Obs_i); 
    
%     [~,SN] = sort(other, 2, 'descend'); %sort for best
    MC_i = zeros(n, m);
    MC_rel = zeros(n, m);
    for i = 1:n
        rel = DCG_SIZE;
        for ii = 1:DCG_SIZE
            MC_i(i, other(i, ii)) = ii;
            MC_rel(i, other(i, ii)) = rel;
            rel = rel - 1;
        end
    end
    matDCG = dcg(MC_rel,Obs_i); 
    r = mean(matDCG ./ ObsDCG);
end