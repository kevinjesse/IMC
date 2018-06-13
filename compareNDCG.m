function [r] = compareNDCG(mat)    

    Obs = mmread("sparseN.mm.mtx");
    [m,n] = size(Obs);
    obsf = Obs';
    
    DCG_SIZE = 10; %m
    DCGV= [0,-2,-1,1,4,5];
    
    [~,ci] = sort(obsf, 2, 'descend');
    Obs_i = zeros(n, m);
    Obs_rel = zeros(n, m);
    
    for i = 1:n
        
        
%         for ii = 1:DCG_SIZE
%             Obs_i(i, ci(i, ii)) = ii;
%             Obs_rel(i, ci(i, ii)) = DCGV(obsf(i,ci(i, ii))+1);
%         end

        %%OLD METHOD (works)
        rel = DCG_SIZE;
        for ii = 1:DCG_SIZE
            Obs_i(i, ci(i, ii)) = ii;
            Obs_rel(i, ci(i, ii)) = rel;
            rel = rel - 1;
        end
    end
    
   ObsDCG = dcg(Obs_rel,Obs_i); 
   

    [~,ci] = sort(mat, 2, 'descend');
    mat_i = zeros(n, m);
    mat_rel = zeros(n, m);
    
    for i = 1:n
%         for ii = 1:DCG_SIZE
%             mat_i(i, ci(i, ii)) = ii;
%             mat_rel(i, ci(i, ii)) = DCGV(round(mat(i,ci(i, ii)))+1);
%         end
        rel = DCG_SIZE;
        for ii = 1:DCG_SIZE
            mat_i(i, ci(i, ii)) = ii;
            mat_rel(i, ci(i, ii)) = rel;
            rel = rel - 1;
        end
    end
   
   matDCG = dcg(mat_rel,Obs_i); 
   
   r = mean(matDCG ./ ObsDCG);
   
   
   