

Obs = mmread("sparseN.mm.mtx");
[m,n] = size(Obs);
obsf = Obs';
obsf = obsf./5;

DCG_SIZE = 400; %m
%DCGV= [0,-2,-1,1,4,5];

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

%without online
X = full(mmread("X-genre.mm.mtx"));
% XR = full(mmread("XR.mm.mtx"));
Y = full(mmread("Y-genre.mm.mtx"));

%     YR = full(mmread("sparseYactor.mm.mtx"));
%     YR = YR(:,1:200);
%M = full(mmread("M.mm.mtx"));
U = full(mmread("U-genre.mm.mtx"));
V = full(mmread("V-genre.mm.mtx"));

% TX = X/XR; %QR = N
% TX(isnan(TX)) = 0;
% TX(isinf(TX)) = 0;
N = X * U*V * Y'; %X * M * Y'
[~,SN] = sort(N, 2, 'descend'); %sort for best
MC_i = zeros(n, m);
MC_rel = zeros(n, m);

for i = 1:n
    rel = DCG_SIZE;
    for ii = 1:DCG_SIZE
        MC_i(i, SN(i, ii)) = ii;
        MC_rel(i, SN(i, ii)) = rel;
        rel = rel - 1;
    end
end
matDCG = dcg(MC_rel,Obs_i); 
r = mean(matDCG ./ ObsDCG);
fprintf("AVG NDCG with just MC: %f\n", r);


%temp
[rate, rank, U, V] = online_train(obsf, U, V, X,Y);

[~,SNO] = sort(rate, 2, 'descend');
MCO_i = zeros(n, m);
MCO_rel = zeros(n, m);

for i = 1:n
    rel = DCG_SIZE;
    for ii = 1:DCG_SIZE
        MCO_i(i, SNO(i, ii)) = ii;
        MCO_rel(i, SNO(i, ii)) = rel;
        rel = rel - 1;
    end
end

matDCG = dcg(MCO_rel,Obs_i); 
r = mean(matDCG ./ ObsDCG);
fprintf("AVG NDCG with Online MC: %f\n", r);
   
   