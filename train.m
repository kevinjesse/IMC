%load variables
format long
%load data and perturb the data
TmpXg = full(mmread("sparseXgenre.mm.mtx"));
TmpYg = full(mmread("sparseYgenre.mm.mtx"));
TmpXm = full(mmread("sparseXmpaa.mm.mtx"));
TmpYm = full(mmread("sparseYmpaa.mm.mtx"));
TmpXa = full(mmread("sparseXactor.mm.mtx"));
TmpYa = full(mmread("sparseYactor.mm.mtx"));
N = mmread("sparseN.mm.mtx");
N = N';
%N = full(Obs')./5;
[~,SN] = sort(N, 2, 'descend');

[Xg, Yg] = perturb(TmpXg, TmpYg);
[Xm, Ym] = perturb(TmpXm, TmpYm);
[Xa, Ya] = perturb(TmpXa, TmpYa);

Xfeature = zeros(5000,15028);
Yfeature = zeros(1188,15028);
for i = 1:5000
    gm = Xg(i,:)'*Xm(i,:);
    ga = Xg(i,:)'*Xa(i,:);
    ma = Xm(i,:)'*Xa(i,:);
    gm = gm(:);
    ga = ga(:);
    ma = ma(:);
    Xfeature(i,:) = [gm' ga' ma'];
end
for i = 1:1188
    gm = Yg(i,:)'*Ym(i,:);
    ga = Yg(i,:)'*Ya(i,:);
    ma = Ym(i,:)'*Ya(i,:);
    gm = gm(:);
    ga = ga(:);
    ma = ma(:);
    Yfeature(i,:) = [gm' ga' ma'];
end

% X = [Xg Xm Xa]./3;
% Y = [Yg Ym Ya]./3;
%[Xgm, Ygm, ~] = perturb([TmpXg TmpXm]./2, [TmpYg TmpYm]./2);
[a,~] = size(Xg);

%create random perm for data that will be used for all operations
folds = 10;
perm = randperm(a);

%matrix factorization
[Wg, Hg] = matrixFact(N, Xfeature, Yfeature, folds, perm);
return
%[Wm, Hm] = matrixFact(N, Xm, Ym, folds, perm);
%[Ua, Va, Ma] = matrixFact(N, Xa, Ya, folds, perm);
%[Ugm, Vgm, Mgm] = matrixFact(N, Xgm, Ygm, folds, perm);


%tempory save
% mmwrite(strcat("U-genre.mm.mtx"),Ug); mmwrite(strcat("U-mpaa.mm.mtx"),Um);
% mmwrite(strcat("V-genre.mm.mtx"),Vg); mmwrite(strcat("V-mpaa.mm.mtx"),Vm);
% mmwrite(strcat("M-genre.mm.mtx"),Mg); mmwrite(strcat("M-mpaa.mm.mtx"),Mm);
% mmwrite(strcat("U-genre-mpaa.mm.mtx"),Ugm);
% mmwrite(strcat("V-genre-mpaa.mm.mtx"),Vgm);
% mmwrite(strcat("M-genre-mpaa.mm.mtx"),Mgm); 
% mmwrite(strcat("XR-genre.mm.mtx"),XRg);
% mmwrite(strcat("XR-mpaa.mm.mtx"),XRm); 

%[rateg, rankg, Ug_on, Vg_on] = online_train_batch(N, Ug, Vg, Xg, Yg, folds, perm);
%[rateg, rankg, Ug_on, Vg_on] = online_train(N, Ug, Vg, Xg, Yg);
%[ratem, rankm, Um_on, Vm_on] = online_train_batch(N, Um, Vm, Xm, Ym, folds, perm);
%[ratem, rankm, Um_on, Vm_on] = online_train(N, Um, Vm, Xm, Ym);


% Ng = Xg * Ug * Vg' * Yg'; %X * M * Y'
% [~,SNg] = sort(Ng, 2, 'descend');
% [~,SNm] = sort(Nm, 2, 'descend');
Nm = Xm * Wm * Hm' * Ym'; %X * M * Y'
[~,SNm] = sort(Nm, 2, 'descend');
temp = perfComp(SN, SNm)
% 
% Ng_on = Xg * Ug_on * Vg_on' * Yg'; %X * M * Y'
% [~,SNg_on] = sort(Ng_on, 2, 'descend');
% 
% Nm_on = Xm * Um_on * Vm_on' * Ym'; %X * M * Y'
% [~,SNm_on] = sort(Nm_on, 2, 'descend');

%ranksSN = rankagg(SNg(1,:), SNm(1,:));
% pGM_base = perfComp(SN(1,:), ranksSN(1,:))
% pNm_base =  perfComp(SN(1,:), SNm(1,:))
% pNg_base =  perfComp(SN(1,:), SNg(1,:))

% table = zeros(5,4);
% table(1,1) = perfComp(SN, SNg); 
% table(2,1) = perfComp(SN, SNm);
% table(4,1) = perfComp(SN, SNgm);
% 
% table(1,2) = perfComp(SN(1,:), rankagg(SNg(1,:))); 
% table(2,2) = perfComp(SN(1,:), rankagg(SNm(1,:)));
% 
% tic
% table(4,2) = perfComp(SN(1,:), rankagg(SNg(1,:), SNm(1,:)));
% toc
% 
% table(1,3) = perfComp(SN(1,:), rankagg(rankg(1,:)));
% table(2,3) = perfComp(SN(1,:), rankagg(rankm(1,:)));
% table(4,3) = perfComp(SN(1,:), rankagg(rankg(1,:), rankm(1,:)));








