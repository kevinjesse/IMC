function [rate, rank, Um, Vm] = online_test(R, U, V, X, Y, folds, perm)
%R = N; U = Ug; V = Vg; X = Xg; Y = Yg;
% [TrainR, TestR] = cvpart(R, folds, perm);
% [TrainX, TestX] = cvpart(X, folds, perm);
[p,o] = size(R);
[j,k] = size(U);


l1=.01; %mpaa
l2=.01; %mpaa
lr = .000001;
% iters = 10;
%cost_history = zeros(iters,1);


% for i = 1:iters

% Uf = zeros(folds, j, k);
% Vf = zeros(folds, j, k);
cost_history = zeros(1, o);
Ut = U;
Vt = V;

[~,ideal] = sort(R, 2, 'descend');
for  i = 1:o
    [Utemp] = online_update_weightsU(X, Y(1:i, :), R(:,1:i), Ut, Vt, l1, lr);    
    [Vtemp] = online_update_weightsV(X, Y(1:i, :), R(:,1:i), Ut, Vt, l2, lr);    



    N=X*Ut*Vt*Y';
    [~,other] = sort(N, 2, 'descend');
%     disp(size(ideal))
%     disp(size(other))
    %cost = costFunction(X, Y, R, Utemp, Vtemp, l1, l2);
    cost = perfComp_test(ideal, other, i);
    
    
    cost_history(i) = cost;
    
    if mod(i,10) == 0
        fprintf("iter: %d | cost: %d\n", i, cost);
    end
    Ut = Utemp;
    Vt = Vtemp; 
end

figure
plot(1:1188, cost_history)
title("DCG of Online learning with Incrementally New Movie Information")
ylabel('DCG')
xlabel('Number of Movies User Has Rated')

% [z, x] = min(final_cost);
% disp(final_cost)
% x
% % Um = mean(Uf);
% Um = squeeze(Uf(x, :, :));
% % Vm = mean(Vf);
% Vm = squeeze(Vf(x, :, :));
% 
% 
% 
% rate = X*Um*Vm*Y';
% [~,rank] = sort(rate, 2, 'descend');
