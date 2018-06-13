function [rate, rank, Um, Vm] = online_train_batch(R, U, V, X, Y, folds, perm)
%R = N; U = Ug; V = Vg; X = Xg; Y = Yg;
[TrainR, TestR] = cvpart(R, folds, perm);
[TrainX, TestX] = cvpart(X, folds, perm);
[j,k] = size(U);

l1=.01; %mpaa
l2=.01; %mpaa
lr = .000001;
precision = .0000001;
iters = 100000;
% iters = 10;
%cost_history = zeros(iters,1);


% for i = 1:iters

Uf = zeros(folds, j, k);
Vf = zeros(folds, j, k);
final_cost = zeros(1, 10);
for f = 1:folds    
    XTr = squeeze(TrainX(f, :, :));
    XTs = squeeze(TestX(f, :, :));
    RTr = squeeze(TrainR(f, :, :));
    RTs = squeeze(TestR(f, :, :));

    i = 1;
    prev_cost = 1/precision;
    cost = Inf;
    Ut = U;
    Vt = V;
    while abs(prev_cost - cost) > precision && (i < iters)
        prev_cost = cost;
        
        
        [Utemp] = online_update_weightsU(XTr, Y, RTr, Ut, Vt, l1, lr);    
        [Vtemp] = online_update_weightsV(XTr, Y, RTr, Ut, Vt, l2, lr);    

        cost = costFunction(XTs, Y, RTs, Utemp, Vtemp, l1, l2);
        %cost_history(i) = cost;
        if mod(i,10) == 0
            fprintf("iter: %d | cost: %d\n", i, cost);
        end
        if cost > prev_cost
            %final_cost(f) = prev_cost;
            break
        end
        Ut = Utemp;
        Vt = Vtemp;
        i= i+ 1;    
    end
    final_cost(f) = prev_cost;
    i
    Uf(f,:,:) = Ut;
    Vf(f,:,:) = Vt;
    
end


% tcost = zeros(1,10);
% for i = 1:folds
% 
% % Um = mean(Uf);
%     Um = squeeze(Uf(i, :, :));
% % Vm = mean(Vf);
%     Vm = squeeze(Vf(i, :, :));
%     cost = costFunction(X, Y, R, Um, Vm, l1, l2);
%     tcost(i) = cost;
% end
% 
% 
% [z, x] = min(tcost);
% disp(tcost)
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


[z, x] = min(final_cost);
disp(final_cost)
x
% Um = mean(Uf);
Um = squeeze(Uf(x, :, :));
% Vm = mean(Vf);
Vm = squeeze(Vf(x, :, :));



rate = X*Um*Vm*Y';
[~,rank] = sort(rate, 2, 'descend');
