function [W, H] = matrixFact(N, X, Y, folds, perm)
    [TrainN, TestN] = cvpart(N, folds, perm);
    [TrainX, TestX] = cvpart(X, folds, perm);
    [a, b] = size(X)
    lambda = logspace(-3,-1,3); %good params
    lambdaTest_ratio = zeros(folds,length(lambda));
    DCG_SIZE = 400;
    
    for i = 1:folds
        XTr = squeeze(TrainX(i, :, :));
        XTs = squeeze(TestX(i, :, :));
        NTr = squeeze(TrainN(i, :, :));
        NTs = squeeze(TestN(i, :, :));
        [N_R, N_i, ~] = NDCG(NTs, DCG_SIZE);
        lfold = zeros(1,length(lambda));
        for g = 1:length(lambda)
            %[UU, SS, VV, ~, ~, ~] = dirtyIMC(NTr ,XTr , Y, lambda(g), lambda1);
            [Wt, Ht, time] = IMC(NTr, XTr, Y, 1188, lambda(g), 50);
            CT =  XTs*Wt'*Ht*Y';
            loss = norm(CT-NTs).^2;
            %if ~any(CT) continue; end
%             [CT_R, ~, ~] = NDCG(CT, DCG_SIZE, N_i)
           % lambdaTest_ratio(i,g) = mean(CT_R./N_R);
            %lambdaTest_ratio(i,g) = loss;
            lfold(g) = loss;
            fprintf("PART %d on going.\n", i);
        end
        lambdaTest_ratio(i,:) = lfold;
        fprintf("PART %d Completed.\n", i);
    end
    average_lambda_ratio_test = zeros(1, length(lambda));
    for ind = 1:length(lambda)
        average_lambda_ratio_test(ind) = mean(lambdaTest_ratio(:,ind));
    end
    [best_ratioTest, best_indexTest] = max(average_lambda_ratio_test);
    fprintf("\n\nBest Test Ratio: %f with lambda %f\n", best_ratioTest, lambda(best_indexTest));
%     [UU, SS, VV, ~, ~, ~,] = dirtyIMC(N, X, Y, lambda(best_indexTest), lambda1);
%     M = UU*SS*VV';
    [W, H, time] = IMC(N, X, Y, 1188, lambda(best_indexTest), 50);

%     CC = X*M*Y';
%     disp(CC(1,1))
%     mmwrite(strcat("CC-",side,".mm.mtx"),CC);
%     mmwrite(strcat("U-",side,".mm.mtx"),UU);
%     mmwrite(strcat("V-",side,".mm.mtx"),VV);
%     mmwrite(strcat("M-",side,".mm.mtx"),M);
%     mmwrite(strcat("XR-",side,".mm.mtx"),XR);
%     mmwrite(strcat("YR-",side,".mm.mtx"),YR);
%     mmwrite(strcat("Y-",side,".mm.mtx"),Y);
%     mmwrite(strcat("X-",side,".mm.mtx"),X);
%     mmwrite(strcat("S-",side,".mm.mtx"),SS);
end