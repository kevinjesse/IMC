function [Train, Test] = cvpart(N, folds, perm)
    start = 1;
    N = N(perm,:);
    [a, b] = size(N);
    partsize = a/folds;
    Train = zeros(folds,(folds-1)*partsize,b);
    Test = zeros(folds, partsize, b);
    for part = start:folds
        if part == start
            train = N(start*partsize+1:folds*partsize,:);
            test = N(start:start*partsize,:);
        elseif part == folds
            train = N(start:(folds-1)*partsize, :);
            test = N((folds-1)*partsize+1:folds*partsize,:);
        else
            train = vertcat(N(start:(part-1)*partsize, :), ...
                N((part*partsize)+1:folds*partsize,:));
            test = N((part-1)*partsize+1:(part)*partsize, :);
        end
        Train(part,:,:) = train;
        Test(part,:,:) = test;
    end
end