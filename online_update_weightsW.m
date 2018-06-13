function [g1] = online_update_weightsW(X, Y, R, W, H, l1,lr)
% theta1 is U
% theta2 is V
% Linear Regression
% grad1 = 0;
% grad2 = 0;

%lr = .0000001; iX
% lr = .00001; %this is good
% lr = .000001;

%lr = .1;

% h = X*U*V'*Y';
% size(X)
% size(W)
% size(H)

h = X * W' * H * Y'; 
grad1 =  2*H*X'*(h-R)*Y + l1*W;
g1 = W - (grad1*lr);

end

