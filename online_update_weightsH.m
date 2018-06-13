function [g2] = online_update_weightsH(X, Y, R,W, H, l1, lr)
% theta1 is U
% theta2 is V
% Linear Regression
% grad1 = 0;
% grad2 = 0;

%lr = .0000001; iX %this is good
% lr = .00001;
% lr = .000001;
% 
% h = X*U*V'*Y';
h = X * W' * H * Y'; 
[m,n] = size(h);



%correct
grad2 =  2*W*X'*(h-R)*Y + l1*H; 
%grad2 =  2*U*X'*(h-R)*Y;




% g1 = U - (grad1*lr);
g2 = H - (grad2*lr);
% g1 = grad1;
% g2 = grad2;
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

end

