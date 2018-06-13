function [X, Y] = perturb(X, Y)
    X(isnan(X))=0; X(isinf(X))=0;
    Y(isnan(Y))=0; Y(isinf(Y))=0;
%     [X, XR] = qr(TmpX,0);
%     [Y, ~] = qr(TmpY,0);
end