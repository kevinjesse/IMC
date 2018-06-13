function [MCR, rank] = predict(Xg, Xm, Xa, W_path, H_path, Y_path)
    load(W_path);
    load(H_path);
    load(Y_path);
    gm = Xg(1,:)'*Xm(1,:);
    ga = Xg(1,:)'*Xa(1,:);
    ma = Xm(1,:)'*Xa(1,:);
    gm = gm(:);
    ga = ga(:);
    ma = ma(:);
    Xf(1,:) = [gm' ga' ma'];
    Xf(isnan(Xf)) = 0; Xf(isinf(Xf)) = 0;

    disp(size(W'));
    disp(size(H));
    disp(size(Y'));


    MCR = Xf * W' * H * Y'; %X * M * Y'
    [~,rank] = sort(MCR, 2, 'descend'); %sort for best  
end

