function [avgN] = avgComp()
    gX = full(mmread("sparseXgenre.mm.mtx"));
    XR = full(mmread("XR-genre.mm.mtx"));
    YR = full(mmread("Y-genre.mm.mtx"));
    M = full(mmread("M-genre.mm.mtx"));

    TX = gX/XR; %QR = N
    TX(isnan(TX)) = 0;
    TX(isinf(TX)) = 0;
    Ng = TX * M * YR'; %X * M * Y'
    [~,SNg] = sort(Ng, 2, 'descend'); %sort for best 
    filename = strcat("SN-genre",".mm.mtx");
%     mmwrite(filename,SN);

% 
    mX = full(mmread("sparseXmpaa.mm.mtx"));
    XR = full(mmread("XR-mpaa.mm.mtx"));
    YR = full(mmread("Y-mpaa.mm.mtx"));
    M = full(mmread("M-mpaa.mm.mtx"));

    TX = mX/XR; %QR = N
    TX(isnan(TX)) = 0;
    TX(isinf(TX)) = 0;
    Nm = TX * M * YR'; %X * M * Y'
    [~,SNm] = sort(Nm, 2, 'descend'); %sort for best 
    filename = strcat("SN-mpaa",".mm.mtx");

    
%     aX = X(29:561);
%     XR = full(mmread("XR-actor.mm.mtx"));
%     YR = full(mmread("Y-actor.mm.mtx"));
%     M = full(mmread("M-actor.mm.mtx"));
% 
%     TX = aX/XR; %QR = N
%     TX(isnan(TX)) = 0;
%     TX(isinf(TX)) = 0;
%     N = TX * M * YR'; %X * M * Y'
%     [~,SNa] = sort(N, 2, 'descend'); %sort for best 
%     filename = strcat("SN-actor",".mm.mtx");
      z = cat(3,Ng,Nm);
      
      avgN = mean(z,3);
%     avg = zeros(1,length(SNg));
%       for i = 1:length(SNg)
%           gi = find(SNg==i);
%           mi = find(SNm==i);
%         %ai = find(SNa==i);
%         %avg(i) = mean([gi,mi, ai]);
%           avg(i) = mean([gi,mi]);
%       end
    
end

