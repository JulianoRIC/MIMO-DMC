function [ FTnew ] = reshapeMatrix(G, FT, i)

%redimensionando linhas
maxL = max(G(i,:));
FTold =  zeros((maxL - size(FT,1)), size(FT,2));
FTnew =  vertcat(FT, FTold);

 end

