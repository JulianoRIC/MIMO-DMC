function [ FTnew ] = reshapeCol( G, FT, i )

%redimensionando colunas
maxC = max(G(:,i));
FTnew = [FT zeros(size(FT,1), (maxC - size(FT,2)))];


end

