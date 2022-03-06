function [ c ] = numberOfEmptyRows( FT )

c=0;
for i=1:length(FT)
    if(any(FT(i,:)) == false)
        c = c+ 1;
    end
    if(any(FT(i,:)) == true)
        break
    end
end


