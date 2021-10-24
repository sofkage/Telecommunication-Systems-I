function[X] = bits_to_4PAM(b)
    j=1;
    for i=1:length(b)/2
        if(b(i)==0 && b(i+1)==0)
            X(j)=3;
        elseif (b(i)==0 && b(i+1)==1)
            X(j)=1;
        elseif (b(i)==1 && b(i+1)==1)
            X(j)= -1;
        else
            X(j)= -3;
        end
        j=j+1;
    end
end