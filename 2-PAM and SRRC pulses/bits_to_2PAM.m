function X = bits_to_2PAM(b)
  
    for i = 1:length(b)
        if (b(i)==0)
            X(i)=1;
        else
            X(i)=-1;
        end
    end
end