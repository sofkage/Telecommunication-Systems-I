function X = bits_to_4_PAM(bit_sec, A)
j=1;
for i=1:2:length(bit_sec)-1
    if(bit_sec(i)==0 && bit_sec(i+1)==0)
        X(j)=3*A;
    elseif (bit_sec(i)==0 && bit_sec(i+1)==1)
        X(j)=A;
    elseif (bit_sec(i)==1 && bit_sec(i+1)==1)
        X(j)=-A;
    else
        X(j)=-3*A;
    end
j=j+1;
end
end

