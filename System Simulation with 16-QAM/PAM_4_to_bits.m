function [est_bit] = PAM_4_to_bits(X,A)
j=1;
for i=1:length(X)
    if(X(i)==3*A)
        est_bit(j)=0;
        est_bit(j+1)=0;
    elseif (X(i)==A)
        est_bit(j)=0;
        est_bit(j+1)=1;
    elseif (X(i)==-A)
        est_bit(j)=1;
        est_bit(j+1)=1;
    else
        est_bit(j)=1;
        est_bit(j+1)=0;
    end
j=j+2;
end
end