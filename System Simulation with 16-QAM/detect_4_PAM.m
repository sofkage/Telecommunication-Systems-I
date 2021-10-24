function est_X = detect_4_PAM(Y,A)

Y_size = size(Y,2);
est_X = zeros(1,Y_size(1));
X = zeros(1,4);

X(1)=3*A;
X(2)=A;
X(3)=-A;
X(4)=-3*A;

distance=zeros(1,4);

for i=1:Y_size(1)
    distance = 0.*distance;
    for j=1:4
        distance(j) = ((X(j)-Y(i)).^2);
    end
    pos = min(distance);
    if (distance(1) == pos)
        est_X(i)= X(1);
    elseif (distance(2) == pos)
        est_X(i)= X(2);
    elseif (distance(3) == pos)
        est_X(i)= X(3);
    else
        est_X(i)= X(4);
    end
end
end