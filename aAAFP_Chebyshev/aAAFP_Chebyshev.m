clear all
% July 23, 2025
% this code is used to compute: C(a,b,m)*b^(m+1) and epsi(a,b,m)b^(m+1)
% A(:,1) stores a and A(:,2) stores b
% m stores the values of m
A=[0.3,0.9;1.5,3;2,5;10,30;20,50];
m=[2, 4, 10, 15];
for i=1:4
    t=m(i);
    for j=1:5

        a=A(j,1);
        b=A(j,2);

        d=abs((2*a.*b-a-b)./(b-a));
        if a<1
            tau=sqrt((b.*(1-a))./(a.*(1-b)));
        else
            tau=sqrt((a.*(b-1))./(b.*(a-1)));
        end
        y1=d+sqrt(d.^2-1);
        y2=d-sqrt(d.^2-1);

        eps=2./(y1.^t+y2.^(t)); % Chebyshev


        feps=2*((tau-1)./(tau+1)).^(t).*(b.^(t+1));

        fc=eps.*(b.^(t+1));
        C(i,j)=fc;
        E(i,j)=feps;
    end
end
C
E
