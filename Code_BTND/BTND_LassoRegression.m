function F = BTND_LassoRegression(Y,D,lambda,A)
%Performs lasso regression with FISTA acceleration

if nargin <4
A = (pinv(D)*Y)';

end


epsilon =10^-5;
maxIt=10000;
It=1;
diff=100;

Cd = D'*D;
l = norm(Cd,2);
Yd=Y'*D;
emax = score(D,A,Y,lambda);
e(1)=score(zeros(size(D)),zeros(size(A)),Y,0);




U=A;
theta=1;
F=A;
while(and(diff>epsilon,It<maxIt))


        A0=A;
        theta0=theta;
    
        U =U+1/l*(Yd-U*Cd); 
        
        A = SF2(U,1/l*lambda);
        theta = (1+sqrt(1+4*theta^2))/2;
        U = A+theta0/theta*(A-A0);


 
    It=It+1;
    e(It) = score(D,A,Y,lambda);
    
    if e(It) < emax
        F= A;
        emax=e(It);
    end
    
    diff=norm(A0-A,'fro')/numel(A);


end
end


function d = SF2(u,lambda)

d=max(u-lambda,0);

end


function s=score(D,A,Y,lambda)
s=1/2*norm(D*A'-Y,'fro')^2+lambda*sum(sum(abs(A)));
end