function V = BTND_FusedRegression(Y,D,gamma,eta,A)
   
if nargin <4
   A = (pinv(D)*Y)';
end

epsilon =10^-6;
maxIt=100;
It=1;
diff=100;


Cd = D'*D;
Bd = (D'*Y)';

while(and(diff>epsilon,It<maxIt))
    
    oA=A;
    
    for i=1:size(D,2)
        
        u = Cd(i,i)^(-1)*(Bd(:,i)-A*Cd(:,i))+A(:,i);

        A(:,i) = BCOPA_FusedProjection(u,gamma,eta);%Use algorithme from https://lcondat.github.io/
         
    end
   
    It=It+1;
    diff=norm(oA-A,'fro');

end

V=A;

end




function nu = BCOPA_FusedProjection(u,gamma,eta)

ScoFirst=Score(SF2(u,0),gamma, eta, 1);
gammaMax=10000;
tau=1;
if ScoFirst <=tau
     nu=SF2(u,0);
else
nu=SF2(TV_Condat_v2(u,0),0)/(1);%Use algorithme from https://lcondat.github.io/
ResNow=Score(nu,gamma, eta, 1);
   if ResNow <= tau

   elseif isnan(ResNow)
       
   nu=zeros(size(u));
       
   else
       i = 0;
       while(ResNow>tau)
           ResBefore=ResNow;
           i=i+1;
           nu=SF2(TV_Condat_v2(u,eta*i),gamma*i)/(1+i);
           ResNow=Score(nu,gamma, eta, 1); 
           if(i>gammaMax)
               'Problem gammaMax too small'
               ResNow=tau+1;
               break
           end
       end
       
       
       if i>=gammaMax
           
       else
          
          coeff=-1;
          while(coeff>-12)
              i=i-10^(coeff+1);
              ResNow=ResBefore;
           while(ResNow>tau)
               ResBefore=ResNow; 
               i=i+10^(coeff);
               nu=SF2(TV_Condat_v2(u,eta*i),gamma*i)/(1+i);
               ResNow=Score(nu,gamma, eta, 1); 
           end
           coeff=coeff-1;
          end
          
          
       end
       
   end
    
    
end

end

function d = SF2(u,lambda)
d=max(u-lambda,0);
end

function sco = Score(F,gamma, eta, t)
    sco=gamma*t*l1(F)+eta*t*TV(F)+t/2*l2(F);
end
        
function n = TV(u)
n=sum(abs(u(2:end)-u(1:(end-1))));
end   

function n = l1(u)
n=sum(sum(abs(u)));
end

function n = l2(u)
n=sum(sum(u.^2));
end
