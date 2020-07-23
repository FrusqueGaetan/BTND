function [F,V,e] = BTND_Method(T,K,param,InitType)
%Dimension of the list of matrices
S{1} = size(T{1},1);
S{3} = length(T);
S{2} = zeros(1,S{3});
for i =1:S{3}
  S{2}(i) = size(T{i},2);   
end


%hyperparameters
lambda = param(1);
gamma = param(2);
eta = param(3);

Zeta=[];
ZetaMatrix=[];
for i =1:S{3}
    Zeta(i) = (norm(T{1},'fro')^2)/(norm(T{i},'fro')^2);
    ZetaMatrix = [ZetaMatrix, sqrt(Zeta(i))*ones(1,S{2}(i))];
end
lambda = S{3}*lambda;

TT = horzcat(T{:}).*repmat(ZetaMatrix,S{1},1);


%Initialisation of the variable F and V
if strcmp(InitType,'SVD')
   [F,~,~]= svds(TT,K);
else
   F= rand(size(TT,1),K);
end

for i =1:S{3}
V{i} = rand(S{2}(i),K);
end 

%Initialisation of the algorithme
epsilon =10^-4;
maxIt=300;
It=1;
diff=100;
e(1) = Score(T,F,V,lambda,S{3},Zeta);

%Algorithme

    while(and(diff>epsilon,It<maxIt))
%Alternate two steps   
        
        %First : fix F, perform a regression over Fused Lasso constraint to
        %find V
        for i =1:S{3}
             V{i} = BTND_FusedRegression(T{i},F,...
                 gamma*(size(V{i},1)/size(V{1},1))^(-1),eta*(size(V{i},1)/size(V{1},1))^(-1),V{i});
        end  
        
        %Second : fix V, perform a lasso regression to find F              
        VV = vertcat(V{:}).*repmat(ZetaMatrix',1,K);
        F = BTND_LassoRegression(TT',VV,lambda,F);
        
        %Compute the cost function
        It=It+1;
        e(It) = Score(T,F,V,lambda,S{3},Zeta);
        diff=(-e(It)+e(It-1))/e(2);
    end
        
end


function Sco = Score(X,F,V,lambda,S,Cb)
Sco=lambda*l1(F);
for i =1:S
    Sco = Sco+1/2*Cb(i)*l2(X{i}-F*V{i}');
end
end



function n = l1(u)
n=sum(sum(abs(u)));
end

function n = l2(u)
n=sum(sum(u.^2));
end





