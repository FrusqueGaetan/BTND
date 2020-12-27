function [Fn,Vn,cost] = BTND(FC,K,param,init)
%Perform de BTND decomposition "init" times (the first initialisation is
%always obtained via the SVD decomposition, the other are random).
%F,V and cost are the results of the best decomposition with : 
%F = FC subgraphs
%V = activation profils
%cost = cost function over iteration of the BTND Decomposition




if init > 1
    [F,V,cost] = BTND_Method(FC,K,param,'SVD');
   for i =1:(init-1)
    [Fbis,Vbis,costbis] = BTND_Method(FC,K,param,'RANDOM');
       if costbis(end)<cost(end)%Keep only the initialisation with the solution tha minimize the cost function
           %of the BTND problem
          F=Fbis;
          V=Vbis;
          cost=costbis;
       end
   end
else
   [F,V,cost] = BTND_Method(FC,K,param,'RANDOM'); 
end


%Normalisation step (enforce the FC subgraphs value to be between 0 and 1)
Fn = F./repmat(max(F),size(F,1),1);
Vn=[];
 for i =1:length(V)
    Vn{i} = V{i}.*repmat(max(F),size(V{i},1),1);
 end
        

end