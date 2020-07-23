%%%%


%%%Main example for the matlab algorithme that perform BTND Decomposition

%Entry : list of FC matrices for one patient. 
load DataFC.mat


%Decomposition
%FC = list of FC matrices
%K = number of factors
%param = hyperparameters lambda, gamma=lambda and eta 
%init = number of different initialisation (the first initialisation is always obtained via the SVD decomposition, the other are random)

FC = DataFC;
K = 6;
param = [0.4,0.4,0.1];
init = 21;

[F,V,cost]=BTND(FC,K,param,init);

%F,V and cost are the results of the best decomposition with : 
%F = FC subgraphs
%V = activation profils
%cost = cost function over iteration of the BTND Decomposition




%Display results
load NameCoordinate.mat
ElectrodesNames = NameCoordinate;

%Display Subgraph
for i =1:K
f1=figure('Renderer', 'painters', 'Position', [1 1 700 1000]);
DISPLAY_graphcreate(F,i,0.2,ElectrodesNames,['Subgraph ',num2str(i)],f1)%%%
end


%Display activation profil
for i =1:length(V)
f2=figure('Renderer', 'painters', 'Position', [1 1 700 400])
imagesc(V{i}');
colormap(flip(pink))
set(gca,'Ytick',1:size(V{i},2),...
      'YTickLabel', 1:size(V{i},2),'FontSize',16);
set(gca,'Xtick',10:10:(size(V{i},1)),...
       'XTickLabel', 10:10:(size(V{i},1)),'FontSize',16);
ylabel(['Activation of the subgraphs'],'FontSize',15)
xlabel('Time (s)','FontSize',20)
grid on
ax = gca;
ax.YGrid = 'off';    
colorbar
title(['Seizure ',num2str(i)],'FontSize',25)

end


