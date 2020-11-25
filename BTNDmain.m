
%%
%Main example for the matlab algorithme that perform BTND Decomposition

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
init = 2;

[F,V,~]=BTND(FC,K,param,init);

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
ylabel(['Subpgraph numbers'],'FontSize',15)
xlabel('Time (s)','FontSize',20)
grid on
ax = gca;
ax.YGrid = 'off';    
cb = colorbar;
ylabel(cb,'Activation level','Fontsize',20)
title(['Seizure ',num2str(i)],'FontSize',25)

end

%%
%Use BTND on your dataset

%Step 1 : Use S seizures recording as input (figure 1.(a) of the paper)
%Step 2 : Compute a dynamical graph for each seizure with a FC measure, they are converted to a list of FC matrix (figure 1.(b) and 1.(c) of the paper)
%Step 3 : Apply BTND method  (figure 1.(d) of the paper)
%Step 4 : Represent the subgraphs and the activation profils  (figure 1.(d) of the paper)


%%%%Step 1 : Use S seizures recording as input, 
%these recordings have to be already filtered in a bandwidth of interest,
% with : 
%Signal{1} = first recording of dimension ExN(1)
%Signal{2} = second recording of dimension ExN(2)
%...
%Signal{S} = last recording of dimension ExN(S)

% E the number of electrode contacts, had to stay the same for every
%recording. N(s), for s=1,...,S, is the number of temporal samples used to
%record each seizure (N(s) can be a different number for each seizure).


%%%%Step 2 : Compute a dynamical graph for each seizure with a FC measure, they are converted to a list of FC matrix (figure 1.(b) and 1.(c) of the paper)
%Three FC measures are already implemented :
% method='COR'; for the Pearson correlation,
% method='PLV'; for Phase Locking Value, 
% method='AEC'; for amplitude envelope correlation.

method = 'PLV';
Freq = 256; % sampling frequency;
TimeSegment = 4;% Size of the temporal segments (in s) used to compute a FC graph
Step = 1;% A graph is computed every 'Step' second
[DataFC] = FC_dynamic(Signal,method,Freq,TimeSegment,Step);



%%%%Step 3 : Apply BTND method
FC = DataFC;
K = 6;%Number of subgraph
lambda=0.4;% parameter for the sparsity level
gamma=0.2;%parameter for the temporal consistency of the activation profil
init = 1;%Number of different initialisation (use init=20 or more if you are sure of the parameters and want an ultimate result)

[F,V,cost]=BTND(FC,K,[lambda,lambda,gamma],init);


%%%%Step 4 : Represent the subgraphs and the activation profils
%Display Subgraph

%Put list of electrodes names or
ElectrodesNames =[];
for i =1:size(Signal{1},1)
ElectrodesNames{i} = num2str(i);
end



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
ylabel(['Subpgraph numbers'],'FontSize',15)
xlabel('Time (s)','FontSize',20)
grid on
ax = gca;
ax.YGrid = 'off';    
cb = colorbar;
ylabel(cb,'Activation level','Fontsize',20)
title(['Seizure ',num2str(i)],'FontSize',25)

end


