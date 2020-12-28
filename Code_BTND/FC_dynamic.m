function [FC] = FC_dynamic(Signal,Method,Fr,Seg,Step)

if size(Signal,1)==1
Nseiz = length(Signal);
SignalS=Signal;
else
Nseiz=1;
SignalS{1}=Signal;
end

FC=[];

for i =1:Nseiz

Por = floor(Seg*Fr);
Pas = (Por+1):(floor(Step*Fr)):(size(SignalS{i},2)-Por);%sample location to apply the window

    if strcmp(Method,'PLV')
    
        FC{i} = PLV(SignalS{i},Por,Pas);

    elseif strcmp(Method,'COR')
    
        FC{i} = Corr(SignalS{i},Por,Pas);
    
    elseif strcmp(Method,'AEC')
    
        FC{i} = Aec(SignalS{i},Por,Pas);

    else
        warning('Use a valid FC measure, PLV was used by default')
        FC{i} = PLV(SignalS{i},Por,Pas);
    end


end

end



function FC=PLV(S,Po,Pa)
Nelec = size(S,1);
nC = nchoosek( Nelec,2);
iA = squeeze(zeros(1,nC));
iB = squeeze(zeros(1,nC));
a = 1;
b = ( Nelec-1);
for i = 1:( Nelec-1)  
    iA(a:b) = repmat(i,( Nelec-i),1); 
    iB(a:b) = (i+1): Nelec;
    a = a+ Nelec-i;
    b = b+ Nelec-i-1;
end
FC = zeros(nchoosek(size(S,1),2),length(Pa));
Psi = angle(hilbert(S')');
for i =1:length(Pa)
    d1=Pa(i)-Po;
    d2=Pa(i)+Po;
    FC(:,i) = abs(1/(2*Po)*sum(exp(1i*(Psi(iA,d1:d2)-Psi(iB,d1:d2))),2));
end
end

function FC=Corr(S,Po,Pa)
FC = zeros(nchoosek(size(S,1),2),length(Pa));
for i =1:length(Pa)
    d1=Pa(i)-Po;
    d2=Pa(i)+Po;
    M = bsxfun(@minus,S(:,d1:d2),mean(S(:,d1:d2),2) )*...
        bsxfun(@minus,S(:,d1:d2),mean(S(:,d1:d2),2) )';
    D = diag(diag(M).^(-1/2));
    C = abs(D*M*D);
    FC(:,i) = GetUnderDiag(C,0);
end
end

function FC=Aec(S,Po,Pa)
FC = zeros(nchoosek(size(S,1),2),length(Pa));
Amp = abs(hilbert(S')');
for i =1:length(Pa)
    d1=Pa(i)-Po;
    d2=Pa(i)+Po;
  M = bsxfun(@minus,Amp(:,d1:d2),mean(Amp(:,d1:d2),2) )*...
        bsxfun(@minus,Amp(:,d1:d2),mean(Amp(:,d1:d2),2) )';
    D = diag(diag(M).^(-1/2));
    C = abs(D*M*D);
    FC(:,i) = GetUnderDiag(C,0);
end
end

