function [] = DISPLAY_graphcreate(F,nn,lambda,NamePosNow,ttext,ff)

N = length(NamePosNow);

%Select only value above lambda
F = F.*(F>lambda);
zf = (F(:,nn)-lambda)*(1/(1-lambda));

%Generate graph adjency matrix
GNow =  ConvertToAdj(zf,N);
%Generate location of the nodes of the circular graph
[xcircle,ycircle] = circ(N);


%Nodes location and nodes names
h = plot(xcircle, ycircle,'g.','MarkerSize',25);
axis square;
hold on
t = atan2(xcircle,ycircle);
textt=[];
longlabel = 1.1;
for i =1:N
    NameReduced = NamePosNow{i};
    textt{i} = text(xcircle(i).*longlabel, ycircle(i).*longlabel, NameReduced, 'FontSize',14,'Color','b');
    if xcircle(i) > 0
        textt{i}.Rotation = fix(180*(-t(i)/pi +0.5 ));
        textt{i}.HorizontalAlignment = 'right';
    else
      textt{i}.Rotation = fix(180*(-t(i)/pi -0.5 ));
    end
end
xlim([-longlabel,longlabel])
ylim([-longlabel,longlabel])
set(gca,'XTick',[],'YTick',[])
    set(gca,'visible','off')

    
       
%Link representation on a Poincar� disc
color = flip(hot(125));
Ncolorbar = 100;
color = color(20:120,:);
[~,Nord] = sort(F(:,nn));
[~, iF] = FPLV_infoG(N);

for j = 1:size(F,1)
    i=Nord(j);
        if(GNow(iF(i,1),iF(i,2))>0)

            u = [xcircle(iF(i,1)),ycircle(iF(i,1))];
            v = [xcircle(iF(i,2)),ycircle(iF(i,2))];
            
            if and(abs(u(1)+v(1))<0.0000001,abs(u(2)+v(2))<0.0000001)
               line([u(1),v(1)],[u(2),v(2)],'color',color(min(floor(Ncolorbar*GNow(iF(i,1),iF(i,2)))+1,Ncolorbar-1),:),'linewidth',2)
            else
               x0 = -(u(2)-v(2))/(u(1)*v(2)-u(2)*v(1));
               y0 = (u(1)-v(1))/(u(1)*v(2)-u(2)*v(1));
               r0= sqrt(x0^2 + y0^2 - 1);
               thetaLim(1) = atan2(u(2)-y0,u(1)-x0);
               thetaLim(2) = atan2(v(2)-y0,v(1)-x0);
               if u(1) >= 0 && v(1) >= 0 
                  theta = [linspace(max(thetaLim),pi,100),linspace(-pi,min(thetaLim),100)];
               else
                  theta = linspace(thetaLim(1),thetaLim(2),200);
               end
               
               xunit = r0 * cos(theta) + x0;
               yunit = r0 * sin(theta) + y0;
               h = plot(xunit,yunit,'color',color(min(floor(Ncolorbar*GNow(iF(i,1),iF(i,2)))+1,Ncolorbar-1),:),'linewidth',2);

            end
        end

end

    

    ax = gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset; 
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3);
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left bottom ax_width (ax_height-0.05)];  
    
   title(ttext,'FontSize',36)

   annotation(ff,'textbox',...
    [0.0181428571428571 0.941823899371069 0.384714285714286 0.0518867924528302],...
    'String',ttext,'FontSize',36,...
    'FitBoxToText','on','LineStyle','none');


   
hold off
end


function [xcircle,ycircle] = circ(N) % Location of node on a circle


th = 2*(0:pi/N:(pi-pi/N));
r=1;
x=0;
y=0;
xcircle = r * cos(th) + x;
ycircle = r * sin(th) + y;

end


function [iG, iF] = FPLV_infoG(n)%iG: corresponding FC vector index on the Adj matrix
                                 %iF: corresponding Adj matrix index on the FC vector

iG = GetUnderDiag(1:nchoosek(n,2),n);


nC = nchoosek(n,2);
iA = squeeze(zeros(1,nC));
iB = squeeze(zeros(1,nC));
a = 1;
b = (n-1);
for i = 1:(n-1)  
    iA(a:b) = repmat(i,(n-i),1); 
    iB(a:b) = (i+1):n;
    a = a+n-i;
    b = b+n-i-1;
end

iF = [iA',iB'];

end


function V =  ConvertToAdj(M,k) % Convert a vector with nchoosek(N,2) element to an adjacency matrix size NxN

   M = squeeze(M);
   V = zeros(k,k);
   
   a = 1;
   for d = 1:k
       for j = 1:(k-d)
        
               V(d+j,d) = M(a);
               V(d,d+j) = M(a);
 
           a = a+1;
       end
   end
   
end



function V =  GetUnderDiag(M,k) %Convert adjacency matrix to vector

    
   M = squeeze(M);
   V = zeros(k,k);
   
   a = 1;
   for d = 1:k
       for j = 1:(k-d)
        
               V(d+j,d) = M(a);
               V(d,d+j) = M(a);

           a = a+1;
       end
   end
   
    
end