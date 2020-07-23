function [] = DISPLAY_graphcreate(F,nn,lambda,NamePos,ttext,ff)

N = length(NamePos);
[~, iF, ~] = FPLV_infoG(N,NamePos);


Sel = sum(F>lambda,2)>=1;
zf = zeros(size(F(:,nn)));
zf(Sel) = (F(Sel,nn)-lambda)*(1/(1-lambda));

Total = iF(Sel,:);
Vectt = reshape(Total,1,numel(Total));
uVectt = unique(Vectt);

Gsel = ismember(1:N,uVectt);

NamePosNow = {NamePos{Gsel}};

Gfc = GetUnderDiag(zf,N);


GNow = Gfc(Gsel,Gsel);


[xcircle,ycircle] = circ(length(uVectt));


h = plot(xcircle, ycircle,'g.','MarkerSize',25);
axis square;
hold on



t = atan2(xcircle,ycircle);
textt=[];
longlabel = 1.1;
for i =1:length(uVectt)
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

    
       

color = flip(hot(125));
Ncolorbar = 100;
color = color(2:Ncolorbar,:);
for i =1:(size(GNow,1)-1)
    for j =(i+1):size(GNow,1)
        if(GNow(i,j)>0)

            u = [xcircle(i),ycircle(i)];
            v = [xcircle(j),ycircle(j)];
            
            if and(abs(u(1)+v(1))<0.0000001,abs(u(2)+v(2))<0.0000001)%diametric points
               line([u(1),v(1)],[u(2),v(2)],'color',color(min(floor(Ncolorbar*GNow(i,j))+1,Ncolorbar-1),:),'linewidth',2)
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
%min(floor(Ncolorbar*GNow(i,j)),Ncolorbar)
               %h = plot(xunit,yunit,'color',color(min(floor(N*GNow(i,j))+1,101),:),'linewidth',2);
               h = plot(xunit,yunit,'color',color(min(floor(Ncolorbar*GNow(i,j))+1,Ncolorbar-1),:),'linewidth',2);

            end
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


function [xcircle,ycircle] = circ(N)


th = 2*(0:pi/N:(pi-pi/N));
r=1;
x=0;
y=0;
xcircle = r * cos(th) + x;
ycircle = r * sin(th) + y;

end


function [iG, iF, iP] = FPLV_infoG(n,Position,Location,FC,col)

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


if nargin >1
    iP=[];
for i =1:nchoosek(n,2)
    iP{i} = [Position{iA(i)},'--',Position{iB(i)}];
end

if nargin > 3
  if nargin==4
   for k = 1:size(FC,1)
    plot3(Location(FC(k,:),1),Location(FC(k,:),2),Location(FC(k,:),3),'r')
    xlabel('x')
    ylabel('y')
    zlabel('z')
    hold on
   end 
  elseif nargin==5
   for k = 1:size(FC,1)
    plot3(Location(FC(k,:),1),Location(FC(k,:),2),Location(FC(k,:),3),col)
    xlabel('x')
    ylabel('y')
    zlabel('z')
    hold on
   end 
  end
    
end

if nargin>2
    
    plot3(Location(:,1),Location(:,2),Location(:,3),'b.')
    text(Location(:,1),Location(:,2),Location(:,3),Position)
end

end

end


function V =  GetUnderDiag(M,k) % M matrice d'adjacence

if(k == 0)
S = size(M);

sV = nchoosek(S(1),2);

V = zeros(sV,1);

a = 1;
for i = 1:(S(1)-1)
    for j = (i+1):S(2)
        
         V(a) = M(i,j);
         a = a +1;
    end
end
else
    
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
end



