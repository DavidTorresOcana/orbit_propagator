function PropagadorKepler()

clc;
clear;
close all;
global mu NN N AT Ntramos A hip AR topo2 si2;
N=600;
AR=1;

mu=3.986E5;
% Posicion inicial
rxo=24000;
ryo=0;
rzo=0;
rvo=[rxo,ryo,rzo];
%Velocidad Inicial
vxo=0.01*sqrt(mu/norm(rvo));
vyo=0.6*sqrt(mu/norm(rvo));
vzo=0.8*sqrt(mu/norm(rvo)); 
vvo=[vxo,vyo,vzo];

%% Calculo de las valores
ro=norm(rvo);
prd=dot(rvo,vvo);

sgo=prd/(sqrt(mu));

vo=norm(vvo);

alfa=abs((2/ro)-(vo.^2)/mu);

Period=2*pi/sqrt(abs(mu*alfa.^3));
AT=Period/N;

h=cross(rvo,vvo);

e=cross(vvo,h)/mu-rvo/ro;
a=norm(h).^2/(mu*(1-norm(e).^2));
nn=sqrt(mu/a.^3);
i=rad2deg(acos(h(3)/norm(h)));
Omega=rad2deg(atan2(h(1),-h(2)));
omega=rad2deg(atan2(e(3)/sin(i),e(1)*cos(Omega)+e(2)*sin(Omega)));
fo=acos((-ro+(norm(h).^2/mu))/(ro*norm(e)));
options=optimset('Display','iter');   % Option to display output
[Eo,fval] = fsolve(@funfo,-pi,options);
TAU=(norm(e)*sin(Eo)-Eo)/nn;

fprintf(' Bienvenido al programa propagador de Orbitas por tramos!!!\n');
fprintf(' Los parametros de su oprbita son: \n');
fprintf('alpha= %f a= %f \n sigma_0= %f  T= %f \n i= %f  Omega= %f  omega=%f \n e= %f TAU= %f  \n',alfa,a,sgo,Period,i,Omega,omega,norm(e),TAU);

%% Discriminacion tipo orbita

if ((2/ro)-(vo.^2)/mu>0.0)
    fprintf(' \n  Elliptic orbit!! \n');
    A=(1/alfa)*1;
    hip=0;
elseif ((2/ro)-(vo.^2)/mu<0.0)
    fprintf(' \n Hyperbolic orbit!! \n');
    A=(1/alfa)*8;
    hip=1;
else
    fprintf('\n   Parabolic orbit!!\n');
    A=(1/alfa)*3;
    hip=2;
end


NN=input('Introduce el numero de orbitas que quieres dar  ');
% Calculo del numero de tramos Ntramos=NN*N/Npropportramo;
Ntramos=NN*N/3;  % Se ve que solo funciona bien si el numero de pasos
                %   propagados en cada tramo es 3
if (floor(NN*N/Ntramos)~=NN*N/Ntramos | floor(Ntramos)~=Ntramos)
    fprintf('   No se puede propagar con pasos o tramos NO enteros \n  Introduzca otro numero de orbitas \n\n\n');
    system('pause');
    orbital();
end


%% Ploteo soluciones 
scrsz = get(0,'ScreenSize');
figure('Position',[5 50 scrsz(3)/1.006 scrsz(4)/1.2],'Name','3D animation','NumberTitle','off');
MainFig= findobj('Name','3D animation');
AxesMain = axes;
set(MainFig,'CurrentAxes',AxesMain)
set(MainFig,'Color',[0.05,0.05,0.05]);

PLOTT=hola(rxo,ryo);
hold all;
PLO3=line([0 rxo],[0 ryo],[0 rzo],'Color',[1 0 0],'Marker','.','LineStyle','-','LineWidth',2)
PLO=plot3(rxo,ryo,rzo);
PLOG=plot3(rxo,ryo,rzo);
PLO2=quiver3(rxo,ryo,rzo,vvo(1),vvo(2),vvo(3),1000,'r','MaxHeadSize',1,'LineWidth',2);
subsatelite=plot3(rxo,ryo,rzo,'r','LineWidth',4);
trazatiera=plot3(rxo,ryo,rzo,'g','LineWidth',4);
% legend('Tierra','Trayectoria','Velocidadx1000');
axis([-AR*A AR*A -AR*A AR*A -AR*A AR*A]);
set(AxesMain,'DrawMode','fast',...
        'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1,1,1],'Interruptible','off','PlotBoxAspectRatioMode','auto');

axis off;
hold off;
%% Mapa de Markov
figure('Position',[5 50 scrsz(3)/1.5 scrsz(4)/2],'Name','Markov','NumberTitle','off');
Markov= findobj('Name','Markov');
AxesMarkov = axes;
set(Markov,'CurrentAxes',AxesMarkov)
phi(1)=0;
lambda(1)=0;
traza=plot(rad2deg(lambda),rad2deg(phi),'r','LineStyle','none','Marker','x','linewidth',2);
hold on

%Medio Giro de la tierra
     for l=1:1:si2(2)
        if l>floor(12*3600*si2(2)/(24*3600))
            topo3(:,l)=topo2(:,ceil(l-12*3600*si2(2)/(24*3600)));
        else % l<=floor(to*si2(2)/(24*3600))
            topo3(:,l)=topo2(:,l+si2(2)-floor(12*3600*si2(2)/(24*3600)));
        end
     end
contour(0:359,-89:90,topo3,[0 0],'b')

set(AxesMarkov,'XLim',[0 360],'YLim',[-90 90], ...
    'XTick',[0 60 120 180 240 300 360], ...
    'Ytick',[-90 -60 -30 0 30 60 90]);
axis equal;
axis([0 360 -90 90 0 0.1 0 1])

%% Calculo de la posicion
desp=0;
to=0;
RX(1)=0;
RY(1)=0;
RZ(1)=0;
for tramo=1:1:Ntramos

   fprintf(['\n    TRAMO ',num2str(tramo),'\n']);
    [r2,v2,to2]=Propagador(rvo,vvo,to);
    
    %Definicion condiciones iniciales para el siguiente tramo
    si=size(r2);
    rvo=r2(si(2),:);
    vvo=v2;
    to=to2;
    % Hacemos to2 menor que un dia para usarlo en las graficas
    while to2>24*3600
        to2=to2-24*3600;
    end
    %Concatenamos matrices
    for p=1:1:si(2)
        r3(p+desp,:)=r2(p,:);
    end
    desp=desp+si(2);%Despues de usarlo !!
    
    if norm(rvo)<=6528 & norm(rvo)>6378
        fprintf(' REENTRADA INMINENTE!!\n\n\n');
    end
    if norm(rvo)<=6378
            fprintf( '   IMPACTO!!!!!!!!! \n\n\n');
            system('pause');
    end

    % Calculo de la traza
    phi(tramo)=atan2(rvo(3),sqrt(rvo(1).^2+rvo(2).^2));
    lambda(tramo)=-2*pi/(24*3600)*to+atan2(rvo(2),rvo(1));
    while lambda(tramo)>=2*pi
        lambda(tramo)=lambda(tramo)-2*pi;
    end
    while lambda(tramo)<0
        lambda(tramo)=lambda(tramo)+2*pi;
    end
    set(traza,'YData',rad2deg(phi),'XData',rad2deg(lambda));
    [RXp(tramo),RYp(tramo),RZ(tramo)]=sph2cart(lambda(tramo),phi(tramo),6378);
    % Cambio de jes a ejes tierra que giran con ella
    for p=1:1:tramo
        RX(p)=cos(-2*pi/(24*3600)*to)*RXp(p)+sin(-2*pi/(24*3600)*to)*RYp(p);
        RY(p)=cos(-2*pi/(24*3600)*to)*RYp(p)-sin(-2*pi/(24*3600)*to)*RXp(p);
    end

    %Calculo de la trayectoria en ejes tierra
    [RXp2(tramo),RYp2(tramo),RZ2(tramo)]=sph2cart(lambda(tramo),phi((tramo)),norm(rvo));
        for p=1:1:tramo
        RX2(p)=cos(-2*pi/(24*3600)*to)*RXp2(p)+sin(-2*pi/(24*3600)*to)*RYp2(p);
        RY2(p)=cos(-2*pi/(24*3600)*to)*RYp2(p)-sin(-2*pi/(24*3600)*to)*RXp2(p);

        end
    
    % Traza subsatelite
    lambda2(tramo)=atan2(rvo(2),rvo(1));
    while lambda2(tramo)>=2*pi
        lambda2(tramo)=lambda2(tramo)-2*pi;
    end
    while lambda2(tramo)<0
        lambda2(tramo)=lambda2(tramo)+2*pi;
    end
    [RX3(tramo),RY3(tramo),RZ3(tramo)]=sph2cart(lambda2(tramo),phi(tramo),6378);
    set(PLOG,'XData',RX2,'YData',RY2,'ZData',RZ2);
    
    set(PLO,'XData',r3(:,1),'YData',r3(:,2),'ZData',r3(:,3));
    set(PLO2,'XData',rvo(1),'YData',rvo(2),'ZData',rvo(3),'UData',vvo(1),'VData',vvo(2),'WData',vvo(3));
    set(PLO3,'XData',[0 rvo(1)],'YData',[0 rvo(2)],'ZData',[0 rvo(3)]);
    set(subsatelite,'XData',RX3,'YData',RY3,'ZData',RZ3);
    set(trazatiera,'XData',RX,'YData',RY,'ZData',RZ);
    
    %Giro de la tierra
     for l=1:1:si2(2)
        
        if l>floor(to2*si2(2)/(24*3600))
            topo3(:,l)=topo2(:,ceil(l-to2*si2(2)/(24*3600)));
        else % l<=floor(to*si2(2)/(24*3600))
            topo3(:,l)=topo2(:,l+si2(2)-floor(to2*si2(2)/(24*3600)));
        end
     end
    set(PLOTT,'CData',topo3);
    
    set(AxesMain,'CameraPosition',[2*rvo(1),2*rvo(2),0.6*A],'CameraTarget',rvo); % Change teh camera position
    
    
%     campos([2*rvo(1),2*rvo(2),0.6*A]);
%     campos([cos(2*pi/(24*3600)*to-pi)*A,sin(2*pi/(24*3600)*to-pi)*A,0.6*A]);

% pause(0.00005)
drawnow

end

pause(2);




%ploteamos finalmente la traza
Traza_final=figure('Position',[5 50 scrsz(3)/1.006 scrsz(4)/1.2],'Name','Traza final','NumberTitle','off');

hold on
I= imread('BigEarth2.jpg');
imagesc(0:360,-90:90,I);
plot(rad2deg(lambda),rad2deg(phi),'r','LineStyle','none','Marker','x','linewidth',2);
hold off;
axis([0 360 -90 90]);
view(2);
end


function [r,v,t]=Propagador(rvo,vvo,to)
%% Propagador de orbitas
global mu NN N Ntramos AT;

% Calculo parametro iniciales
ro=norm(rvo);
vo=norm(vvo);
%% Que hacemos con estos valores??
prd=dot(rvo,vvo);
sgo=prd /(sqrt(mu));
alfa=abs((2/ro)-(vo.^2)/mu);

Period=2*pi/sqrt(abs(mu*alfa.^3));

%% Los actualizamos o no?

X=0;
t=to;

for i=1:1:NN*N/Ntramos
    t=t+AT;
    while((t-to)>Period)
        t0=t0+Period;
    end
    
    % Determinacion de X
    xa=X;
	xb=(t-to)*sqrt(mu)/ro;
    X=SolveX(xa,xb,t-to,ro,sgo,alfa);
    fprintf([' X ',num2str(X),'\n']);
    if X>2*pi/sqrt(mu)
        X=X-2*pi/sqrt(mu);
    end
    
    
    % Determinacion del vector posicion
    U=sub_U(X,alfa);
    
    r(i,:)=rvo*(1-U(2)/ro)+vvo*((ro*U(1)+sgo*U(2))/sqrt(mu));
 
    v=-sqrt(mu)*U(1)*rvo/(ro*norm(r(i,:)))+(1-U(2)/norm(r(i,:)))*vvo;
    
    fprintf([num2str(r(i,:)),' \t  tiempo ',num2str(t),'\n']);
    
end



end

%% Funciones usadas
% function plottraza(phi,lambda,traza)
%     
%     campos('auto');
%     refreshdata;
%     drawnow;
% end
function X=SolveX(xa,xb,t2,ro,sgo,alfa)
   global mu NN N AT Ntramos; 
    
    %Hacemos xb sea mayor 
    j=1;
    U=sub_U(xb,alfa);
    F=ro*U(1)+sgo*U(2)+U(3)-t2*sqrt(mu);
    while(F<0)
        j=j+1;
        xb=xb*j;
        U=sub_U(xb,alfa);
        F=ro*U(1)+sgo*U(2)+U(3)-t2*sqrt(mu);
    end
    %Hacemos que xa sea menor
    U=sub_U(xa,alfa);
    F=ro*U(1)+sgo*U(2)+U(3)-t2*sqrt(mu);
    while(F>0)
        xa=xa-0.2*xa;
        U=sub_U(xa,alfa);
        F=ro*U(1)+sgo*U(2)+U(3)-t2*sqrt(mu);
    end
    
    while(j>0)
        X=(xa+xb)/2;
        U=sub_U(X,alfa);
        
        F=ro*U(1)+sgo*U(2)+U(3)-t2*sqrt(mu);
        if(X==xa || X==xb)
            return
        elseif (F==0)
            return
        elseif F<0
            xa=X;
        elseif F>0
            xb=X;
        end
        
    end
    
    
end
function U=sub_U(x,alfa)
    global hip;
 if hip==0
%   U(1)=x-(1/6)*alfa.*x.^3+(1/120)*alfa.^2*x.^5-(1/5040)*alfa.^3*x.^7+(1/362880)*alfa.^4*x.^9;
%   U(2)=(1/2)*x.^2-(1/24)*alfa.*x.^4+(1/720)*alfa.^2*x.^6-(1/40320)*alfa.^3*x.^8+(1/3628800)*alfa.^4*x.^10;
%   U(3)=(1/6)*x.^3-(1/120)*alfa.*x.^5+(1/5040)*alfa.^2*x.^7-(1/362880)*alfa.^3*x.^9+(1/39916800)*alfa.^4*x.^11;
%     
  U(1)=sin(sqrt(alfa)*x)/sqrt(alfa);
  U(2)=-(-1+cos(sqrt(alfa)*x))/alfa;
  U(3)=(sqrt(alfa)*x-sin(sqrt(alfa)*x))/alfa.^(3/2);
elseif hip==1
   U(1)=x+(1/6)*alfa.*x.^3+(1/120)*alfa.^2*x.^5+(1/5040)*alfa.^3*x.^7;
  U(2)=(1/2)*x.^2+(1/24)*alfa.*x.^4+(1/720)*alfa.^2*x.^6+(1/40320)*alfa.^3*x.^8;
  U(3)=(1/6)*x.^3+(1/120)*alfa.*x.^5+(1/5040)*alfa.^2*x.^7+(1/362880)*alfa.^3*x.^9;
    
%   U(1)=sinh(sqrt(alfa)*x)/sqrt(alfa);
%   U(2)=(-1+cosh(sqrt(alfa)*x))/alfa;
%   U(3)=-(sqrt(alfa)*x-sinh(sqrt(alfa)*x))/alfa.^(3/2);
else
  U(1)=x;
  U(2)=(1/2)*x.^2;
  U(3)=(1/6)*x.^3;
end
end
function s=hola(rxo,ryo)
global A AR topo2 si2 topomap1
cla reset;
load('topo.mat','topo','topomap1');
topo2=topo;
si2=size(topo);

[x y z] = sphere();
s = surface(x*6378,y*6378,z*6378,'FaceColor','texturemap','CData',topo,'LineStyle','none');

colormap(topomap1);

% axis([-AR*A AR*A -AR*A AR*A -AR*A AR*A]);
 

end
function F = funfo(x)
F=tan(2.7642/2)-sqrt((1+0.710862)/(1-0.710862))*tan(0.5*x);
end
