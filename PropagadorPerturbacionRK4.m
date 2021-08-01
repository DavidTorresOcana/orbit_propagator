function PropagadorPerturbacionRK4()

clc;
clear;
close all;
global mu NN N AT Ntramos A hip AR topo2 si2;
N=1200;
AR=1.5;

mu=3.986E5;
% Posicion inicial
rxo=1000+6378;
ryo=0;
rzo=0;
rvo=[rxo,ryo,rzo];
%Velocidad Inicial
vxo=0;
Incl=45;%deg

vyo=cos(deg2rad(Incl))*sqrt(mu/norm(rvo));
vzo=sin(deg2rad(Incl))*sqrt(mu/norm(rvo));  % definimos una circunferencia
vvo=[vxo,vyo,vzo];

%% Calculo de los valores
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
% Calculo del numero de tramos
Ntramos=NN*N/3;

%% Ploteo soluciones 
scrsz = get(0,'ScreenSize');
figure('Position',[5 50 scrsz(3)/1.006 scrsz(4)/1.2],'Name','3D animation','NumberTitle','off');
MainFig= findobj('Name','3D animation');
AxesMain = axes;
set(MainFig,'CurrentAxes',AxesMain)
set(MainFig,'Color',[0.05,0.05,0.05]);



PLOTT=hola(rxo,ryo);
hold all;
PLO=plot3(rxo,ryo,rzo);
PLO2=quiver3(rxo,ryo,rzo,vvo(1),vvo(2),vvo(3),1000,'r','MaxHeadSize',1,'LineWidth',2);
PLO3=plot3(rxo,ryo,rzo);
PLO4=quiver3(rxo,ryo,rzo,vvo(1),vvo(2),vvo(3),1000,'g','MaxHeadSize',1,'LineWidth',2);
PLO5=line([0 rxo],[0 ryo],[0 rzo],'Color',[1 0 0],'Marker','.','LineStyle','-','LineWidth',2);
PLOG=plot3(rxo,ryo,rzo,'g');
subsatelite=plot3(rxo,ryo,rzo,'r','LineWidth',4);
trazatiera=plot3(rxo,ryo,rzo,'g','LineWidth',4);
% legend('Tierra','Trayectoria','Velocidadx1000');
axis([-AR*A AR*A -AR*A AR*A -AR*A AR*A]);
set(AxesMain,'DrawMode','fast',...
        'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1,1,1],'Interruptible','off','PlotBoxAspectRatioMode','auto');

axis off;
hold off;
%% Mapa de Traza
figure('Position',[5 50 scrsz(3)/1.5 scrsz(4)/2],'Name','Markov','NumberTitle','off');
TraceMap= findobj('Name','TraceMap');
AxesTraceMap = axes;
set(TraceMap,'CurrentAxes',AxesTraceMap)

% % % I= imread('BigEarth2.jpg');
% % % imagesc(0:360,-90:90,I);

hold on

phi(1)=0;
lambda(1)=0;
traza=plot(rad2deg(lambda),rad2deg(phi),'r','LineStyle','none','Marker','x','linewidth',2);


%Medio Giro de la tierra
     for l=1:1:si2(2)
        if l>floor(12*3600*si2(2)/(24*3600))
            topo3(:,l)=topo2(:,ceil(l-12*3600*si2(2)/(24*3600)));
        else % l<=floor(to*si2(2)/(24*3600))
            topo3(:,l)=topo2(:,l+si2(2)-floor(12*3600*si2(2)/(24*3600)));
        end
     end
contour(0:359,-89:90,topo3,[0 0],'b')


set(AxesTraceMap,'XLim',[0 360],'YLim',[-90 90], ...
    'XTick',[0 60 120 180 240 300 360], ...
    'Ytick',[-90 -60 -30 0 30 60 90]);

axis equal;
axis([0 360 -90 90 0 0.1 0 1])


to=0;
desp=0;
Deltaan=[0 0 0];   % VEctor!!!!!!!!!!!!
Deltapan=[0 0 0];
rreal=rvo;
vreal=vvo;
RREAL(1,:)=rreal;
NUM=1;
for tramo=1:1:Ntramos
    %% Problema de Kepler
%    fprintf(['\n    TRAMO ',num2str(tramo),'\n']);
    [r2,v2,to2]=Propagador(rvo,vvo,to);
    
    % Guardamos posicion anterior para el problema Delta
    rvoan=rvo;
    vvoan=vvo;
    %Definicion condiciones iniciales para el siguiente tramo
    si=size(r2);
    
    At=to2-to; % At del tramo en concreto
    
    
    %Concatenamos matrices
    for p=1:1:si(2)
        r3(p+desp,:)=r2(p,:);
    end
    desp=desp+si(2);%Despues de usarlo !!


     %% Problema diferencia Delta. Segun el metodo de Encke . 
    % Se basa en decir que en cada tramo se produce un incremento de velocidad
    % finito y modelado segun una fuerza y discretizacion. La posicion se
    % modela como suma del problema de Kepler y el problema diferencia:
    %  V(tramo)+=V(tramo)-  + gamma; Donde gamma es el termino de aceleracion
    %  de perturbacion por el At del tramo.  El problema diferencia se modela:
    %  Delta=r(real)-r(kepler); Quedando la siguiente EDO:
    %  Delta''=mu/r(kepler)^3*([1-r(kepler)^3/r(real)^3]*(r(kepler)+delta)-delta)+gamma
    

   %% Integracion de la Orb diferencia. Met RK4
   if tramo <= 1
       gamma=Gamma_acell(norm(rreal),phi(tramo));
%        gamma=sph2cart(0,-3.986*10^5*cos(phi(tramo))*(1.321205405*10^5*sin(phi(tramo))*cos(phi(tramo))/norm(rreal).^2-6.569669248*10^5*((15/2)*sin(phi(tramo))^2*cos(phi(tramo))-(3/2)*cos(phi(tramo)))/norm(rreal).^3-2.665815540*10^9*((35/2)*sin(phi(tramo))^3*cos(phi(tramo))-(15/2)*sin(phi(tramo))*cos(phi(tramo)))/norm(rreal).^4)/norm(rreal).^2,3.986*10^5*((44040.18018*(-1/2+(3/2)*sin(phi(tramo))^2))/norm(rreal).^2-6.569669248*10^5*((5/2)*sin(phi(tramo))^3-(3/2)*sin(phi(tramo)))/norm(rreal).^3-2.665815540*10^9*(3/8+(35/8)*sin(phi(tramo))^4-(15/4)*sin(phi(tramo))^2)/norm(rreal).^4)/norm(rreal).^2-3.986*10^5*(-(88080.36036*(-1/2+(3/2)*sin(phi(tramo))^2))/norm(rreal).^3+1.970900774*10^6*((5/2)*sin(phi(tramo))^3-(3/2)*sin(phi(tramo)))/norm(rreal).^4+1.066326216*10^10*(3/8+(35/8)*sin(phi(tramo))^4-(15/4)*sin(phi(tramo))^2)/norm(rreal).^5)/norm(rreal));
   end
    n=tramo+1;
%             % Delta
%             k1=Deltapan;
%             k2=Deltapan+AT*k1/2;
%             k3=Deltapan+AT*k2/2;
%             k4=Deltapan+AT*k3;
%             Delta=Deltaan+AT/6*(k1+2*k2+2*k3+k4);
%             
%             % Deltap
%                 q=( dot(rvoan,Deltaan)+0.5*norm(Deltaan).^2 )/(norm(rvoan).^2);
%                 Aq=3*q*(1-5/2*q+35/6*q.^2 -315/24*q.^3);
%             k1=mu/norm(rvoan).^3*( Aq*(rvoan + (Deltaan) )  -(Deltaan) )+gamma;
%                 fprintf(' q=%f \n  Aq=%f \n k1=%f\n',q,Aq,k1);
%                 q=( dot(rvoan,Deltaan+AT*k1/2)+0.5*norm(Deltaan+AT*k1/2).^2 )/(norm(rvoan).^2);
%                 Aq=3*q*(1-5/2*q+35/6*q.^2 -315/24*q.^3);
%             k2=mu/norm(rvoan).^3*( Aq*(rvoan+ (Deltaan+AT*k1/2) )  -(Deltaan+AT*k1/2) )+gamma;
%                 fprintf(' q=%f \n  Aq=%f \n k2=%f\n',q,Aq,k2);
%                 q=( dot(rvoan,Deltaan+AT*k2/2)+0.5*norm(Deltaan+AT*k2/2).^2 )/(norm(rvoan).^2);
%                 Aq=3*q*(1-5/2*q+35/6*q.^2 -315/24*q.^3);
%             k3=mu/norm(rvoan).^3*( Aq*(rvoan+ (Deltaan+AT*k2/2) )  -(Deltaan+AT*k2/2) )+gamma;
%                 fprintf(' q=%f \n  Aq=%f\n  k3=%f\n',q,Aq,k3);
%                 q=( dot(rvoan,Deltaan+AT*k3)+0.5*norm(Deltaan+AT*k3).^2 )/(norm(rvoan).^2);
%                 Aq=3*q*(1-5/2*q+35/6*q.^2 -315/24*q.^3);
%             k4=mu/norm(rvoan).^3*( Aq*(rvoan+ (Deltaan+AT*k3) )  -(Deltaan+AT*k3) )+gamma;
%                 fprintf(' q=%f \n  Aq=%f \n k4=%f \n',q,Aq,k4);
%                 
%             Deltap=Deltapan+AT/6*(k1+2*k2+2*k3+k4);
            
%      fprintf(' Delta=%f %f %f \n  Deltap=%f %f %f \n',Delta,Deltap);
    % Integracion ODE45
    [~,u] = ode45(@(t,u) DIFF(t,u,rvoan,gamma) ,[0 AT],[Deltaan , Deltapan],odeset('RelTol',1e-6,'AbsTol',1e-6,'NormControl','on'));
    Delta=u(end,1:3);
    Deltap=u(end,4:6);
%     pause
%% Posicion real

rrealan=rreal; % Hay que ponerla antes del updateo de rreal
rreal=Delta+r2(si(2),:);% posicion real
%  fprintf(' rreal=%f %f %f \n',rreal);
RREAL(tramo+1,:)=rreal;
% Velocidad real
% vreal=vreal+gamma*(to2-to) %%%%%%%%%%%Verdaderos problemas para calcular la velocidad de la orbita perturbada :(
% vreal=v2+gamma*(to2-to) %%%%%%%%%%%Verdaderos problemas para calcular la velocidad de la orbita perturbada :(
[rrealprop,vreal,to2prop]=Propagador(rrealan,vreal,to);   %% Se estima la velocidad propagando la orbita perturbada desde la rrealan y con la vreal



 
    rvo=r2(si(2),:);
    vvo=v2;
    to=to2;
    Deltaan=Delta;
    Deltapan=Deltap;
%% Calculo de graficos
% Hacemos to2 menor que un dia para usarlo en las graficas
    while to2 >=24*3600
        to2=to2-24*3600;
    end
    
    if norm(rreal)<=6528 & norm(rreal)>6378
        fprintf(' REENTRADA INMINENTE!!\n\n\n');
    end
    if norm(rreal)<=6378
            fprintf( '   IMPACTO!!!!!!!!! \n\n\n');
            system('pause');
    end

    % Calculo de la traza
    phi(tramo)=atan2(rreal(3),sqrt(rreal(1).^2+rreal(2).^2));
    lambda(tramo)=-2*pi/(24*3600)*to+atan2(rreal(2),rreal(1));
    while lambda(tramo)>=2*pi
        lambda(tramo)=lambda(tramo)-2*pi;
    end
    while lambda(tramo)<0
        lambda(tramo)=lambda(tramo)+2*pi;
    end
    set(traza,'YData',rad2deg(phi),'XData',rad2deg(lambda));
    [RXp(tramo),RYp(tramo),RZ(tramo)]=sph2cart(lambda(tramo),phi((tramo)),6378);
    % Cambio de jes a ejes tierra que giran con ella
    for p=1:1:tramo
        RX(p)=cos(-2*pi/(24*3600)*to)*RXp(p)+sin(-2*pi/(24*3600)*to)*RYp(p);
        RY(p)=cos(-2*pi/(24*3600)*to)*RYp(p)-sin(-2*pi/(24*3600)*to)*RXp(p);
    end

    %Calculo de la trayectoria en ejes tierra
    [RXp2(tramo),RYp2(tramo),RZ2(tramo)]=sph2cart(lambda(tramo),phi((tramo)),norm(rreal));
        for p=1:1:tramo
        RX2(p)=cos(-2*pi/(24*3600)*to)*RXp2(p)+sin(-2*pi/(24*3600)*to)*RYp2(p);
        RY2(p)=cos(-2*pi/(24*3600)*to)*RYp2(p)-sin(-2*pi/(24*3600)*to)*RXp2(p);

        end
    
    % Traza subsatelite
    lambda2(tramo)=atan2(rreal(2),rreal(1));
    while lambda2(tramo)>=2*pi
        lambda2(tramo)=lambda2(tramo)-2*pi;
    end
    [RX3(tramo),RY3(tramo),RZ3(tramo)]=sph2cart(lambda2(tramo),phi((tramo)),6378);
    
    
%% Ploteo de graficos
    set(PLOG,'XData',RX2,'YData',RY2,'ZData',RZ2);
    set(PLO,'XData',r3(:,1),'YData',r3(:,2),'ZData',r3(:,3));
    set(PLO2,'XData',rvo(1),'YData',rvo(2),'ZData',rvo(3),'UData',vvo(1),'VData',vvo(2),'WData',vvo(3));
    set(PLO3,'XData',RREAL(:,1),'YData',RREAL(:,2),'ZData',RREAL(:,3));
    set(PLO4,'XData',rreal(1),'YData',rreal(2),'ZData',rreal(3),'UData',vreal(1),'VData',vreal(2),'WData',vreal(3));
    set(PLO5,'XData',[0 rreal(1)],'YData',[0 rreal(2)],'ZData',[0 rreal(3)]);
    set(subsatelite,'XData',RX3,'YData',RY3,'ZData',RZ3);
    set(trazatiera,'XData',RX,'YData',RY,'ZData',RZ);
    %Giro de la tierra (ploteo)
%     fprintf('    to2 (dias)= %f  \n',to2/(24*3600));
     for l=1:1:si2(2)
        if l>floor(to2*si2(2)/(24*3600))
            topo3(:,l)=topo2(:,ceil(l-to2*si2(2)/(24*3600)));
        else % l<=floor(to*si2(2)/(24*3600))
            topo3(:,l)=topo2(:,l+si2(2)-floor(to2*si2(2)/(24*3600)));
        end
     end
    set(PLOTT,'CData',topo3);
    
%         set(AxesMain,'CameraPosition',[2*rvo(1),2*rvo(2),0.6*A],'CameraTarget',rvo); % Change the camera position
%         set(AxesMain,'CameraPosition',[1*A*cos(  to2*2*pi/(24*3600) - deg2rad( 180 )  ),1*A*sin(  to2*2*pi/(24*3600) - deg2rad( 180 ) )  ,0.6*A],'CameraTarget',[0,0,0]); % Change teh camera position

        
%     campos([2*rvo(1),2*rvo(2),0.6*A]);
% campos([cos(2*pi/(24*3600)*to-pi)*A,sin(2*pi/(24*3600)*to-pi)*A,0.6*A]);
    drawnow;
%% Modelo de aceleracion para el siguiente tramo. Usando l
%   Usamos la perturbacion terrestre con terminos J hasta 4
       gamma=Gamma_acell(norm(rreal),phi(tramo));

% gamma=gamma; % Disable the preturbation
% fprintf('  GAMMMAAAAAA= %14f %14f %14f',gamma);
end

pause(3)

% % % % %  Ploteamos finalmente la traza
% % % % Traza_final=figure('Position',[5 50 scrsz(3)/1.006 scrsz(4)/1.2],'Name','Traza final','NumberTitle','off');
% % % % 
% % % % hold on
% % % % I= imread('BigEarth2.jpg');
% % % % imagesc(0:360,-90:90,I);
% % % % plot(rad2deg(lambda),rad2deg(phi),'r','LineStyle','none','Marker','x','linewidth',2);
% % % % hold off;
% % % % axis([0 360 -90 90]);
% % % % view(2);
end

function DIF=DIFF(t,u,rvoan,gamma)
global mu
if size(u,1)==6
    u=u';
end
    q=( dot(rvoan,u(1:3))+0.5*norm(u(1:3)).^2 )/(norm(rvoan).^2);
    Aq=3*q*(1-5/2*q+35/6*q.^2 -315/24*q.^3);
    k1=mu/norm(rvoan).^3*( Aq*(rvoan + (u(1:3)) )  -(u(1:3)) )+gamma;
    DIF=[u(4:6),k1]';
end
function [r,v,t]=Propagador(rvo,vvo,to)
%% Propagador de orbitas
global mu NN N Ntramos AT;

% Calculo parametro iniciales
ro=norm(rvo);
vo=norm(vvo);
%% Que hacemos con estos valores??
prd=dot(rvo,vvo);
sgo=prd/(sqrt(mu));
alfa=abs((2/ro)-(vo.^2)/mu);

Period=2*pi/sqrt(abs(mu*alfa.^3));

%% Los actualizamos o no?

X=0;
t=to;

for i=1:1:NN*N/Ntramos
    t=t+AT;
    while((t-to)>Period)
        to=to-Period;
    end
    
    % Determinacion de X
    xa=X;
	xb=(t-to)*sqrt(mu)/ro;
    X=SolveX(xa,xb,t-to,ro,sgo,alfa);
%     fprintf([' X ',num2str(X),'\n']);
    if X>2*pi/sqrt(mu)
        X=X-2*pi/sqrt(mu);
    end
    
    
    % Determinacion del vector posicion
    U=sub_U(X,alfa);
    
    r(i,:)=rvo*(1-U(2)/ro)+vvo*((ro*U(1)+sgo*U(2))/sqrt(mu));
 
    v=-sqrt(mu)*U(1)*rvo/(ro*norm(r(i,:)))+(1-U(2)/norm(r(i,:)))*vvo;
    
%     fprintf([num2str(r(i,:)),' \t  tiempo ',num2str(t),'\n']);
    
end



end

%% Funciones usadas

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
  U(1)=x-(1/6)*alfa.*x.^3+(1/120)*alfa.^2*x.^5-(1/5040)*alfa.^3*x.^7+(1/362880)*alfa.^4*x.^9;
  U(2)=(1/2)*x.^2-(1/24)*alfa.*x.^4+(1/720)*alfa.^2*x.^6-(1/40320)*alfa.^3*x.^8+(1/3628800)*alfa.^4*x.^10;
  U(3)=(1/6)*x.^3-(1/120)*alfa.*x.^5+(1/5040)*alfa.^2*x.^7-(1/362880)*alfa.^3*x.^9+(1/39916800)*alfa.^4*x.^11;
    
%   U(1)=sin(sqrt(alfa)*x)/sqrt(alfa);
%   U(2)=-(-1+cos(sqrt(alfa)*x))/alfa;
%   U(3)=(sqrt(alfa)*x-sin(sqrt(alfa)*x))/alfa.^(3/2);
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
load topo;
topo2=topo;
si2=size(topo);

[x y z] = sphere();
s = surface(x*6378,y*6378,z*6378,'FaceColor','texturemap','CData',topo,'LineStyle','none');

colormap(topomap1);
end
function F = funfo(x)
F=tan(2.7642/2)-sqrt((1+0.710862)/(1-0.710862))*tan(0.5*x);
end
function gamma=Gamma_acell(r,phi)
A(2)=-0.10e2 * (0.2633162373e10 * r ^ 2 - 0.7899487119e10 * r ^ 2 * sin(phi) ^ 2 + 0.2625554512e12 * sin(phi) ^ 3 * r - 0.1575332707e12 * sin(phi) * r + 0.2047995935e15 + 0.2389328591e16 * sin(phi) ^ 4 - 0.2047995935e16 * sin(phi) ^ 2) / r ^ 6;
A(1)=-0.3986e6 * cos(phi) / r ^ 2 * (0.1321205405e6 / r ^ 2 * sin(phi) * cos(phi) - 0.6586940570e6 / r ^ 3 * (0.15e2 / 0.2e1 * sin(phi) ^ 2 * cos(phi) - 0.3e1 / 0.2e1 * cos(phi)) - 0.2740252131e10 / r ^ 4 * (0.35e2 / 0.2e1 * sin(phi) ^ 3 * cos(phi) - 0.15e2 / 0.2e1 * sin(phi) * cos(phi)));

gamma=sph2cart(0,A(1),A(2));
end