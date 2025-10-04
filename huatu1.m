%%%%%%%plot%%%%%%%%%%%%%
clear
clc

%% mesh

DX=[2000 2000 2000 2000 2000 1000 1000 1000 500 250 250 ...
    250 250 500 500 500 500 500 500 500 500 500 500 250 250 ...
    250 250 500 1000 1000 1000 2000 2000 2000 2000 2000];
DY=[2000 2000 2000 2000 2000 1000 1000 1000 500 250 250 ...
    250 250 500 500 500 500 500 500 500 500 500 500 250 250 ...
    250 250 500 1000 1000 1000 2000 2000 2000 2000 2000];
DZ1=[4000 2000 2000 1000 500 200 200 100];
DZ2=[100 200 200 500 500 500 500 500 ...
    500 500 500 500 500 500 ....
    500 500 1000 1000 1000 2000 3000 3000 3000 3000 3000 3000];
DZ=[DZ1 DZ2];

X=[0 cumsum(DX)]-sum(DX)/2;
Y=[0 cumsum(DY)]-sum(DY)/2;

NX=size(DX,2);
NY=size(DY,2);

%%%%%%%%
xa1=12;xa2=25;
ya1=12;ya2=25;
za1=9;za2=14;
%%%%%%%%
xx=[X(xa1) X(xa2+1) X(xa2+1) X(xa1) X(xa1)];
yy=[Y(ya1) Y(ya1) Y(ya2+1) Y(ya2+1) Y(ya1)];

%%%%%%%%%%%%%%%%%
X1=(X(2:end)-DX/2);
Y1=(Y(2:end)-DY/2);


xt=4:NX-3;
yt=4:NY-3;

X2=X1(xt)./1000;
Y2=Y1(yt)./1000;

%% data
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%FD%%%%
load rhoxxtananiso01hz
load rhoxytananiso01hz
load rhoyxtananiso01hz
load rhoyytananiso01hz

rhoxxtan=(rhoxx(xt,yt,1)).';
rhoxytan=(rhoxy(xt,yt,1)).';
rhoyxtan=(rhoyx(xt,yt,1)).';
rhoyytan=(rhoyy(xt,yt,1)).';

%%%%%%the first boundary%%%%
load rhoxxDF
load rhoxyDF
load rhoyxDF
load rhoyyDF

%%%%%%the mixed boundary%%%%
load rhoxxDT
load rhoxyDT
load rhoyxDT
load rhoyyDT


load rhoxxDT2
load rhoxyDT2
load rhoyxDT2
load rhoyyDT2




%% plot

figure(1)
subplot(2,2,1)
pcolor(X2,Y2,rhoxyDT);
shading interp
ylabel('Y(km)')
 xlabel('X(km)')
title('\lambda=1')
h1=colorbar('vert');
colormap('Jet')
 set(h1,'Position',[0.94 0.62 0.015 0.24])
h2=get(h1,'title');
set(h2,'string','\rho_x_y(\Omegam)','fontsize',10);
rhoscalexy=[50 120];
caxis(rhoscalexy)
hold on;
plot(yy/1000,xx/1000,'w--','linewidth',1.5)

subplot(2,2,3)
pcolor(X2,Y2,rhoyxDT);
shading interp
 ylabel('Y(km)')
xlabel('X(km)')
h1=colorbar('vert');
colormap('Jet')
 set(h1,'Position',[0.94 0.12 0.015 0.24])
h2=get(h1,'title');
set(h2,'string','\rho_y_x(\Omegam)','fontsize',10);
rhoscaleyx=[50 120];
caxis(rhoscaleyx)
hold on;
plot(yy/1000,xx/1000,'w--','linewidth',1.5)
% 


subplot(2,2,2)
pcolor(X2,Y2,rhoxyDT2);
shading interp
title('\lambda=2')
 ylabel('Y(km)')
 xlabel('X(km)')
h1=colorbar('vert');
colormap('Jet')
set(h1,'Position',[0.94 0.62 0.015 0.24])
h2=get(h1,'title');
set(h2,'string','\rho_x_y(\Omegam)','fontsize',10);
rhoscalexx=[50 120];
caxis(rhoscalexx)
hold on;
plot(yy/1000,xx/1000,'w--','linewidth',1.5)
% 
subplot(2,2,4)
pcolor(X2,Y2,rhoyxDT2);
shading interp
xlabel('X(km)')
 ylabel('Y(km)')
h1=colorbar('vert');
colormap('Jet')
set(h1,'Position',[0.94 0.12 0.015 0.24])
h2=get(h1,'title');
set(h2,'string','\rho_y_x(\Omegam)','fontsize',10);
rhoscaleyy=[50 120];
caxis(rhoscaleyy)
hold on;
plot(yy/1000,xx/1000,'w--','linewidth',1.5)



