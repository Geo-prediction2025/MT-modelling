%%%%%%%%%%%% aniso  model%%%%%%%

clear
clc
tic

%Unit side length
DX=[2000 2000 2000 2000 2000 1000 1000 1000 500 250 250 ...
    250 250 500 500 500 500 500 500 500 500 500 500 250 250 ...
    250 250 500 1000 1000 1000 2000 2000 2000 2000 2000];
DY=[2000 2000 2000 2000 2000 1000 1000 1000 500 250 250 ...
    250 250 500 500 500 500 500 500 500 500 500 500 250 250 ...
    250 250 500 1000 1000 1000 2000 2000 2000 2000 2000];
DZ1=[4000 2000 2000 1000 500 200 200 100];
DZ2=[100 200 200 500 500 500 500 500 ...
    500 500 500 500 500 500 ....
    500 500 1000 1000 1000 2000 3000 3000 3000 3000 3000 3000];%最深处30Km
DZ=[DZ1 DZ2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Air layer, vertical direction 1-8
%Abnormal body, x horizontal 12-25, y horizontal 12-25, vertical 17-22
%%%Parameters
Nair=size(DZ1,2);

NX=size(DX,2);
NY=size(DY,2);
NZ=size(DZ,2);

X=[0 cumsum(DX)]-sum(DX)/2;
Y=[0 cumsum(DY)]-sum(DY)/2;
Z=-[0 cumsum(DZ(Nair+1:end))];

%%%%%%%%%
%%%anomaly%%%%%
xa1=12;xa2=25;
ya1=12;ya2=25;
za1=9;za2=14;
%%%%%%%%
xx=[X(xa1) X(xa2+1) X(xa2+1) X(xa1) X(xa1)];
yy=[Y(ya1) Y(ya1) Y(ya2+1) Y(ya2+1) Y(ya1)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%aniso conductivity
pxyz=[1 5 9];
Pair=1e12;

%%%%anomaly aniso%%%%%%%%%%%%%
% rhoax=10;rhoay=15;rhoaz=20;
% srhoa1=[1/rhoax 0 0;0 1/rhoay 0;0 0 1/rhoaz];
% thetaS=10/180*pi;
% thetaD=20/180*pi;
% thetaL=30/180*pi;

rhoax=5;rhoay=5;rhoaz=20;%The anisotropy coefficient is 2.
srhoa1=[1/rhoax 0 0;0 1/rhoay 0;0 0 1/rhoaz];
thetaS=0/180*pi;
thetaD=0/180*pi;
thetaL=0/180*pi;


RzS1=[cos(-thetaS) sin(-thetaS) 0;-sin(-thetaS) cos(-thetaS) 0;0 0 1];
RxD1=[1 0 0;0 cos(-thetaD) sin(-thetaD);0 -sin(-thetaD) cos(-thetaD)];
RzL1=[cos(-thetaL) sin(-thetaL) 0;-sin(-thetaL) cos(-thetaL) 0;0 0 1];

RzS2=[cos(thetaS) sin(thetaS) 0;-sin(thetaS) cos(thetaS) 0;0 0 1];
RxD2=[1 0 0;0 cos(thetaD) sin(thetaD);0 -sin(thetaD) cos(thetaD)];
RzL2=[cos(thetaL) sin(thetaL) 0;-sin(thetaL) cos(thetaL) 0;0 0 1];

sigani=RzS1*RxD1*RzL1*srhoa1*RzL2*RxD2*RzS2;



%%%%surrounding rock aniso%%%%%%%%%%
% rhoax=80;rhoay=120;rhoaz=100;
% srhoa2=[1/rhoax 0 0;0 1/rhoay 0;0 0 1/rhoaz];
% thetaS=20/180*pi;
% thetaD=45/180*pi;
% thetaL=30/180*pi;

rhoax=100;rhoay=100;rhoaz=100;
srhoa2=[1/rhoax 0 0;0 1/rhoay 0;0 0 1/rhoaz];
thetaS=0/180*pi;
thetaD=0/180*pi;
thetaL=0/180*pi;

RzS1=[cos(-thetaS) sin(-thetaS) 0;-sin(-thetaS) cos(-thetaS) 0;0 0 1];
RxD1=[1 0 0;0 cos(-thetaD) sin(-thetaD);0 -sin(-thetaD) cos(-thetaD)];
RzL1=[cos(-thetaL) sin(-thetaL) 0;-sin(-thetaL) cos(-thetaL) 0;0 0 1];

RzS2=[cos(thetaS) sin(thetaS) 0;-sin(thetaS) cos(thetaS) 0;0 0 1];
RxD2=[1 0 0;0 cos(thetaD) sin(thetaD);0 -sin(thetaD) cos(thetaD)];
RzL2=[cos(thetaL) sin(thetaL) 0;-sin(thetaL) cos(thetaL) 0;0 0 1];

siganiD=RzS1*RxD1*RzL1*srhoa2*RzL2*RxD2*RzS2;

%%%%%%%%%%%%%
sig1d=zeros(2,9);
sig1d(1,pxyz)=1/Pair;
sig1d(2,:)=[siganiD(1,:),siganiD(2,:),siganiD(3,:)];

sig3d=[sigani(1,:),sigani(2,:),sigani(3,:)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
srho=zeros(9,NX*NY*NZ);
sigmani1d=zeros(NZ,9);

%%%%%%%%air%%%%%
for z=1:Nair
    for y=1:NY
        for x=1:NX
             index=(z-1)*NX*NY+(y-1)*NX+x;
             srho(:,index)=sig1d(1,:);
             sigmani1d(z,:)=sig1d(1,:);
        end
    end
end

%%%%%%%surrounding%%%%%%
for z=Nair+1:NZ
    for y=1:NY
        for x=1:NX
             index=(z-1)*NX*NY+(y-1)*NX+x;
             srho(:,index)=sig1d(2,:);
             sigmani1d(z,:)=sig1d(2,:);
        end
    end
end

%%%%%%%%anomaly%%%%%%%%
for z=za1+Nair:za2+Nair
    for y=ya1:ya2
        for x=xa1:xa2
             index=(z-1)*NX*NY+(y-1)*NX+x;
             srho(:,index)=sig3d(1,:);
        end
    end
end


%%%%%%%%%%%calculate boundarry%%%%%%%%
fre=0.1;

%%%%%%%%%%%%%%%%%
eTop1=[1,0];
eTop2=[0,1];
zNode=[1,cumsum(DZ)];

%%%%%%%%%%%1D analytical%%%%%%%%%%
for fm=1:size(fre,2)
    [Ex1t,Ey1t,Ez1t,Hx1t,Hy1t,Ex2t,Ey2t,Ez2t,Hx2t,Hy2t]=mt1DAniAnalyticxyEH(fre(fm),sigmani1d,zNode,eTop1,eTop2);
    
    Ex1tf(fm,:)=Ex1t(1,:);
    Ey1tf(fm,:)=Ey1t(1,:);
    Ez1tf(fm,:)=Ez1t(1,:);
   
    Ex2tf(fm,:)=Ex2t(1,:);
    Ey2tf(fm,:)=Ey2t(1,:);
    Ez2tf(fm,:)=Ez2t(1,:);  
end

[Ex1,Ey1,Ez1,Hx1,Hy1,Hz1,Ex2,Ey2,Ez2,Hx2,Hy2,Hz2]=MT3DvectanisoforT(srho,sig1d,fre,Nair,NX,NY,NZ,DX,DY,DZ,Ex1tf,Ey1tf,Ez1tf,Ex2tf,Ey2tf,Ez2tf);
toc

%%%cal Z%%
Temp=Hx1.*Hy2-Hx2.*Hy1;
Zxx=(Ex1.*Hy2-Ex2.*Hy1)./Temp;
Zxy=(Ex2.*Hx1-Ex1.*Hx2)./Temp;
Zyx=(Ey1.*Hy2-Ey2.*Hy1)./Temp;
Zyy=(Ey2.*Hx1-Ey1.*Hx2)./Temp;

%计cal apprarent resistivity
mu0=4*pi*1e-7;
for ff=1:size(fre,2)
    rhoxxT(:,:,ff)=abs((Zxx(:,:,ff)).^2*sqrt(-1)/(2*pi*fre(ff)*mu0));
    rhoxyT(:,:,ff)=abs((Zxy(:,:,ff)).^2*sqrt(-1)/(2*pi*fre(ff)*mu0));
    rhoyxT(:,:,ff)=abs((Zyx(:,:,ff)).^2*sqrt(-1)/(2*pi*fre(ff)*mu0));
    rhoyyT(:,:,ff)=abs((Zyy(:,:,ff)).^2*sqrt(-1)/(2*pi*fre(ff)*mu0));
    phxxT(:,:,ff)=-atan(imag(Zxx(:,:,ff))./real(Zxy(:,:,ff))).*180/pi;
    phxyT(:,:,ff)=-atan(imag(Zxy(:,:,ff))./real(Zxy(:,:,ff))).*180/pi;
    phyxT(:,:,ff)=-atan(imag(Zyx(:,:,ff))./real(Zyx(:,:,ff))).*180/pi;
    phyyT(:,:,ff)=-atan(imag(Zyy(:,:,ff))./real(Zxy(:,:,ff))).*180/pi;
end

save MT3DanisoDT2
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%surveying points%%%%%%%%%
X1=(X(2:end)-DX/2);
Y1=(Y(2:end)-DY/2);
%%%%%%%plot%%%%%%%%%%
%%%%%%%rangge-10km to 10km%%%%%%%%%%%
xt=4:NX-3;
yt=4:NY-3;

rhoxxDT2=rhoxxT(yt,xt,1);
rhoxyDT2=rhoxyT(yt,xt,1);
rhoyxDT2=rhoyxT(yt,xt,1);
rhoyyDT2=rhoyyT(yt,xt,1);

save rhoxxDT2 rhoxxDT2
save rhoxyDT2 rhoxyDT2
save rhoyxDT2 rhoyxDT2
save rhoyyDT2 rhoyyDT2

X2=X1(xt)./1000;
Y2=Y1(yt)./1000;

