function [up,vp,cp] = convectlin(s,dt)
%modelo de documentacao a partir de:
%http://www.engin.umd.umich.edu/CIS/course.des/cis400/matlab/oop.html

%TERRAIN terrain class constructor.
%   t = Terrain(m) creates a terrain object from the matrix m,
%   containing the terrain data

%Name: Terrain
%Location: <path>/@Terrain
%Purpose: create the mother class Terrain

% modificado em 17/05/2006
% manter o centroide, criar uma matriz de conveccao aplicando nos vetores
% velocidade e pressao

[uc vc pc cc]=getBC(s.m);
IEN = getIEN(s.m);
X=getX(s.m);
Y=getY(s.m);
Z=getZ(s.m);

nele=size(IEN,1);
nnodes=size(X,1);
nvert=nnodes-nele;

us=s.us;
vs=s.vs;
cs=s.cs;
xp=X-us*dt;
yp=Y-vs*dt;
up=us;
vp=vs;
cp=cs;

convlin=sparse(nvert,nvert);
convlin2=sparse(nvert,nvert);
found=zeros(nvert,1);
found(s.m.idbcu)=1;
vetvert=find(found==0);
vetmele=tsearch(X(1:nvert),Y(1:nvert),IEN(:,1:3),xp(vetvert),yp(vetvert));
subvetvert1=find(vetmele>0);

for iii=1:size(subvetvert1,1)
    ii=vetvert(subvetvert1(iii));
    mele=vetmele(subvetvert1(iii));

    v1=IEN(mele,1);
    v2=IEN(mele,2);
    v3=IEN(mele,3);

    v=[v1;v2;v3];

    n1=[-(Y(v3)-Y(v2));X(v3)-X(v2)];
    n2=[-(Y(v1)-Y(v3));X(v1)-X(v3)];
    %n3=[-(Y(v2)-Y(v1));X(v2)-X(v1)];
    t12=[X(v1)-X(v2);Y(v1)-Y(v2)];
    t23=[X(v2)-X(v3);Y(v2)-Y(v3)];
    %t31=[X(v3)-X(v1);Y(v3)-Y(v1)];

    %vp1=[xp(ii)-X(v1);yp(ii)-Y(v1)];
    vp2=[xp(ii)-X(v2);yp(ii)-Y(v2)];
    vp3=[xp(ii)-X(v3);yp(ii)-Y(v3)];

    l1=(n1'*vp2)/(n1'*t12);
    l2=(n2'*vp3)/(n2'*t23);
    l3=1-l2-l1; %l3=(n3'*vp1)/(n3'*t31);

    % if((l1<=1)&(l1>=0)&(l2<=1)&(l2>=0)&(l3<=1)&(l3>=0))

    found(ii)=1;

    convlin(ii,v1)=l1;
    convlin(ii,v2)=l2;
    convlin(ii,v3)=l3;

end;

subvetvert2=find(found==0);
ii=subvetvert2;
jj=dsearch(X(1:nvert),Y(1:nvert),IEN(:,1:3),xp(ii),yp(ii));
convlin2=sparse(ii,jj,1,nvert,nvert);
found(ii)=1;
convlin=convlin+convlin2;

if(sum(found)-nvert~=0)
    pause;
end;

up=convlin*us(1:nvert);
vp=convlin*vs(1:nvert);
up(s.m.idbcu)=s.m.uc(s.m.idbcu);
vp(s.m.idbcv)=s.m.vc(s.m.idbcv);
ua = setCentroid(s,[up;zeros(nnodes-nvert,1);vp;zeros(nnodes-nvert,1)]);
up=ua(1:nnodes);
vp=ua(nnodes+1:2*nnodes);

cp=convlin*cs;


