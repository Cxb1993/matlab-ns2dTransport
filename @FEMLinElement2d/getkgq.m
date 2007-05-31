function [kxxc,kyyc,nglec,v]=getkgq(t,mele,IEN,X,Y,Z)

%TELEMENT element class constructor.
%   t = TElement(m) creates a element object

%Name: getkgq
%Location: <path>/@FEMLinElement/
%Purpose: get element matrix with numerical integration - Gauss Quadrature

%ngleu: numero de graus de liberdade do elemento associados a u
%nglep: idem associados a p

nglec=3;
quadpoints=4;

% fator de correcao da matriz para dar o divergente correto
% refazer a matriz do elemento
v1=IEN(mele,1);
v2=IEN(mele,2);
v3=IEN(mele,3);
v=[v1;v2;v3];

A=0.5*det ([X(v2)-X(v1) X(v3)-X(v1);Y(v2)-Y(v1) Y(v3)-Y(v1)]);
jacobian = X(v3) * (Y(v1) - Y(v2)) + X(v1) * (Y(v2) - Y(v3)) + X(v2) * (-Y(v1) + Y(v3));


gqPoints= [ 0.3333333333 0.3333333333 0.3333333333;...
            0.6000000000 0.2000000000 0.2000000000;... 
            0.2000000000 0.6000000000 0.2000000000;... 
            0.20000000000 0.20000000000 0.6000000000];

gqWeights = [-0.5625000000 0.5208333333 0.5208333333 0.5208333333 0.5208333333];

phiJ = [ 0.3333333333 0.3333333333 0.3333333333;...
     	 0.6000000000 0.2000000000 0.2000000000;... 
         0.2000000000 0.6000000000 0.2000000000;... 
         0.20000000000 0.20000000000 0.6000000000];
    
dphiJdl1 = [1.0 0.0 -1.0;...
		   	1.0 0.0 -1.0;...
		    1.0 0.0 -1.0;...
            1.0 0.0 -1.0];

dphiJdl2 = [0.0 1.0 -1.0;...
		   	0.0 1.0 -1.0;...
		    0.0 1.0 -1.0;...
            0.0 1.0 -1.0];
        
  for k=1:quadpoints
    valx = 0.0;
    valy = 0.0;
    for i=1:nglec
      valx = valx + X(v(i)) * phiJ(k,i);
      valy = valy + Y(v(i)) * phiJ(k,i);
    end;
    localx(k) = valx;
    localy(k) = valy;
  end;
  

  for k=1:quadpoints
    valxl1 = 0.0;
    valxl2 = 0.0;
    valyl1 = 0.0;
    valyl2 = 0.0;
    for i=1:nglec
      valxl1 = valxl1 + X(v(i))* dphiJdl1(k,i);
      valxl2 = valxl2 + X(v(i))* dphiJdl2(k,i);
      valyl1 = valyl1 + Y(v(i))* dphiJdl1(k,i);
      valyl2 = valyl2 + Y(v(i))* dphiJdl2(k,i);
      end;
    dxdl1(k) = valxl1;
    dxdl2(k) = valxl2;
    dydl1(k) = valyl1;
    dydl2(k) = valyl2;
    end;

for k=1:quadpoints
	for i=1:nglec
		dphiJdx(k,i)=(dphiJdl1(k,i)*dydl2(k)-dphiJdl2(k,i)*dydl1(k))/jacobian;
		dphiJdy(k,i)=(-dphiJdl1(k,i)*dxdl2(k)+dphiJdl2(k,i)*dxdl1(k))/jacobian;
	end;
end;

kxxc=zeros(nglec,nglec);
kyyc=zeros(nglec,nglec);
for i=1:nglec
    for j=1:nglec
        for k=1:quadpoints
           kxxc(i,j)=kxxc(i,j)+dphiJdx(k,i)*dphiJdx(k,j)*jacobian*gqWeights(k);
           kyyc(i,j)=kyyc(i,j)+dphiJdy(k,i)*dphiJdy(k,j)*jacobian*gqWeights(k);
        end;
    end;
end;
