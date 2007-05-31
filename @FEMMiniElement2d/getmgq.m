function [massele,kxx,kyy,kxy,kyx,gxele,gyele,ngleu,nglep,v]=getmgq(t,mele,IEN,X,Y,Z)

%TELEMENT element class constructor.
%   t = TElement(m) creates a element object

%Name: getmgq
%Location: <path>/@FEMMiniElement/
%Purpose: get matrix element with numerical integration - Gauss Quadrature

rule = 6;

%ngleu: numero de graus de liberdade do elemento associados a u
%nglep: idem associados a p

ngleu=4;
nglep=3;


% fator de correcao da matriz para dar o divergente correto
% refazer a matriz do elemento
v1=IEN(mele,1);
v2=IEN(mele,2);
v3=IEN(mele,3);
v4=IEN(mele,4);
v=[v1;v2;v3;v4];


A=0.5*det ([X(v2)-X(v1) X(v3)-X(v1);Y(v2)-Y(v1) Y(v3)-Y(v1)]);
jacobian = X(v3) * (Y(v1) - Y(v2)) + X(v1) * (Y(v2) - Y(v3)) + X(v2) * (-Y(v1) + Y(v3));

if ( rule == 1 )
    
   % formula         : 3 pontos
   % grau de precisao: 2
    
    r = 0.666666666666667;
    s = 0.166666666666667;
    t = 0.166666666666667; % t = 1 - (r + s)
    w = 0.333333333333333;

    xtab(1:3) =   [ r, s, s ];
    ytab(1:3) =   [ s, r, t ];
    ztab(1:3) =   [ t, t, r ];
    weight(1:3) = [ w, w, w ];

    norder = 3;


elseif ( rule == 2 )
    
   % formula         : 4 pontos
   % grau de precisao: 3

    a = 0.333333333333333;
    b = 0.600000000000000;
    c = 0.200000000000000;
    w = -0.562500000000000;
    q = 0.520833333333333;

    xtab(1:4) =   [ a, b, c, c ];
    ytab(1:4) =   [ a, c, b, c ];
    ztab(1:4) =   [ a, c, c, b ];
    weight(1:4) = [ w, q, q, q ];

    norder = 4;

elseif ( rule == 3 )
    
   % formula         : 6 pontos
   % grau de precisao: 3

    a = 0.659027622374092;
    b = 0.231933368553031;
    c = 0.109039009072877;
    w = 0.166666666666667;

    xtab(1:6) =   [ a, b, c, b, a, c ];
    ytab(1:6) =   [ b, a, b, c, c, a ];
    ztab(1:6) =   [ c, c, a, a, b, b ];
    weight(1:6) = [ w, w, w, w, w, w ];

    norder = 6;


elseif ( rule == 4 )
    
   % formula         : 7 pontos
   % grau de precisao: 5

    a = 0.333333333333333;
    b = 0.797426985353087;
    c = 0.470142064105115;
    d = 0.101286507323456;
    e = 0.059715871787770;
    w = 0.225000000000000;
    q = 0.125939180544827;
    r = 0.132394152788506;

    xtab(1:7) =   [ a, b, d, d, c, c, e ];
    ytab(1:7) =   [ a, d, b, d, c, e, c ];
    ztab(1:7) =   [ a, d, d, b, e, c, c ];
    weight(1:7) = [ w, q, q, q, r, r, r ];

    norder = 7;

   elseif ( rule == 5 )
       
   % formula         : 9 pontos
   % grau de precisao: 5  

    a = 0.124949503233232;
    b = 0.797112651860071;
    c = 0.437525248383384;
    d = 0.165409927389841;
    e = 0.037477420750088;
    w = 0.205950504760887;
    q = 0.063691414286223;
    
    xtab(1:9) =   [ a, c, c, b, d, d, e, e, b ];
    ytab(1:9) =   [ c, a, c, d, b, e, d, b, e ];
    ztab(1:9) =   [ c, c, a, e, e, b, b, d, d ];
    weight(1:9) = [ w, w, w, q, q, q, q, q, q ];

    norder = 9; 
    

elseif ( rule == 6 )
    
   % formula         : 12 pontos
   % grau de precisao: 6

    a = 0.873821971016996;
    b = 0.501426509658179;
    c = 0.636502499121399;
    d = 0.063089014491502;
    e = 0.249286745170910;
    f = 0.310352451033785;
    g = 0.053145049844816;
    w = 0.050844906370207;
    q = 0.116786275726379;
    r = 0.082851075618374;

    xtab(1:12) =   [ a, d, d, b, e, e, c, f, f, c, g, g ];
    ytab(1:12) =   [ d, a, d, e, b, e, f, c, g, g, c, f ];
    ztab(1:12) =   [ d, d, a, e, e, b, g, g, c, f, f, c ];
    weight(1:12) = [ w, w, w, q, q, q, r, r, r, r, r, r ];

    norder = 12;

    elseif ( rule == 7 )
    
   % formula         : 13 pontos
   % grau de precisao: 7

    a = 0.333333333333333;
    b = 0.479308067841923;
    c = 0.260345966079038;
    d = 0.869739794195568;
    e = 0.065130102902216;
    f = 0.638444188569809;
    g = 0.312865496004875;
    h = 0.086903154253160;
    w = -0.149570044467670;
    q = 0.175615257433204;
    r = 0.053347235608839;
    s = 0.077113760890257;
    
    xtab(1:13) =   [ a, b, c, c, d, e, e, f, f, g, g, h, h ];
    ytab(1:13) =   [ a, c, b, c, e, d, e, g, h, h, f, f, g ];
    ztab(1:13) =   [ a, c, c, b, e, e, d, h, g, f, h, g, f ];
    weight(1:13) = [ w, q, q, q, r, r, r, s, s, s, s, s, s ];    
   
    norder = 13;
    
    
else

    fprintf ( 1, '\n' );
    fprintf ( 1, 'TETRA_UNIT_SET - Fatal error!\n' );
    fprintf ( 1, '  Illegal value of RULE = %d\n', rule );
    error ( 'TETRA_UNIT_SET - Fatal error!' );

end;

quadpoints=norder;
L1 = xtab;
L2 = ytab;
L3=1-(L1+L2);
gqWeights = (1./2)*weight;

gqPoints = [ L1' L2' L3' ];

N4=27*L1'.*L2'.*L3';

phiJ =  [ L1'-N4/3 L2'-N4/3 L3'-N4/3 N4];

dphiJdl1 = [ 1-9*L2'.*L3'+9*L1'.*L2' -9*L2'.*L3'+9*L1'.*L2' -9*L2'.*L3'-(1-9*L1'.*L2') 27*L2'.*L3'-27*L1'.*L2' ];

dphiJdl2 = [ -9*L1'.*L3'+9*L1'.*L2' 1-9*L1'.*L3'+9*L1'.*L2' -9*L1'.*L3'-(1-9*L1'.*L2') 27*L1'.*L3'-27*L1'.*L2' ];

  for k=1:quadpoints
    valx = 0.0;
    valy = 0.0;
    for i=1:ngleu
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
    for i=1:ngleu
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

 %% formulacao proposta para o caso 3D - correta!   
    
  for k=1:quadpoints
    dxdl=[dxdl1(k) dxdl2(k);...
          dydl1(k) dydl2(k)];
    for i=1:ngleu
        dphiJdl = [dphiJdl1(k,i) dphiJdl2(k,i)];
        dphiJ = dphiJdl*inv(dxdl);
        dphiJdx(k,i) = dphiJ(1);
        dphiJdy(k,i) = dphiJ(2); 
    end;
  end;
    
%% formulacao proposto por Sao Carlos - correta!
%% tanto faz usar uma ou outra. Caso 3D validado!

%   for k=1:quadpoints
%     for i=1:ngleu
%       dphiJdx(k,i) = (dphiJdl1(k,i)*dydl2(k) - dphiJdl2(k,i)*dydl1(k))/jacobian;
%       dphiJdy(k,i) = (-dphiJdl1(k,i)*dxdl2(k) + dphiJdl2(k,i)*dxdl1(k))/jacobian;
%     end;
%   end;

kxx=zeros(ngleu,ngleu);
kyy=zeros(ngleu,ngleu);
kyx=zeros(ngleu,ngleu);
kxy=zeros(ngleu,ngleu);
massele=zeros(ngleu,ngleu);
for i=1:ngleu
    for j=1:ngleu
        for k=1:quadpoints
           massele(i,j)=massele(i,j)+phiJ(k,i)*phiJ(k,j)*jacobian*gqWeights(k);
           kxx(i,j)=kxx(i,j)+dphiJdx(k,i)*dphiJdx(k,j)*jacobian*gqWeights(k);
           kyy(i,j)=kyy(i,j)+dphiJdy(k,i)*dphiJdy(k,j)*jacobian*gqWeights(k);
           kyx(i,j)=kyx(i,j)+dphiJdy(k,i)*dphiJdx(k,j)*jacobian*gqWeights(k);
           kxy(i,j)=kxy(i,j)+dphiJdx(k,i)*dphiJdy(k,j)*jacobian*gqWeights(k);
        end;
    end;
end;

% montando operador GRADIENTE
gxele=zeros(ngleu,nglep);
gyele=zeros(ngleu,nglep);
for i=1:ngleu
    for j=1:nglep
        for k=1:quadpoints
            gxele(i,j)=gxele(i,j)-gqPoints(k,j)*dphiJdx(k,i)*jacobian*gqWeights(k);
            gyele(i,j)=gyele(i,j)-gqPoints(k,j)*dphiJdy(k,i)*jacobian*gqWeights(k);
        end;
    end;
end;



