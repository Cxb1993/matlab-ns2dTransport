function s = init(s)
%modelo de documentacao a partir de:
%http://www.engin.umd.umich.edu/CIS/course.des/cis400/matlab/oop.html

%SIMULATOR simulator class constructor.
%   s = Simulator(m) creates a simulator object from the mesh object

%Name: init
%Location: <path>/@Simulator
%Purpose: initialize the global matrix from elementary matrix

% modificado em 13/01/2006
% revisado   em 18/04/2006

[uc vc pc cc]=getBC(s.m);

IEN = getIEN(s.m);
X= getX(s.m);
Y=getY(s.m);
Z=getZ(s.m);
B=s.m.B;

nele=size(IEN,1);
nnodes=size(X,1);
nvert=nnodes-nele;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% alocacao de memoria para matrizes e vetores                   %
% K, M, G e D montados para u e v                               %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
K11 = sparse(nnodes,nnodes);
K12 = sparse(nnodes,nnodes);
K22 = sparse(nnodes,nnodes);
Kc = sparse(nvert,nvert);
M = sparse(nnodes,nnodes);
Mz = sparse(nnodes,nnodes);
Mc = sparse(nvert,nvert);
G1 = sparse(nnodes,nvert);
G2 = sparse(nnodes,nvert);

%%% Galerkin
% dx = sparse(nnodes,nvert);
% dy = sparse(nnodes,nvert);
% dxc = sparse(nvert,nvert);
% dyc = sparse(nvert,nvert);


s.uant=sparse(nnodes*2+nvert,1);

%%% determinacao do tipo de elemento na classe TElement
element=FEMMiniElement2d();
elementc=FEMLinElement2d();

%%% parametros obtidos a partir de dados experimentais
s.muzero=2.255;
s.eme=0.81315;
s.cs=cc;

for mele = 1:nele

    v1=IEN(mele,1);
    v2=IEN(mele,2);
    v3=IEN(mele,3);
    c=(s.cs(v1)+s.cs(v2)+s.cs(v3))/3;
    b=(B(v1)+B(v2)+B(v3))/3;
    mu=s.muzero*exp(s.eme*c);
    mub=mu*b;
    Difb=1/mub;
    mele/nele

    [massele,kxx,kyy,kxy,kyx,gxele,gyele,ngleu,nglep,v]=getmgq(element,mele,IEN,X,Y,Z);
    
    vp = v(1:3);
    %ngleu: numero de graus de liberdade do elemento associados a velocidade
    %nglep: idem associados a p

    % K= [11 12]
    % M= [21 22]

    % bloco 11
    K11(v,v) = K11(v,v)+mub*(2*kxx+kyy);
    M(v,v) = M(v,v)+massele;
    % bloco 12
    K12(v,v) = K12(v,v)+mub*kyx;
    % bloco 22
    K22(v,v) = K22(v,v)+mub*(2*kyy+kxx);

    G1(v,vp)=G1(v,vp)+gxele;
    G2(v,vp)=G2(v,vp)+gyele;
    
    %%% Galerkin
    % dx(v,vp) = dx(v,vp)+gxele;
    % dy(v,vp) = dy(v,vp)+gyele;
    
    [masselec,kxxc,kyyc,gxelec,gyelec,nglec,v]=getmgq(elementc,mele,IEN,X,Y,Z);

    %ngleu: numero de graus de liberdade do elemento associados a velocidade
    %nglep: idem associados a p
    Kc(vp,vp) = Kc(vp,vp)+Difb*(kxxc+kyyc);
    Mc(vp,vp) = Mc(vp,vp)+masselec;    
    
    %%% Galerkin
    % dxc(vp,vp) = dxc(vp,vp)+gxelec;  
    % dyc(vp,vp) = dyc(vp,vp)+gyelec; 

end;

s.nnodes=nnodes;
s.nvert=nvert;
s.M=[ M Mz;Mz M ];
s.Mc=Mc;
s.K=[ K11 K12; K12' K22 ];
s.Kc=Kc;
s.G=[ G1;G2 ];
s.D=s.G';

%%% Galerkin
% s.dx=dx;
% s.dy=dy;
% s.dxc=dxc;
% s.dyc=dyc;

s.us=uc;
s.vs=vc;
s.ps=pc;

s.convlin=sparse(nvert,nvert);

