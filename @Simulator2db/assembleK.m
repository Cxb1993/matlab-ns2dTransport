function s = assembleK(s)

%modelo de documentacao a partir de:
%http://www.engin.umd.umich.edu/CIS/course.des/cis400/matlab/oop.html

%SIMULATOR simulator class constructor.
%   s = Simulator(m) creates a simulator object from the mesh object,
%   containing the simulator 01-05-2007

%Name: assembleK
%Location: <path>/@Simulator
%Purpose: matrix K and Kc update

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
% K, M, G e D montados para u, v e w                            %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

K11 = sparse(nnodes,nnodes);
K12 = sparse(nnodes,nnodes);
K22 = sparse(nnodes,nnodes);
Kc = sparse(nvert,nvert);

%%% determinacao do tipo de elemento na classe TElement
element=FEMMiniElement2d();
elementc=FEMLinElement2d();

%%% parametros obtidos a partir de dados experimentais
s.muzero=2.255;
s.eme=0.81315;

for mele = 1:nele
    v1=IEN(mele,1);
    v2=IEN(mele,2);
    v3=IEN(mele,3);
    c=(s.cs(v1)+s.cs(v2)+s.cs(v3))/3;
    b=(B(v1)+B(v2)+B(v3))/3;
    mu=s.muzero*exp(s.eme*c);
    mub=mu*b;
    Difb=1/mub;
    mele/nele;

    [kxx,kyy,kxy,kyx,ngleu,nglep,v]=getkgq(element,mele,IEN,X,Y,Z);
    
    vp = v(1:3);
    %ngleu: numero de graus de liberdade do elemento associados a velocidade
    %nglep: idem associados a p


    % K= [11 12]
    %    [21 22]

    % bloco 11
    K11(v,v) = K11(v,v)+mub*(2*kxx+kyy);

    % bloco 12
    K12(v,v) = K12(v,v)+mub*kxy;

    % bloco 22
    K22(v,v) = K22(v,v)+mub*(kxx+2*kyy);

    [kxxc,kyyc,nglec,v]=getkgq(elementc,mele,IEN,X,Y,Z);

    % bloco 1
    Kc(vp,vp) = Kc(vp,vp)+Difb*(kxxc+kyyc);

end;

s.K=[K11 K12;K12' K22];
s.Kc=Kc;
