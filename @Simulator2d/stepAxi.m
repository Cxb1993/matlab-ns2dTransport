function s = stepAxi(s,dt,comp,steptype,Re,Sc)

%modelo de documentacao a partir de:
%http://www.engin.umd.umich.edu/CIS/course.des/cis400/matlab/oop.html

%SIMULATOR simulator class constructor.
%   s = Simulator(m) creates a simulator object from the mesh object

%Name: stepAxi
%Location: <path>/@Simulator
%Purpose: this is the main program,

% modificado em 13/03/2006
% revisado   em 09/04/2007

IEN = getIEN(s.m);
X= getX(s.m);
Y=getY(s.m);
Z=getZ(s.m);

nelem=size(IEN,1);
nvert=size(X,1)-nelem;
nnodes=size(X,1);

velu=s.us;
velv=s.vs;
velc=s.cs;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%  alpha = 0   -> explicito                                     %
%  alpha = 0.5 -> crank-nicholson                               %
%  alpha = 1   -> implicito                                     %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

alpha=1;

%% coeficiente do laplaciano - viscosidade e concentracao
k=1/Re;
kc=1/(Re*Sc);

%%% calculo do termo convectivo para o met semi-lagrangeano
[up,vp,cp] = convectLin(s,dt);
Mclump=diag(sparse(sum(s.Mc,2)));

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% Lagrangiano  --  montagem dos vetores e matriz                %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% no pressure correction
%va=((1/dt)*s.M-(1-alpha)*k*s.K)*[velu;velv];
% pressure correction - Lagrangian
%va=((1/dt)*s.M-(1-alpha)*k*s.K)*[velu;velv]-s.G*s.ps; 

%vc=((1/dt)*Mclump-(1-alpha)*kc*s.Kc)*velc;


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% semi-lagrangiano  --  montagem dos vetores e matriz           %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% no pressure correction
%va=((1/dt)*s.M-(1-alpha)*k*s.K)*[up;vp]; 
% pressure correction - Semi-Lagrangian
va=((1/dt)*s.M-(1-alpha)*k*s.K)*[up;vp]-s.G*s.ps;

vc=((1/dt)*Mclump-(1-alpha)*kc*s.Kc)*cp;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% metodo acoplado                                               %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

if(strcmp(steptype,'coupled'))

    mat=1/dt*s.M+alpha*k*s.K;
    % A=[K G]=[Ku  0 Gx][ ]=[us]  [us]=[ ]
    %   [D 0] [ 0 Kv Gy][u]=[vs]  [vs]=[b]
    %         [Dx Dy  0][ ]=[ps]  [ps]=[ ]

    A=sparse([mat s.G; s.D 0*(s.D*s.G)]);
    b=[va;zeros(nvert,1)];
    [A b] = setCoupledBC(s,A,b);
    u=A\b;
    us=u(1:nnodes,1);
    vs=u(1+nnodes:nnodes*2,1);
    ps=u(1+2*nnodes:2*nnodes+nvert,1);
    s.uant=u;

    Mclump=diag(sparse(sum(s.Mc,2)));
    matc=1/dt*Mclump+alpha*kc*s.Kc;
    Ac=matc;
    b1c=vc;
    [Atc b1c] = setCoupledCBC(s,Ac,b1c);
    c=Atc\b1c;
    cs=c;

end;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% Metodo da Projecao discreto baseado em decomposicao LU        %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

if(strcmp(steptype,'uncoupled'))

    if(comp)
        s=assembleK(s);
        mat=1/dt*s.M+alpha*k*s.K;
        Mclump=diag(sparse(sum(s.Mc,2)));
        matc=1/dt*Mclump+alpha*kc*s.Kc;

        [At G D E b1 b2 ip] = setUncoupledBC(s,mat,va*0);
        [Atc b1c ipc] = setUncoupledCBC(s,matc,vc*0);
        inva=diag(sparse(1./sum(mat,2)));

        % uzawa
        %inva=inv(At);

        E=E-D*inva*G;
        s.At=At;
        s.Atc=Atc;
        s.Gt=G;
        s.Dt=D;
        s.Et=E;
        s.b1=b1;
        s.b1c=b1c;
        s.b2=b2;
        s.ip=ip;
        s.ipc=ipc;
        s.inva=inva;
    else
        inva=s.inva;
        At=s.At;
        Atc=s.Atc;
        G=s.Gt;
        D=s.Dt;
        E=s.Et;
        b1=s.b1;
        b1c=s.b1c;
        b2=s.b2;
        ip=s.ip;
        ipc=s.ipc;
    end;

    b1=b1+va.*ip;
    b1c=b1c+vc.*ipc;

    % solution of Concentration equation
    cs=Atc\b1c;
    %cs = pcg(Atc,b1c,1e-6,10,Rc',Rc);

    ut=At\b1;

    % [A A G][]=[b1]
    % [A A G][]=[b1]
    % [D D E][]=[b2]

    %E=E-D*inva*G;
    b2=b2-D*ut;

    pt=E\b2;

    ua=ut-inva*G*pt;

    u=[ua;pt];
    us=u(1:nnodes,1);
    vs=u(1+nnodes:nnodes*2,1);
    ps=u(1+2*nnodes:2*nnodes+nvert,1);
    s.uant=u;
end;

s.us=us;
s.vs=vs;

% no pressure correction
%s.ps=ps;

% pressure correction
% SETUNCOUPLEDBC must be set to b2 = 0
s.ps=s.ps+ps; 

s.cs=cs;
s.time=s.time+dt;

