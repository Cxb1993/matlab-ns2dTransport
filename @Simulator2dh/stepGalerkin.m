function s = stepGalerkin(s,dt,comp,steptype,Re,Sc)

%modelo de documentacao a partir de:
%http://www.engin.umd.umich.edu/CIS/course.des/cis400/matlab/oop.html

%SIMULATOR simulator class constructor.
%   s = Simulator(m) creates a simulator object from the mesh object 

%Name: step
%Location: <path>/@Simulator
%Purpose: this is the main program, 

% modificado em 13/03/2006
% revisado   em 09/04/2007

IEN = getIEN(s.m);
X= getX(s.m);
Y=getY(s.m);
Z=getZ(s.m);
H=s.m.H;

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

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% % 
% Galerkin                                                      % 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% % 


    conv=(diag(sparse(velu))*s.dx+diag(sparse(velv))*s.dy);
    convu=-conv*s.us(1:nvert,1);
    convv=-conv*s.vs(1:nvert,1);
    convu=(convu.*s.m.outflow);
    convv=(convv.*s.m.outflow);
    veta=(((1/dt)*s.M-(1-alpha)*k*s.K)*s.uant(1:nnodes*2,1));
    va=veta+[convu;convv];
    convc=(diag(sparse(velu(1:nvert)))*s.dxc+diag(sparse(velv(1:nvert)))*s.dyc);
    convcc=-convc*s.cs(1:nvert,1);
    convcc=(convcc.*s.m.outflow(1:nvert));
    vetc=(((1/dt)*s.Mc-(1-alpha)*kc*s.Kc)*s.cs);
    vc=vetc+[convcc];


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
        fman=0;
        lambda=spdiags([H;H],0,2*nnodes,2*nnodes);
        %lambda1=spdiags([H.^(-1);H.^(-1)],0,2*nnodes,2*nnodes);
        
        s=assembleK(s);
        mat=1/dt*s.M+fman*s.M+alpha*k*s.K;
        Mclump=diag(sparse(sum(s.Mc,2)));
        matc=1/dt*Mclump+alpha*kc*s.Kc;    
        
        [At G D E b1 b2 ip] = setUncoupledBC(s,mat,H,va*0);
        [Atc b1c ipc] = setUncoupledCBC(s,matc,vc*0);
        inva=diag(sparse(1./sum(mat,2)));

        % uzawa
        %inva=inv(At);

        E=E-D*lambda*inva*G;
        r1 = symrcm(-At);
        r2 = symrcm(-E);
        rc = symrcm(-Atc);
        At=At(r1,r1);
        E= E(r2,r2);
        Atc=Atc(rc,rc);
        R1 = cholinc(-At,1e-4);
        R2 = cholinc(-E,1e-5);
        Rc = cholinc(-Atc,1e-4);
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
        s.R1=R1;
        s.R2=R2;
        s.Rc=Rc;
        s.r1=r1;
        s.r2=r2;
        s.rc=rc;
        s.inva=inva;
        s.lambda=lambda;
    else
        inva=s.inva;
        lambda=s.lambda;
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
        R1=s.R1;
        R2=s.R2;
        Rc=s.Rc;
        r1=s.r1;
        r2=s.r2;
        rc=s.rc;
    end;

    b1=b1+va.*ip;
    b1c=b1c+vc.*ipc;

    %ut=At\b1;
    ut(r1,1) = pcg(At,b1(r1),1e-6,20,R1',R1);
    %ut = pcg(At,b1,1e-6,10,R1',R1);dd

    cs(rc,1) = pcg(Atc,b1c(rc),1e-6,20,Rc',Rc);
    %cs = pcg(Atc,b1c,1e-6,10,Rc',Rc);

    % [A A G][]=[b1]
    % [A A G][]=[b1]
    % [D D E][]=[b2]

    %E=E-D*inva*G;
    b2=b2-D*lambda*ut;


    %pt=E\b2;
    %pt = pcg(-E,-b2,1e-6,20,R2',R2);
    pt(r2,1) = pcg(-E,-b2(r2),1e-6,20,R2',R2);

    ua=ut-inva*G*pt;

    u=[ua;pt];
    us=u(1:nnodes,1);
    vs=u(1+nnodes:nnodes*2,1);
    ps=u(1+2*nnodes:2*nnodes+nvert,1);
    s.uant=u;
end;

s.us=us;
s.vs=vs;

%s.ps=ps;
% correcao na pressao
% SETUNCOUPLEDBC deve ser ajustado para b2 = 0
s.ps=s.ps+ps;

s.cs=cs;
s.time=s.time+dt;
