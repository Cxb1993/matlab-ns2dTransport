function [At G D E b1 b2 ip] = setUncoupledBC(s,At,B,b1)

%modelo de documentacao a partir de:
%http://www.engin.umd.umich.edu/CIS/course.des/cis400/matlab/oop.html

%SIMULATOR simulator class constructor.
%   s = Simulator(m) creates a simulator object from the mesh object 

%Name: setUncoupledBC
%Location: <path>/@Simulator
%Purpose: set uncoupled boundary condiction                                

% modificado em 25/04/2006
% revisado   em 09/04/2007

G=s.G;
D=s.D;
b2=sparse(s.nvert,1);
ip=sparse(ones(s.nnodes*2,1));
% nbc -> numero de condicao de contorno para u
% idbcu -> indice condicao de contorno de u
nbc=size(s.m.idbcu,2);
for ii=1:nbc
    i=s.m.idbcu(ii);
    b1=b1-s.m.uc(i)*At(:,i);
    b2=b2-s.m.uc(i)*B(i)*D(:,i);
    At(:,i)=0*b1;
    D(:,i)=0*b2;
    At(i,:)=0*b1';
    G(i,:)=0*b2';
    At(i,i)=1;
    b1(i)=s.m.uc(i);
    ip(i)=0;
end;

% nbc -> numero de condicao de contorno para v
% idbcv -> indice condicao de contorno de v
nbc=size(s.m.idbcv,2);
for ii=1:nbc
    i=s.m.idbcv(ii);
    b1=b1-s.m.vc(i)*At(:,i+s.nnodes);
    b2=b2-s.m.vc(i)*B(i)*D(:,i+s.nnodes);
    At(:,i+s.nnodes)=0*b1;
    D(:,i+s.nnodes)=0*b2;
    At(i+s.nnodes,:)=0*b1';
    G(i+s.nnodes,:)=0*b2';
    At(i+s.nnodes,i+s.nnodes)=1;
    b1(i+s.nnodes)=s.m.vc(i);
    ip(i+s.nnodes)=0;
end;


% nbc -> numero de condicao de contorno para p
% idbcp -> indice condicao de contorno de p
nbc=size(s.m.idbcp,2);
E=sparse(s.nvert,s.nvert);

for ii=1:nbc
    i=s.m.idbcp(ii);
    b1=b1-s.m.pc(i)*G(:,i);
    G(:,i)=0*b1;
    D(i,:)=0*b1';
    E(i,i)=-1;
    b2(i)=-s.m.pc(i);
end;

