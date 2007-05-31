function [A b] = setCoupledBC(s,A,b)

%modelo de documentacao a partir de:
%http://www.engin.umd.umich.edu/CIS/course.des/cis400/matlab/oop.html

%SIMULATOR simulator class constructor.
%   s = Simulator(m) creates a simulator object from the mesh object 

%Name: setCoupledBC
%Location: <path>/@Simulator
%Purpose: set coupled boundary condiction                                

% modificado em 13/04/2006
% revisado   em 09/04/2007

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% ESTRATEGIA DE APLICACAO DA CONDICAO DE CONTORNO                    %
%	1-localizar o vertice onde ha condicao de contorno               % 
%	2-somar a coluna onde o vertice se encontra no vetor do lado     %
%     direto                                                         %
%	3-zerar a coluna e a linha colocando 1 no vertice da matriz      %
%	4-repetir a operacao para outros vertices condicao de contorno   %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% nbc -> numero de condicao de contorno para u
% idbcu -> indice condicao de contorno de u
nbc=size(s.m.idbcu,2);
for ii=1:nbc
    i=s.m.idbcu(ii);
    b=b-s.m.uc(i)*A(:,i);
    A(:,i)=0*b;
    A(i,:)=0*b';
    A(i,i)=1;
    b(i)=s.m.uc(i);
end;

% nbc -> numero de condicao de contorno para v
% idbcv -> indice condicao de contorno de v
nbc=size(s.m.idbcv,2);
for ii=1:nbc
    i=s.m.idbcv(ii);
    b=b-s.m.vc(i)*A(:,i+s.nnodes);
    A(:,i+s.nnodes)=0*b;
    A(i+s.nnodes,:)=0*b';
    A(i+s.nnodes,i+s.nnodes)=1;
    b(i+s.nnodes)=s.m.vc(i);
end;


% nbc -> numero de condicao de contorno para p
% idbcp -> indice condicao de contorno de p
nbc=size(s.m.idbcp,2);
for ii=1:nbc
    i=s.m.idbcp(ii);
    b=b-s.m.pc(i)*A(:,i+2*s.nnodes);
    A(:,i+2*s.nnodes)=0*b;
    A(i+2*s.nnodes,:)=0*b';
    A(i+2*s.nnodes,i+2*s.nnodes)=-1;
    b(i+2*s.nnodes)=-s.m.pc(i);
end;
    
