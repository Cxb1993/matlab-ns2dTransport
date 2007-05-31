function [Atc b1c ipc] = setUncoupledCBC(s,Atc,b1c)
%modelo de documentacao a partir de:
%http://www.engin.umd.umich.edu/CIS/course.des/cis400/matlab/oop.html

%SIMULATOR simulator class constructor.
%   s = Simulator(m) creates a simulator object from the mesh object 

%Name: setUncoupledCBC
%Location: <path>/@Simulator
%Purpose: set uncoupled boundary condiction of mass transport                                
% modificado em 01/05/2006
% revisado   em 09/04/2007

ipc=sparse(ones(s.nvert,1));
% nbcc -> numero de condicao de contorno para c
% idbcc -> indice condicao de contorno de c
nbc=size(s.m.idbcc,2);
for ii=1:nbc
	i=s.m.idbcc(ii);
	b1c=b1c-s.m.cc(i)*Atc(:,i);
	Atc(:,i)=0*b1c;
	Atc(i,:)=0*b1c';
	Atc(i,i)=1;
	b1c(i)=s.m.cc(i);
	ipc(i)=0;
end;
