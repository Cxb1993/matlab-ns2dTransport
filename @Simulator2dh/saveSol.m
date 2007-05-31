function n = saveSol(s,dir,name,number)

%modelo de documentacao a partir de:
%http://www.engin.umd.umich.edu/CIS/course.des/cis400/matlab/oop.html

%SIMULATOR simulator class constructor.
%   s = Simulator(m) creates a simulatort object from the mesh object

%Name: saveSol
%Location: <path>/@Simulator
%Purpose: save solution file                          

% modificado em 26/11/2006
% revisado   em 09/04/2007

fname=sprintf('%s%s-%.4f.sol',dir,name,number);

us=s.us;
vs=s.vs;
ps=s.ps;
cs=s.cs;
K=s.K;
Kc=s.Kc;

save(fname,'us','vs','ps','cs','K','Kc');

