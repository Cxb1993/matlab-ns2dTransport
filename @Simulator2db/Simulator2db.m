function t = Simulator2db(a)

%modelo de documentacao a partir de:
%http://www.engin.umd.umich.edu/CIS/course.des/cis400/matlab/oop.html

%SIMULATOR simulator class constructor.
%   s = Simulator(m) creates a simulator object from the mesh object 

%Name: Simulator
%Location: <path>/@Simulator
%Purpose: create the mother class Simulator, initializing the matrix

% modificado em 28/01/2007
% revisado   em 09/04/2007

t.us=[];
t.vs=[];
t.ps=[];
t.cs=[];
t.K=[];
t.M=[];
t.inva=[];
t.G=[];
t.D=[];
t.uant=[];
t.dx=[];
t.dy=[];
t.dxc=[];
t.dyc=[];
t.nnodes=[];
t.nvert=[];
t.At=[];
t.Atc=[];
t.Gt=[];
t.Dt=[];
t.Et=[];
t.b1=[];
t.b1c=[];
t.b2=[];
t.ip=[];
t.ipc=[];
t.R1=[];
t.R2=[];
t.Rc=[];
t.r1=[];
t.r2=[];
t.rc=[];
t.convlin=[];
t.eme=[];
t.muzero=[];
t.Kc=[];
t.Mc=[];
t.lambda=[];
t.time=0;

if nargin == 1
    if isa(a,'Model2db');
        t.m=a;
    end;
end;

t = class(t,'Simulator2db');
