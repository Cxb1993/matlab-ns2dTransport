function s = vort(s)
%modelo de documentacao a partir de:
%http://www.engin.umd.umich.edu/CIS/course.des/cis400/matlab/oop.html

%SIMULATOR terrain class constructor.
%   s = vort(m) creates a simulator object from the mesho object

%Name: vort
%Location: <path>/@Simulator
%Purpose: plot vorticy field (curl)                               

% modificado em 01/05/2007

%%% VORTICIDADE : rot = dv/dx - du/dy

[uc vc pc cc]=getBC(s.m);
IEN = getIEN(s.m);
X= getX(s.m);
Y=getY(s.m);
Z=getZ(s.m);

nelem=size(IEN,1);
nvert=size(X,1)-nelem;
nnodes=size(X,1);

mat=s.M;

velu=s.us;
velv=s.vs;


IntVort=s.dx*s.vs(1:nvert,1)-s.dy*s.us(1:nvert,1) ;
A=mat(1:nnodes,1:nnodes);

%%%%% esta eh a matriz A
vort=A\IntVort;
vort=vort.*s.m.outflow;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% plotar grafico de vorticidade                                 %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

subplot(4,2,4)
trisurf(IEN(:,1:3),X(1:nvert),Y(1:nvert),vort(1:nvert))
shading interp
view(2)

