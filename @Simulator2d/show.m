function n = show(s)

%modelo de documentacao a partir de:
%http://www.engin.umd.umich.edu/CIS/course.des/cis400/matlab/oop.html

%SIMULATOR simulator class constructor.
%   s = Simulator(m) creates a simulator object from the mesh object 

%Name: show
%Location: <path>/@Simulator
%Purpose: show the velocity, pressure and scalar fields                         

% modificado em 28/01/2007
% revisado   em 09/04/2007

uc=getUC(s.m);
IEN = getIEN(s.m);
X= getX(s.m);
Y=getY(s.m);
Z=getZ(s.m);
nvert=size(X,1)-size(IEN,1);

subplot(4,2,1)
trisurf(IEN(:,1:3),X(1:nvert),Y(1:nvert),s.us(1:nvert))
shading interp
view(2)

subplot(4,2,3)
trisurf(IEN(:,1:3),X(1:nvert),Y(1:nvert),s.vs(1:nvert))
shading interp
view(2)

subplot(4,2,5)
trisurf(IEN(:,1:3),X(1:nvert),Y(1:nvert),s.ps(1:nvert))
shading interp
view(2)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% plotar graficos de c e t (escalares)                          %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

subplot(4,2,7)
trisurf(IEN(:,1:3),X(1:nvert),Y(1:nvert),s.cs(1:nvert))
shading interp
view(2)

drawnow

