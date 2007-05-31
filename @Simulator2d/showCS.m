function n = showCS(s)

%modelo de documentacao a partir de:
%http://www.engin.umd.umich.edu/CIS/course.des/cis400/matlab/oop.html

%SIMULATOR simulator class constructor.
%   s = Simulator(m) creates a simulator object from the mesh object 

%Name: show
%Location: <path>/@Simulator
%Purpose: show the velocity, pressure and scalar fields                         

% modificado em 28/01/2007
% revisado   em 09/04/2007

IEN = getIEN(s.m);
X= getX(s.m);
Y=getY(s.m);
Z=getZ(s.m);
nvert=size(X,1)-size(IEN,1);

trisurf(IEN(:,1:3),X(1:nvert),Y(1:nvert),s.cs(1:nvert))
shading interp
view(2)

drawnow