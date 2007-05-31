function n = show(m)

%modelo de documentacao a partir de:
%http://www.engin.umd.umich.edu/CIS/course.des/cis400/matlab/oop.html

%MODEL model class constructor.
%   m = Model() creates a mesho object

%Name: show
%Location: <path>/@Model
%Purpose: plot the mesh object                               

% modificado em 13/01/2006
% revisado   em 09/04/2007

IEN=m.IEN(:,1:3);
nele=size(IEN,1);
nvert=size(m.X,1);
X=m.X(1:nvert);
Y=m.Y(1:nvert);
Z=m.Z(1:nvert);

trimesh(IEN,X,Y,Z);

n=0;


