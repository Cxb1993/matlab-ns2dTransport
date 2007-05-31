function m = Model2d()
%modelo de documentacao a partir de:
%http://www.engin.umd.umich.edu/CIS/course.des/cis400/matlab/oop.html

%MODEL model class constructor.
%   m = Model() creates a mesh object 

%Name: Model2d
%Location: <path>/@Model2d
%Purpose: create the mother class Model

% modificado em 05/04/2007
% revisado   em 09/04/2007

m.X=[];
m.Y=[];
m.Z=[];
m.idbcu=[];
m.idbcv=[];
m.idbcp=[];
m.idbcc=[];
m.uc=[];
m.vc=[];
m.pc=[];
m.cc=[];
m.IEN=[];
m.outflow=[];

m = class(m,'Model2d');

if nargin ~= 0
 	disp('classe Model does not need arguments');
end;
