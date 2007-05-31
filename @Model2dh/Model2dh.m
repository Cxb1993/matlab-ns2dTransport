function m = Model2dh()
%modelo de documentacao a partir de:
%http://www.engin.umd.umich.edu/CIS/course.des/cis400/matlab/oop.html

%MODEL model class constructor.
%   m = Model2dh() creates a mesh object 

%Name: Model
%Location: <path>/@Model
%Purpose: create the mother class Model

% modificado em 05/04/2007
% revisado   em 09/04/2007

m.X=[];
m.Y=[];
m.Z=[];
m.H=[];
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

m = class(m,'Model2dh');

if nargin ~= 0
 	disp('classe Model does not need arguments');
end;
