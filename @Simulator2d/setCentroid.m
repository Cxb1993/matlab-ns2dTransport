function ut = setCentroid(s,ut)

%modelo de documentacao a partir de:
%http://www.engin.umd.umich.edu/CIS/course.des/cis400/matlab/oop.html

%SIMULATOR simulator class constructor.
%   s = Simulator(m) creates a simulator object from the mesh object

%Name: setCentroid
%Location: <path>/@Simulator
%Purpose: define mini element centroid 

% modificado em 05/04/2007
% revisado   em 09/04/2007

IEN = getIEN(s.m);
X= getX(s.m);
Y=getY(s.m);
Z=getZ(s.m);
nele=size(IEN,1);
nnodes=size(X,1);


for mele = 1:nele

    v1=IEN(mele,1);
    v2=IEN(mele,2);
    v3=IEN(mele,3);
    v4=IEN(mele,4);
    v=[v1;v2;v3;v4];
   
    um=0;
    vm=0;
    for i=1:3
		um=um+ut(v(i));
		vm= vm+ut(v(i)+nnodes);
    end;
    ut(v4)=um/3;
    ut(v4+nnodes)=vm/3;
end;
