function n = vtkVOut(s,dir,name,number)
%modelo de documentacao a partir de:
%http://www.engin.umd.umich.edu/CIS/course.des/cis400/matlab/oop.html

%SIMULATOR simulator class constructor.
%   s = Simulator(m) creates a simulator object from the mesh object

%Name: vtkVOut
%Location: <path>/@Simulator
%Purpose: save velocity V solution file in VTK format                          

% modificado em 20/05/2007
% revisado   em 20/05/2007

IEN = getIEN(s.m);
X= getX(s.m);
Y=getY(s.m);
Z=getZ(s.m);
ps=s.ps;
us=s.us;

nvert=size(ps,1);
nelem=size(IEN,1);
nnodes=nvert+nelem;

for i=1:size(IEN,1)
	for j=1:3
		IEN(i,j)=IEN(i,j)-1;
	end;
end;

fname = sprintf('%s%s-%.4d.vtk',dir,name,number);
fid = fopen(fname, 'wt');

fprintf(fid, '# vtk DataFile Version 1.0\n');
fprintf(fid, 'Velocity V 2dh\n');
fprintf(fid, 'ASCII\n');
fprintf(fid, 'DATASET UNSTRUCTURED_GRID\n');
fprintf(fid, 'POINTS %d float\n',nvert);

for k=1:nvert
    fprintf(fid, '%2.5f %2.5f %2.2f\n', X(k), Y(k), 0);
end;

fprintf(fid, '\n',i);
fprintf(fid, 'CELLS %d %d\n',nelem,4*nelem);
for i=1:nelem
    fprintf(fid, '3 ');
    for j=1:3
        fprintf(fid, '%d ', IEN(i,j));
    end;
    fprintf(fid, '\n',i);
end;

fprintf(fid, '\n',i);
fprintf(fid, 'CELL_TYPES %d\n',nelem);
for i=1:nelem
    fprintf(fid, '5 ');
end;
fprintf(fid, '\n',i);

fprintf(fid, '\n',i);
fprintf(fid, 'POINT_DATA %d\n',nvert);
fprintf(fid, 'SCALARS scalars float\n');
fprintf(fid, 'LOOKUP_TABLE default\n');
for i=1:nvert
    fprintf(fid, '%2.10f\n', vs(i));
end;

fclose(fid);

