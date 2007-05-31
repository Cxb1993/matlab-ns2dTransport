function s = loadDump(s,dir,name,number)
%modelo de documentacao a partir de:
%http://www.engin.umd.umich.edu/CIS/course.des/cis400/matlab/oop.html

%TERRAIN terrain class constructor.
%   t = Terrain(m) creates a terrain object from the matrix m,
%   containing the terrain data

%Name: Terrain
%Location: <path>/@Terrain
%Purpose: create the mother class Terrain                               

% modificado em 01/05/2007

fname=sprintf('%s%s-%.4d.dmp',dir,name,number);
ret=load('-mat',fname,'s');
s=ret.s;

