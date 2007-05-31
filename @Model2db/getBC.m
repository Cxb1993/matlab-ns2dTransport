function [uc vc pc cc]=getBC(m)

%modelo de documentacao a partir de:
%http://www.engin.umd.umich.edu/CIS/course.des/cis400/matlab/oop.html

%MODEL model class constructor.
%   m = Model() creates a mesh object

%Name: getBC
%Location: <path>/@Model
%Purpose: get all boundary condiction                                

% modificado em 13/01/2006
% revisado   em 09/04/2007

uc=m.uc;
vc=m.vc;
pc=m.pc;
cc=m.cc;

