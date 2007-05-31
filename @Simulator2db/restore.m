function t = restore(a,a1,a2)
%modelo de documentacao a partir de:
%http://www.engin.umd.umich.edu/CIS/course.des/cis400/matlab/oop.html

%TERRAIN terrain class constructor.
%   t = Terrain(m) creates a terrain object from the matrix m,
%   containing the terrain data

%Name: Terrain
%Location: <path>/@Terrain
%Purpose: create the mother class Terrain                               

%para qualquer numero de argumentos de entrada menor que 3 sao definidos:
%a origem do sistema de coordenadas
%o tipo de sistema de coordenadas
%um objeto t da classe Terrain
% modificado em 13/01/2006

if nargin == 0
    t.datum='Corrego Alegre';
    t.coorSys = 'metri';
    t.projType = 'UTM';
    t.projCat='cylindric';
    
    t = class(t,'Terrain');
elseif isa(a,'Terrain')%se o argumento de entrada for um objeto de Terrain entao a variavel t recebe este argumento
    t = a;
else
    if nargin == 3
        t.datum='Corrego Alegre';
    t.coorSys = 'metri';
    t.projType = 'UTM';
    t.projCat='cylindric';
    
        t = class(t,'Terrain');
    else
        t.datum='Corrego Alegre';
    t.coorSys = 'metri';
    t.projType = 'UTM';
    t.projCat='cylindric';
    
        t = class(t,'Terrain');
    end;
end