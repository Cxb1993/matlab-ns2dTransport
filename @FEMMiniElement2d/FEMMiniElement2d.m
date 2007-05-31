function t = FEMMiniElement2d()
%TELEMENT2d TElement2d class constructor.
%   t = TElement() creates a element object;

t.kxx=[];
t.kxy=[];
t.kyx=[];
t.kyy=[];

t.massele=[];
    
t.dxele=[];
t.dyele=[];
     
telement2d=TElement2d(2);
    
t = class(t,'FEMMiniElement2d',telement2d);

