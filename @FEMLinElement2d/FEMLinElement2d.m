function t = FEMLinElement2d()
%TELEMENT TElement class constructor.
%   t = TElement() creates a element object;

t.kxxc=[];
t.kyyc=[];

t.masselec=[];
    
telement2d=TElement2d(2);
    
t = class(t,'FEMLinElement2d',telement2d);

