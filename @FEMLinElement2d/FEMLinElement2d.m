function t = FEMLinElement2d()
%TELEMENT TElement class constructor.
%   t = TElement() creates a element object;

t.masselec=[];
t.kxxc=[];
t.kyyc=[];
t.kyc=[];
t.kxc=[];
t.gxelec=[];
t.gyelec=[];

telement2d=TElement2d(2);
    
t = class(t,'FEMLinElement2d',telement2d);

