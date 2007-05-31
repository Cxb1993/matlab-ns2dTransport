function m = test(m,a,a1,problem)

%modelo de documentacao a partir de:
%http://www.engin.umd.umich.edu/CIS/course.des/cis400/matlab/oop.html

%MODEL model class constructor.
%   m = Model() creates a mesho object

%Name: test
%Location: <path>/@Model
%Purpose: model method to run test problems

% modificado em 01/05/2007
% revisado   em 09/04/2007

m.IEN=[];

[X,Y] = meshgrid(1:a,1:a1);

m.X = reshape(X,1,[])';
m.Y = reshape(Y,1,[])';
m.Z = 0*m.X;

m.IEN = delaunay(m.X,m.Y,{'Qt','QbB','Qc'});

nvert=size(m.X,1);
nelem=size(m.IEN,1);
nnodes=nvert+nelem;
IEN = zeros(nelem,4);
IEN(:,1:3)=m.IEN;

X=zeros(nnodes,1);
X(1:nvert)=m.X;

Y=zeros(nnodes,1);
Y(1:nvert)=m.Y;

Z=zeros(nnodes,1);
Z(1:nvert)=m.Z;

%%%%% verificacao do sentido de numeracao do triangulo %%%%%

for i=1:nelem
    v1=IEN(i,1);
    v2=IEN(i,2);
    v3=IEN(i,3);


    A=0.5*det ([X(v2)-X(v1) X(v3)-X(v1);Y(v2)-Y(v1) Y(v3)-Y(v1)]);
    if(abs(A)<1e-10)

        disp('triangulo singular');
    end
    if(A<0)

        % disp('triangulo negativo');

        IEN(i,2)=v3;
        IEN(i,3)=v2;
    end;

    v4=nvert+i;
    IEN(i,4)=v4;
    X(v4)=(X(v1)+X(v2)+X(v3))/3;
    Y(v4)=(Y(v1)+Y(v2)+Y(v3))/3;
    Z(v4)=(Z(v1)+Z(v2)+Z(v3))/3;

end;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%         condicoes de contorno e condicoes inicias             %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

uc=sparse(nnodes,1);
vc=sparse(nnodes,1);
pc=sparse(nvert,1);
cc=sparse(nvert,1);
idbcu=[];
idbcv=[];
idbcp=[];
idbcc=[];

if(strcmp(problem,'cavity'))

    for i=1:nnodes
        if((X(i)==1)||(X(i)==a)||(Y(i)==1)||(Y(i)==a1))
            idbcu=[idbcu i];
            idbcv=[idbcv i];
            uc(i)=0;
            vc(i)=0;
            if(Y(i)==a1)
                uc(i)=1;
            end;
        end;
    end;

    for i=1:nvert
        if((X(i)==a)&&(Y(i)==1))
            idbcp=[idbcp i];
            pc(i)=0;
            outflow(i)=0;
        end;
        if(Y(i)==a1)
            idbcc=[idbcc i];
            cc(i)=1;
        end;
    end;
end;

if(strcmp(problem,'step'))
    outflow=ones(size(X));
    j=0;
    for i=1:nnodes
        if((X(i)==1)||(Y(i)==1)||(Y(i)==a1))
            idbcu=[idbcu i];
            idbcv=[idbcv i];

            uc(i)=0;
            vc(i)=0;
            if((X(i)==1)&&(Y(i)>a1/2)&&(Y(i)<a1))
                uc(i)=1;
            end;
        end;
    end;

    for i=1:nvert

        if(X(i)==a)
            idbcp=[idbcp i];
            pc(i)=0;
            outflow(i)=0;
        end;

        if((X(i)==1)||(Y(i)==1)||(Y(i)==a1))
            idbcc=[idbcc i];
            if((X(i)==1)&&(Y(i)>a1/2)&&(Y(i)<a1))
                cc(i)=1;
            end;
            if((X(i)==1)&&(Y(i)<=a1/2))
                cc(i)=0;
            end;
        end;
    end;
end;

if(strcmp(problem,'symmetric'))
    outflow=ones(size(X));
    j=0;
    for i=1:nnodes
        if((X(i)==1)||(Y(i)==1))
            idbcu=[idbcu i];
            idbcv=[idbcv i];

            vc(i)=0;
            if((X(i)==1)&&(Y(i)>a1/2))
                uc(i)=1;

            end;
            if((X(i)==1)&&(Y(i)<=a1/2))
                uc(i)=0;

            end;
        end;
        if(Y(i)==a1)
            %idbcu=[idbcu i];
            idbcv=[idbcv i];
            vc(i)=0;
        end;
    end;

    j=0;
    for i=1:nvert

        if(X(i)==a)
            idbcp=[idbcp i];
            pc(i)=0;
            outflow(i)=0;
        end;


        if((X(i)==1)||(Y(i)==1)||(Y(i)==a1))
            idbcc=[idbcc i];
            if((X(i)==1)&&(Y(i)>a1/2)&&(Y(i)<a1))
                cc(i)=1;
            end;
            if((X(i)==1)&&(Y(i)<=a1/2))
                cc(i)=0;
            end;
        end;
    end;
end;

if(strcmp(problem,'couette'))
    outflow=ones(size(X));

    for i=1:nnodes
        if((X(i)==1)||(Y(i)==1)||(Y(i)==a1))
            idbcu=[idbcu i];
            idbcv=[idbcv i];

            uc(i)=0;
            vc(i)=0;
            if((X(i)==1) &&(Y(i)>1)&&(Y(i)<a1))
                uc(i)=1;
            end;
        end;
    end;

    for i=1:nvert

        if(X(i)==a)
            idbcp=[idbcp i];
            pc(i)=0;
            outflow(i)=0;
        end;

        if((X(i)==1) &&(Y(i)>1)&&(Y(i)<a1))
            idbcc=[idbcc i];
            cc(i)=1;
        end;

    end;
end;

m.idbcu=idbcu;
m.idbcv=idbcv;
m.idbcp=idbcp;
m.idbcc=idbcc;

m.uc=uc;
m.vc=vc;
m.pc=pc;
m.cc=cc;
m.outflow=outflow;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% adimensionalizacao da malha                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

factor=1/(max(Y)-min(Y));
Y=(Y-min(Y))*factor;
X=(X-min(X))*factor;

% profundidade adimensional
for i=1:nnodes
    if X(i)>1
        B(i)=2;
    else
        B(i)=1;
    end;
end;
% for i=1:nnodes
%     if X(i)<1.5
%         if X(i)<1
%             B(i)=2;
%         else
%         B(i)=3;
%         end;
%     else
%         B(i)=1;
%     end;
% end;
B=B';

% B=2+0.1+(X.*(X-2))+4*(Y.*(Y-1));
% B=ones(size(X));

% L=max(X)-min(X);
% Y=Y-0.5;
% 
% c1=0.2; % sinuoso
% c2=-0.2; % varicoso
% c3=0.05;
% 
% R=X+0.25;
% T=Y*pi/4+sin(2*pi*X/L).*(c1+c2*Y.*X)+cos(3*pi*X/L)*c3;;
% 
% Y=R.*sin(T);
% X=R.*cos(T);

m.IEN=IEN;
m.X=X;
m.Y=Y;
m.Z=Z;
m.B=B;

