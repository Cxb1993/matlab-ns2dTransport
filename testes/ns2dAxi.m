path('~/projects/matlab/ns2dTransport/',path)
workDir='~/projects/matlab/';
workDirVtk='~/projects/matlab/';

Re=100;
Sc=10;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% utilizacao:                                               %             
% test: test (step,cavity,couette)                          %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

m1=Model2d();
m1=test(m1,29,15,'poiseuilleAxi');


%modelOut(m1,workDirVtk,'mesh');
%show(m1);

s1=Simulator2d(m1);
s1=initAxi(s1);

cfl=1;
dt=cfl*sqrt((max(m1.Y)-min(m1.Y))*(max(m1.X)-min(m1.X))/s1.nvert)/max(s1.us);

i=1
s1=stepAxi(s1,dt,true,'uncoupled',Re,Sc);
%saveDump(s1,workDir,'sim',i)
%saveSol(s1,workDir,'sim',i)

for i=2:40
    i
    s1=stepAxi(s1,dt,false,'uncoupled',Re,Sc);
    %saveSol(s1,workDir,'sim',i)
    show(s1)
    vtkCompleteOut(s1,workDirVtk,'field',i)
    %vort(s1);
end;
