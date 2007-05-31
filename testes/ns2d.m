path('/home/gustavo/sandbox/NS2d',path)
workDir='/gesar/LEN/gustavo/simulacao2d/';
workDirVtk='/gesar/LEN/gustavo/simulacao2d/vtk/';

Re=10000;
Sc=2000;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% utilizacao:                                               %             
% test: test (step,cavity,couette)                          %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

m1=Model2d();
m1=test(m1,29,15,'step');

%show(m1);

s1=Simulator2d(m1)
s1=init(s1);

cfl=1;
dt=cfl*sqrt((max(m1.Y)-min(m1.Y))*(max(m1.X)-min(m1.X))/s1.nvert)/max(s1.us);

i=1
s1=step(s1,dt,true,'uncoupled',Re,Sc);
saveDump(s1,workDir,'sim',i)
saveSol(s1,workDir,'sim',i)

for i=2:40
    i
    s1=step(s1,dt,false,'uncoupled',Re,Sc);
    %saveSol(s1,workDir,'sim',i)
    show(s1)
    %vtkCompleteOut(s1,workDirVtk,'field',i)
    %vort(s1);
end;