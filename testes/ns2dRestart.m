path('/home/gustavo/sandbox/NS2d',path)
workDir='/gesar/LEN/gustavo/simulacao2d/';
workDirVtk='/gesar/LEN/gustavo/simulacao2d/vtk/';

Re=10000;
Sc=2000;

m1=Model2d();
s1=Simulator2d(m1);

i=1;
s1=loadDump(s1,workDir,'sim',i);
m1=s1.m;
dt=s1.time/i;
cfl=dt/(sqrt((max(m1.Y)-min(m1.Y))*(max(m1.X)-min(m1.X))/s1.nvert)...
    /max(s1.m.uc));


i=20;
s1=loadSol(s1,workDir,'sim',i);
for i=21:800
    i
    s1=step(s1,dt,false,'uncoupled',Re,Sc);
    %saveSol(s1,workDir,'sim',i);
    show(s1)
    %vtkCompleteOut(s1,workDirVtk,'field',i);
    %vort(s1);
end;
