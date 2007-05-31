path('/home/gustavo/sandbox/NS2dh',path)
workDir='/gesar/LEN/gustavo/simulacao2dh/';
workDirVtk='/gesar/LEN/gustavo/simulacao2dh/vtk/';

Re=10000;
Sc=20;

m1=Model2dh();
s1=Simulator2dh(m1);

i=1;
s1=loadDump(s1,workDir,'sim',i);
m1=s1.m;
dt=s1.time/i;
cfl=dt/(sqrt((max(m1.Y)-min(m1.Y))*(max(m1.X)-min(m1.X))/s1.nvert)...
    /max(s1.m.uc));

i=20;
s1=loadSol(s1,workDir,'sim',i);
for i=21:800
    s1=step(s1,dt,false,'uncoupled',Re,Sc);
    %saveSol(s1,workDir,'sim',i);
    show(s1)
    %vtkCompleteOut(s1,workDirVtk,'field',i);
    %vort(s1);
end;
