clear all
fclose all
%% read file
r=[];
R=[];
rdf=[];
filename='small.txt';
readrdf=fopen(filename,'r');
 while ~feof(readrdf)
     r=fgetl(readrdf);
     R=str2num(r);
     rdf=[rdf;R];
 end
r1=0;
R1=[];
rdf1=[];
k=size(rdf,2);
com=rdf(:,1);
 for i=8:k
%      r1=mean(rdf(:,i));
%      R1=rdf(:,i)-r1;
     R1=rdf(:,i);
     rdf1=[rdf1,R1];
 end
 rdf2=[com,rdf1];
 rdf3=[];
 for i=1:size(rdf2,1)
     if rdf2(i,1)<10000
         continue
     else
         rdf3=[rdf3
               rdf2(i,:)];
     end
 end
 nrdf3=size(rdf3,2);
 for i=1:(nrdf3-1)
     rdf3(:,i+1)=rdf3(:,i+1);
 end
%% corrx

Nevery=1;
Nrepeat=size(rdf3,1);
Nfreq=size(rdf3,1);
kstep=1;
rdf4=[];
rdf4store=[];
aveco2=[2:1:21];
aveco3=[27:1:58];
avecol=[aveco2,aveco3];
navecol=size(avecol,2);
aven=[];
mode=[];
for i=1:navecol
    caln=avecol(i);
    limave=1000;
    dol=mod(Nfreq,limave);
    for n=1:dol
        hlim=n*limave;
        llim=(n-1)*limave+1;
        k=1;
        rdf41=[];
    for p=llim:Nevery:(hlim-2)
        pl=p-llim+1;
        deltat=pl;
        rdf41(k,1)=pl*10;
        ko=0;
        rdf5=[];
        for m=1:(1000-deltat)
             rdf5=[rdf5
                     rdf3(m,1),rdf3(m,caln)*rdf3(m+deltat,caln)];
              ko=ko+1;
        end
    kw=trapz(rdf5(:,1),rdf5(:,2))/ko;
    rdf41(k,2)=kw;
    k=k+1;
    end
    mode=[mode rdf41(:,2)];
    end
    rdf4store(:,i)=mean(mode,2);
end
nrdf4store=size(rdf4store,1);
for i=1:nrdf4store
    rdf4(i,1)=rdf41(i,1);
    rdf4(i,2)=mean(rdf4store(i,:));
end
normalx=rdf4(1,2);
rdf4(:,2)=rdf4(:,2)./normalx;
rdf4(1,2)=rdf4(2,2);
rdf6=[];
for o=1:size(rdf4(:,2))
        if rdf4(o,1)<=100000 && rdf4(o,2)>0.0
            rdf6=[rdf6
                  rdf4(o,:)];
        else
            continue
        end
end
weight1=[];
w=[];
% for i=2:size(rdf6(:,1))
%     weight1=rdf6(2:i,2);
%     vk=var(weight1);
%     w(i)=vk^0.003;
% end
k1=size(rdf6(:,1));
weight1(2:k1)=1./rdf6(2:k1,1).^0.3;
weight1(1)=weight1(2);
w=(weight1./weight1(1))';
am=size(rdf6(:,1));
am1=rdf6((am-200):am,2);
vkm=var(am1);
ft=fittype('(A+vkm).*cos(w.*t).*exp(-t.^(0.24)./alpha)+p3','independent',{'t'},'coefficients',{'A','alpha','w','vkm','p3'});
coe=coeffnames(ft);
sp=[1.23,2.061,-0.005231,2.13,0.3];
options=fitoptions('Method', 'NonlinearLeastSquares','Weight',w,'StartPoint',sp);
fitobject = fit(rdf6(:,1),rdf6(:,2),ft,options);
plot(fitobject,rdf4(:,1),rdf4(:,2))
xlabel('\Delta t')
ylabel('correlation')
axis([0 5*10^4 -0.5 1.0])
figure
plot(rdf2(:,1),rdf2(:,caln))
%% print figure
figHandles = findall(0,'Type','figure');
savefig(figHandles,['test_' datestr(now,30)],'compact');
%% save variable
varifile='test.mat';
save(varifile);