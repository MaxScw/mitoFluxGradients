clear all;
clc;
close all;

x = 36:44;
y = 16:24;

param(1,:)= [40,20,1.5,1.5,0,40,5];
%param(2,:) = [40,26,1.5,1.5,0,40,10];
c = [100,0,1];
%param = [40,20,2,2,0,40,10,40,28,2,2,0,30,10]';

GlobalVarTracking;
ConstructGlobalVar(x,y,c);
z = JointGaussianModel(x,y,param,c);

h1 = figure;
imagesc(z);
axis image;


Nsim = 500;
pfit = zeros(Nsim,size(param,2));
Ntot = zeros(Nsim,1);
cm = zeros(Nsim,2);
xrms = zeros(Nsim,1);
%%
%matlabpool open

tic
for i = 1:Nsim
    
    simz = SimulateData(x,y,param,[100,0,1]);
    Ntot(i) = sum(simz(:));
%     
%     h3 = figure;
%     imagesc(simz)
    
    
    pmin = param;
    pmax = param;
    pmin(1,1:2) = pmin(1:2)-1;
    pmax(1,1:2) = pmax(1:2)+1;
    pmin(1,3:4) = pmin(1,3:4)-1;
    pmax(1,3:4) = pmax(1,3:4)+1;
    pmin(1,6) = max(simz(:))-10;
    pmax(6) = max(simz(:))+10;
    pmin(7)= 0;
    pmax(7)= 15;
    pinit = pmin+rand(1,7).*(pmax-pmin);
    [X,Y]=meshgrid(x,y);
    cmx = sum(X(:).*simz(:))/sum(simz(:));
    cmy = sum(Y(:).*simz(:))/sum(simz(:));
    cm(i,:) = [cmx,cmy];
    pinit(1:2) = cm(i,:); 
    pinit(6) = max(simz(:));
    pinit(5) = 0;
    dp = ones(1,7)*0.01;
    dp(4) = 0;
    dp(5) = 0;
    
    nonzero_z = simz;
    nonzero_z(simz==0)=1;
    sigz = sqrt(nonzero_z);
    weight = 1./sigz;
    [pfit(i,:),X2,sigp,sigy,corr,Rsq,cvg_hst, converged] = lm2(@JointGaussianModel,pinit,x,y,simz,weight,dp,pmin,pmax,[100,1,1]);
    
end
toc
%matlabpool close

ClearGlobalVar;

%%
for i = 1:Nsim
    xrms(i) = sqrt(sum((pfit(i,1:2)-param(1,1:2)).^2));
end

%save('test07.mat','param','x','y','z','Nsim','pfit','Ntot','cm','xrms','c');
%%
% 
% h1 = figure;
% surf(z);
% 
% h2 = figure;
% filteredz = bpass(z,1,6);
% surf(filteredz)
% 
% h3 = figure;
% surf(simz)
% 
% %%
% h4 = figure;
% filtered = bpass(simz,1,6);
% surf(filtered)
% 
% 
% %% 
% placefigure(h1,[2,2,1,1]);
% placefigure(h2,[2,2,1,2]);
% placefigure(h3,[2,2,2,1]);
% placefigure(h4,[2,2,2,2]);


%%
