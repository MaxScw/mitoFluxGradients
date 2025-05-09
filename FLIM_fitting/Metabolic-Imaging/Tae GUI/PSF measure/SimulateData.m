function simz = SimulateData(x,y,param,c)
% 
% x = XY{1};
% y = XY{2};
z = JointGaussianModel(x,y,param,c);

N = sum(z(:));
zp = randp(z(:),1,N);
zp = hist(zp,length(x)*length(y));
simz = reshape(zp,[length(y),length(x)]);
