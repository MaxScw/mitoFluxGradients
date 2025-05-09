
%% Lighthouse problem, p29
clear all;
N=20; % Number of measurements
azims = (rand(N,1)-0.5).*pi; % generate random angles
plot(azims)


% Convert random angles into random shore measurements positions (x_k)
xks = @(al,bet) bet.*tan(azims)+al; 
hist(xks(0,1),30)

% likelihood of uniform azimuthal emission -> cauchy distribution in x_k
% (2.34)
Lkhd_x = @(xk,al,bet) bet./(pi.*(bet.^2+(xk-al).^2));

% Now try to find the best estimate of al, assuming bet=1
Pr_al = 1; % Prior, assume that alpha prob dist is a constant, could be anywhere.

% Use Bayes' theorem to caculate posterior. It's basically equal to the
% likelihood. Still formula for just one measurement. To do all
% measurements, multiply likelihoods of each observed x_k together because
% they are independent events.
Po_al = @(xk,al,bet) bet./(pi.*(bet.^2+(xk-al).^2)); 

% Take the log (from the formula in the book). This step accounts for all
% experimental measurements (or rand angles in this case) by 
L_al = @(xk,al,bet) -sum(log(bet.^2+(xk-al).^2));

% Now we're at the point where we can plug in our data points and
% numerically scan over the al parameter to find the maximum

% At this stage, generate virtual points by setting alpha to the "true"
% value, so we can generate the x's and use the analysis to identify it.
% xks(0,1) is all the x_k's generated from the initial azims.
% Change it around and you can see that the best estimate calculated below
% will be close to the initial alpha value you give.
i = 1;
for al = -10:.01:10
    L_al_num(i,1)= al;
    L_al_num(i,2)= L_al(xks(0,1),al,1);
    i=i+1;    
end
plot(L_al_num(:,1),L_al_num(:,2))
xlim([-4,4])
% Best estimate is given by max of log(posterior)
best = L_al_num(L_al_num(:,2)==max(L_al_num(:,2)),1)

% Standard deviation. Can calculate from 2nd derivative 
% dL_al = @(xk,al,bet) -2.*sum((xk-al)./(bet^2+(xk-al)^2));
% 
% hist(xks(1,1),40)


%% Lighthouse problem varying both beta and alpha
clear all;
N=100; % Number of measurements
azims = (rand(N,1)-0.5).*pi; % generate random angles
plot(azims)


% Convert random angles into random shore measurements positions (x_k)
xks = @(al_o,bet_o) bet_o.*tan(azims)+al_o; 
hist(xks(0,3),30)

% likelihood of uniform azimuthal emission -> cauchy distribution in x_k
% (2.34)
Lkhd_x = @(xk,al,bet) bet./(pi.*(bet.^2+(xk-al).^2));

% Now try to find the best estimate of al and bet. 
Pr_al = 1; % Prior, assume that alpha, beta prob dist is a constant, could be anywhere.

% Use Bayes' theorem to caculate posterior. It's basically equal to the
% likelihood. Still formula for just one measurement. To do all
% measurements, multiply likelihoods of each observed x_k together because
% they are independent events.
Po_al = @(xk,al,bet) bet./(pi.*(bet.^2+(xk-al).^2)); 

% Take the log (from the formula in the book). Contrary to formula in book,
% the beta term in the numerator now varies as we scan. It's no longer a
% constant, so we need to include it in the sum.
L_al = @(xk,al,bet) sum(log(bet)-log(bet.^2+(xk-al).^2));

% Now we're at the point where we can plug in our data points and
% numerically scan over the al parameter to find the maximum

% At this stage, generate virtual points by setting alpha to the "true"
% value, so we can generate the x's and use the analysis to identify it.
% xks(0,1) is all the x_k's generated from the initial azims.
% Change it around and you can see that the best estimate calculated below
% will be close to the initial alpha value you give.

i = 1;
j = 1;
al_o = 0;
bet_o = 4;
for al = -10:.1:10
    j=1;
    for bet = 0:.1:20
        L_al_num(j,i)= L_al(xks(al_o,bet_o),al,bet);
        j=j+1;
    end
    i=i+1;
end
[AL,BET] = meshgrid(-10:.1:10,0:.1:20);
surf(AL,BET,L_al_num,'EdgeColor','none');
view([0 0 1])

% Best estimate is given by max of log(posterior)
[a,b] = find(L_al_num==max(max(L_al_num)));
al_best = AL(a,b)
bet_best = BET(a,b)
% Standard deviation. Can calculate from 2nd derivative 
% dL_al = @(xk,al,bet) -2.*sum((xk-al)./(bet^2+(xk-al)^2));
% 
% hist(xks(1,1),40)

% Cool. Bayesian.... learned.


%% Poisson simulation of gaussian function with Amplitude (A) and 
% Background (B) and Bayesian
% analysis. Center distribution at 0, with a true amplitude of 1 and a
% background of 2.

%% Generate 4 different distributions of random numbers that follow the 
% Gaussian distribution, with A=1, B=2

% Just use randn
D1 = 1*normrnd(0,5,100,1)+2;
subplot(2,1,1) 
plot(D1);
subplot(2,1,2) 
hist(D1,15);




%% Generate Datum according to eq(3.1), with A=1, B=2
rands = rand(10000,1)-.5;
D1 = 1000.*(rands).*(1.*exp(-(1000.*(rands)).^2/(2.*5^2))+1);
subplot(2,1,1) 
plot(D1);
subplot(2,1,2) 
hist(D1,30);













%% Fail

% Measured distribution should be a gaussian distribution, but with poisson
% noise
Poiss = @(N,D) D^N*exp(-D)/factorial(N);

for i = 0:100
    P(i+1) = Poiss(i,30);
end
plot(P);

% Matlab function for generating random poisson numbers
for i = 0:10000
    D1(i+1) = poissrnd(1)+2;
end