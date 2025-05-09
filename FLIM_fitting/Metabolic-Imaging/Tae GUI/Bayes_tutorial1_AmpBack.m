

%% Lighthouse problem

% al = ?, bet = 1

plot(Lkhd_x(0:.1:10,5,1)) % likelihood funtion

% Find al, use prior that lighthouse could be anywhere, Pr_al = A
% Bayes' theorem:
Pr_al = 1;
% Po_al = Lkhd_x.*Pr_al;
% This is 
% L_al = @(xk,al,bet) -sum(log(bet.^2+(xk-al).^2));

%% Lighthouse problem

N=100; % Number of measurements
azims = (rand(N,1)-0.5).*pi; % generate random angles
plot(azims)


% Convert random angles into random shore measurements positions (x_k)
xks = @(al,bet) bet.*tan(azims)+al; 
hist(xks(1,1))

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
% Now we're at the point where we can plug in our 
L_al(xks(1,1),1,1)
dL_al = @(xk,al,bet) -2.*sum((xk-al)./(bet^2+(xk-al)^2));

hist(xks(1,1),40)
















%% Poisson simulation of Amplitude (A) and Background noise (B) and Bayesian
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