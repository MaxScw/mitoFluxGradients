clear all; close all;
% Bayesian tutorial

s=[3;5]; 

%%

figure(1);
h=plot(s(1),s(2),'ro');  % Plot the flag as a red circle.
set(h,'markersize',6,'linewidth',3); % make the circle big.
axis([0,10,0,10]); % Set the scale of the plot so that we can see the origin.
hold on;
n=2*randn(2,100); % Create a 100-sample noise sequence with a standard deviation of 2.
x=zeros(2,100);
for (i=1:100)
    x(:,i)=s+n(:,i);  % Add the noise to the true state to create 100 observations of the true state.
    plot(x(1,i),x(2,i),'k.');
end;
hold off;

%%
sest=mean(x')';  % The ' indicates a transpose.  Because mean takes the
                 % average over the columns, I swap things around to get it to work.
hold on;
plot(sest(1),sest(2),'bs');  % Plot the average.
hold off;

%%
figure(2); % Switch to a new figure window.

sest=x(:,1); % The first estimate is just the first observation.  Draw it.
subplot(211); plot(1,sest(1)); hold on;
line([1,100],[s(1),s(1)]); % Draw a line at the location of the x component.
subplot(212); plot(1,sest(2)); hold on;
line([1,100],[s(2),s(2)]); % Draw a line at the location of the y component.

sold=sest;
for (n=2:100)
    sest = (n-1)/n * sold + 1/n * x(:,n);
    subplot(211);plot(n,sest(1),'k.');
    subplot(212); plot(n,sest(2),'k.');
    sold=sest;
end;
subplot(211); hold off; subplot(212);hold off;

%%
Sa=[2:1:4];
Sb=[4:1:6];
% To make the set of possible states more refined, uncomment the following.
%Sa=[2:0.1:4];
%Sb=[4:0.1:6];

%% Prior/posterior for x(:,1)

L=length(Sa);
Pr=ones(L,L); % Initialize the table to all ones
Po=ones(L,L);
Pr=Pr/sum(sum(Pr)); % Turn the table into a pmf by dividing by the sum.
Po=Po/sum(sum(Po)); % Each value is now 1/9.
%Pr=0*Pr;Pr(2,2)=1;

K=[4,0;0,4]; % covariance matrix.
m=0*Pr;
for (i=1:length(Pr))  % For each entry in my prior table.
    for (j=1:length(Pr))
        me=[Sa(i);Sb(j)];
        m(i,j) = 1/sqrt((2*pi)^2*det(K)) * exp(-(x(:,1)-me)'*inv(K)*(x(:,1)-me)/2); % Compute likelihood
        m(i,j) = m(i,j) * Pr(i,j); % Combine with prior
    end;
end;
Po=m/sum(sum(m));

%% Iterate

figure(3); % Switch to a new figure window.

[a,b]=find(Po==max(max(Po)));  % Pull out the indices at which Po achieves its max.
sest=[Sa(a);Sb(b)];  % The best estimate of the true state.

subplot(211); plot(1,sest(1)); hold on;
line([1,100],[s(1),s(1)]); % Draw a line at the location of the x component.
subplot(212); plot(1,sest(2)); hold on;
line([1,100],[s(2),s(2)]); % Draw a line at the location of the y component.

for (n=2:length(x));
    Pr=Po;
    AllPr{n}=Pr;
    m=0*Pr;
	for (i=1:length(Pr))  % For each entry in my prior table.
        for (j=1:length(Pr))
            me=[Sa(i);Sb(j)];
            m(i,j) = 1/sqrt((2*pi)^2*det(K)) * exp(-(x(:,n)-me)'*inv(K)*(x(:,n)-me)/2); %Compute likelihood
            m(i,j) = m(i,j) * Pr(i,j); % Combine with prior    
        end;
	end;
	Po=m/sum(sum(m));
    [a,b]=find(Po==max(max(Po)));  % Pull out the indices at which Po achieves its max.
    sest=[Sa(a);Sb(b)];  % The best estimate of the true state.
    subplot(211);plot(n,sest(1),'k.');
    subplot(212); plot(n,sest(2),'k.');
end;   
subplot(211); hold off;
subplot(212); hold off;





