function [E,k,E_bar,k_bar,E_,k_] = NLD(pars)
%this models the coupled economic growth- environment quality dynamics Chen and Li 2011 
% Based on difference equations 9 and 11.
% this is a John-Pecchenino model with habit formation of environmental
% quality and consumption tax

% intialize the parameters;
E0 = pars.E0;           %intial environment quality
k0 = pars.k0;           %intial capital per unit labor
neta = pars.neta;       %preference for environmental quality, neta>0
b = pars.b;             %Environment'self rate of degradation; 0<b<1
gamma = pars.gamma;     %Environmental maintainance efficiency
tau = pars.tau;         %tax rate: extra payment to conserve environemtn per unit consumption, tau>=0
delta = pars.delta;     %capital rate of depreciation: how capital corrodes over time, 0<delta<1
phi = pars.phi;         %degree of habit formation of environmental quality, 0<=delta<1
beta = pars.beta;       %degree of environmental quality caused by a unit of consumption,beta>0
alpha = pars.alpha ;    %capital share of output, comes in the prodution fucntion y=Ak^alpha,, 0<alpha<1   
N = pars.N ;            %number of steps of difference equations computed
A = pars.A ;            %total factor of productivity, comes in the prodcution function y=Ak^alpha, A>0
%k: capital per unit labor, k>=0;
%E: environment quality

% define the Environment quality E and capital per unit labour k matrices;
E = zeros(N,1);
k = zeros(N,1);

% intial conditions
E(1) = E0;
k(1) = k0;

%define the needed coefficients a0,a1,a2;
% a0 = (neta/(1+neta)) * (1- b - (beta - gamma*tau)*(1-delta)/(neta*gamma*(1+tau)) + phi/neta) ;
% a1 = phi*(beta-gamma*tau)*(1-delta)/((1+neta)*gamma*(1+tau));
% a2 = (neta^(1-alpha))/(1+neta)*A*(gamma*(1-alpha) - alpha*(beta-gamma*tau)/(1+tau))/(gamma^alpha);

%now the loop and put difference equations within
for i = 2:1:N
    f = A*(k(i-1)^alpha); %production ;
    wt = (1-alpha)*f;  %factor price of labor
    rt = alpha*A*(k(i-1)^(alpha-1)) - delta; %interest rate: factor price of capital ;
    if (k(i-1) == 0)
        rt = 0;
    end 
    Rt = 1+rt; %gross rate of return on savings
    if i>2 
    %E(i) = a0*E(i-1)+a1*E(i-2)+a2*(E(i-1)-phi*E(i-2));
        E(i) = (neta/(1+neta))*((1-b)*E(i-1) - (beta-gamma*tau)*(Rt/(neta*gamma*(1+tau)))*(E(i-1) - phi*E(i-2)) + gamma*(wt + phi*E(i-1)/(neta*gamma)));
    end 
    if i==2 
        %E(i) = a0*E(i-1)+a2*E(i-1);
        E(i) = (neta/(1+neta))*((1-b)*E(i-1) - (beta-gamma*tau)*(Rt/(neta*gamma*(1+tau)))*E(i-1) + gamma*(wt + phi*E(i-1)/(neta*gamma)));
    end    
    if (E(i)<=0) 
        E(i) = 0 ; 
    end 
    k(i) = (E(i)-phi*E(i-1))/(neta*gamma);
    if (k(i)<=0) 
        k(i) = 0 ; 
    end 
end 

%get 
k_ = 0.00:0.01:1;
E_ = 0.00:0.01:1;
E_bar = zeros(length(k_),1);
k_bar = zeros(length(E_),1);
%plot steady state 
for i = 1:1:length(k_)
    f = A*(k_(i)^alpha); %production ;
    wt = (1-alpha)*f;  %factor price of labor
    rt = alpha*A*(k_(i)^(alpha-1)) - delta; %interest rate: factor price of capital ;
    Rt = 1+rt;
    E_bar(i) = 1/(1 - (neta/(neta+1))*((1-b) - (beta-gamma*tau)*(Rt/(neta*gamma*(1+tau)))*(1-phi)+ gamma*phi/(neta*gamma))) ...
        * gamma*wt*(neta/(1+neta)) ;
    k_bar(i)
    k_bar(i) = E_(i)*(1-phi)/(neta*gamma);
end 
