% run the NLD

%% parameter description
% E0 = pars.E0;           %intial environment quality
% k0 = pars.k0;           %intial capital per unit labor
% neta = pars.neta;       %preference for environmental quality, neta>0
% b = pars.b;             %Environment'self rate of degradation; 0<b<1
% gamma = pars.gamma;     %Environmental maintainance efficiency
% tau = pars.tau;         %tax rate: extra payment to conserve environemtn per unit consumption, tau>=0
% delta = pars.delta;     %capital rate of depreciation: how capital corrodes over time, 0<delta<1
% phi = pars.phi;         %degree of habit formation of environmental quality, 0<=delta<1
% beta = pars.beta;       %degree of environmental quality caused by a unit of consumption,beta>0
% alpha = pars.alpha ;    %capital share of output, comes in the prodution fucntion y=Ak^alpha,, 0<alpha<1   
% N = pars.N ;            %number of steps of difference equations computed
% A = pars.A ;            %total factor of productivity, comes in the prodcution function y=Ak^alpha, A>0
% %k: capital per unit labor, k>=0;
% %E: environment quality
%%
% define pars first 
pars.E0 = 0.2;
pars.neta = 0.5;
pars.b = 0.2;
pars.gamma = 0.5;
pars.tau = 0.1;
pars.delta = 0.025;
pars.phi = 0.9;
pars.beta = 0.1;
pars.alpha = 0.1;
pars.N = 50;
pars.k0 = 0.01 ;
pars.A = 1;
%%
% call NLD 
[E,k,E_bar,k_bar,E_,k_] = NLD(pars);
%%
% plot E-k
figure ;
plot(k_,E_bar,'r-','LineWidth',2); hold on ;
plot(k_bar,E_,'k-','LineWidth',2);
xlabel('Capital per unit labor (Economic growth)','FontSize',12);
ylabel('Environment Quality','FontSize',12);
xlim([-0.1,1]);
ylim([-0.1,1]);

plot(k,E,'b-','LineWidth',1);
for i =1:1:length(k)
    plot(k(i),E(i),'b-','MarkerSize',10);
    if (i<=5)
        pause ;
    else 
        
       
        pause(0.1);
    end 
end 

hold off ;
figure; plot(E,'LineWidth',2)
xlabel('Time','FontSize',12)
ylabel('Environment Quality','FontSize',12)
figure; plot(k,'LineWidth',2)
xlabel('Time','FontSize',12)
ylabel('Economic Growth (captial/labor)','FontSize',12)

