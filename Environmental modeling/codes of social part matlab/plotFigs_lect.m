%create functions & plot
c = 0.01:0.01:1;
E = 0.01:0.01:1;
%plot(c,log(c)+4.7,'k-','LineWidth',2);
% xlabel('consumption','FontSize',12);
% ylabel('Utility','FontSize',12);

% meshgrid: utility 2D
[X,Y] = meshgrid(c,E);
% Z = log(X)+log(Y);
% [C,h] = contour(X,Y,Z+8);
% set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2,'Color','Black','LineWidth',2);
% xlabel('consumption 1 (bitterballen)','FontSize',12);
% ylabel('consumption 2 (beer)','FontSize',12);

% Production function 2D
Z = 1.*X.^0.3.*Y.^0.7;      %Cobb-Douglas prduction function
[C,h] = contour(X,Y,Z,4);
set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2,'Color','Black','LineWidth',2);
xlabel('Input 1 (Labor)','FontSize',12);
ylabel('Input 2 (Capital)','FontSize',12);

% %Production function 1D
% plot(c,1.*c.^0.3,'k-','LineWidth',2);
% xlabel('Labour input','FontSize',12);
% ylabel('Production (y)','FontSize',12);

