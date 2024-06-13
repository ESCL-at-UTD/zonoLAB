% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Example:
%   The timed plotting of zonotopes, constrained zonotopes.
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

%% Zonotope
% n = 3;
% nG = 200;
% Z  = randomSet(1,'zono',n,nG);
% 
% figure;
% tStart = tic;
% [v,f] = plot(Z,'b',0.1);
% toc(tStart)

%% Constrained Zonotope
nG = 10;
nC = ceil(nG/2);
Zc  = randomSet(1,'conZono',n,nG,[],nC);

figure;
tStart = tic;
[v,f] = plot(Zc,'b',0.1);
toc(tStart)

% %% Plot sets
% 
% % Define (constrained) hypercubes bounding factors
% B = zono(eye(Z.nG),zeros(Z.nG,1));                          % For Z
% Bc = conZono(eye(Z.nG),zeros(Z.nG,1),A,b);                  % For Zc
% Bh = hybZono(eye(Z.nG),zeros(3,1),zeros(Z.nG,1),A,Ab,b);    % For Zh
% 
% % Define plot settings
% set(0,'defaultLineLineWidth', 2)
% set(0,'defaultAxesFontName' , 'Times')
% set(0,'defaultTextFontName' , 'Times')
% set(0,'defaultAxesFontSize' , 16)
% set(0,'defaultTextFontSize' , 16)
% set(0,'defaulttextinterpreter','latex')
% set(0,'defaultAxesGridLineStyle','-.')
% 
% 
% figure('Position',[-1700,100,600,800]);
% % Plot zonotope
% subplot(3,2,1)
% plot(B,'b',0.1)
% axis square
% axis([-1.1 1.1 -1.1 1.1 -1.1 1.1])
% xlabel('$\xi^c(1)$')
% ylabel('$\xi^c(2)$')
% zlabel('$\xi^c(3)$')
% title('(a)','position', [-2 1 1.75]);
% subplot(3,2,2)
% plot(Z,'b',0.1)
% axis square
% axis([-2.1 2.1 -2.1 2.1])
% xlabel('$z(1)$')
% ylabel('$z(2)$')
% title('(b)','position', [-3.5 1.5]);
% 
% 
% % Plot constrained zonotope
% subplot(3,2,3); hold on
% plot(B,'b',0)
% plot(Bc,'m',0.1)
% axis square
% axis([-1.1 1.1 -1.1 1.1 -1.1 1.1])
% xlabel('$\xi^c(1)$')
% ylabel('$\xi^c(2)$')
% zlabel('$\xi^c(3)$')
% title('(c)','position', [-2 1 1.75]);
% subplot(3,2,4)
% plot(Zc,'m',0.1)
% axis square
% axis([-2.1 2.1 -2.1 2.1])
% xlabel('$z(1)$')
% ylabel('$z(2)$')
% title('(d)','position', [-3.5 1.5]);
% 
% % Plot hybrid zonotope
% subplot(3,2,5); hold on
% plot(B,'b',0)
% plot(Bh,'g',0.1)
% axis square
% axis([-1.1 1.1 -1.1 1.1 -1.1 1.1])
% xlabel('$\xi^c(1)$')
% ylabel('$\xi^c(2)$')
% zlabel('$\xi^c(3)$')
% title('(e)','position', [-2 1 1.75]);
% subplot(3,2,6)
% plot(Zh,'g',0.1)
% axis square
% axis([-2.1 2.1 -2.1 2.1])
% xlabel('$z(1)$')
% ylabel('$z(2)$')
% title('(f)','position', [-3.5 1.5]);
% 
% % % Export plot
% % set(gcf,'color','w');
% % export_fig setDefinition.pdf