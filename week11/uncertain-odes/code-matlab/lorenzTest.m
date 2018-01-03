clear all;

dim=3;

% number of particles
N=100;

% number of frames
nFrames = 8000;

% time step
dt=0.01;

% time line
t=0:dt:nFrames*dt;

% inital position
X0 = [0;1;0];

% inital covariance around X0
sigma=0.1;
Sigma=(sigma^2)*eye(dim);

% initail points:
X=NormalSample([0 1 0],Sigma,N);  % without sampling
XS=X;                           % with sampling

% ODE solver options
% options = odeset('RelTol',1e-3,'AbsTol',[1e-3 1e-3 1e-4]);

% R(1,:) the bigest radius without samplig
% R(2,:) -- with sampling
Rad = zeros(2,nFrames);

Rad(:,1)=max(eig(Sigma).^0.5);

tic % seems too slow
for tau=1:nFrames
    
    for j=1:N
        [~,sol]=ode45(@lorenz,[0,dt],X(:,j));%,options);
        X(:,j)=sol(end,:);
        
        [~,sols]=ode45(@lorenz,[0,dt],XS(:,j));%,options);
        XS(:,j)=sols(end,:);
    end
    
    Sigma=cov(X');

    muS = mean(XS,2);
    SigmaS = cov(XS');
            
    la=eig(Sigma);
    laS=eig(SigmaS); % redundant
    
    Rad(:,tau+1)=[max(la.^0.5);max(laS.^0.5)];
    
    XS=NormalSample(muS,SigmaS,N);

end
toc

plot(t,Rad(1,:),t, Rad(2,:));
h = legend('no sampl.','sampl',2);
set(h,'Interpreter','none')
