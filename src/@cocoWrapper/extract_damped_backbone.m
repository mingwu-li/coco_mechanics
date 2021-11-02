function bds = extract_damped_backbone(obj,parRange,optdof,varargin)
%EXTRACTFRC This function extracts forced response curve with po toolbox of
% COCO
% BD = EXTRACT_FRC(OBJ, PARRANGE, OPTDOF, VARARGIN)
%
% PARRANGE: computational domain of parameters {[freq1,freq2],[eps1,eps2]}
% OPTDOF:   amplitude for a set of degrees-of-freedom
% VARARGIN: initSol.(t0,x0,p0) initial solution guess

%% initial setup
N = obj.system.N;
n = obj.system.n;
omega0 = parRange{1}(1); % start from left end point if initial guess is not specified
T0     = 2*pi/omega0*obj.periodsRatio;
tf     = obj.nCycles*T0;
epsilon = obj.system.fext.epsilon;
optnorm = obj.optnorm;
if strcmp(optnorm,'linf')
    assert(numel(optdof)==1,'L_infty supports only for a single dof');
end
assert(numel(obj.system.Omega)<=1, 'coco run assumes single freq component');
assert(obj.system.order == 2, 'fnl avaliable only for second-order systems')
for i=1:numel(obj.system.fnl)
    fnli = obj.system.fnl{i};
    if ~isempty(fnli)
        assert(size(fnli,2)==n, 'current implementation assumes f(x) instead of f(x,xd)');
    end
end

obj.fnlTensor2Multi();
odedata.fnl = obj.multiFnl;
switch obj.initialGuess
    case 'forward'
        %% initial solution by forward simulation
        % ode45 is used here. Integration option may be added in future
        x0_init = zeros(N,1);
        odefun = @(t,x) obj.ode_het(t,x,[omega0;epsilon;0],odedata);
        [~, x0_po] = ode15s(odefun , [0 tf], x0_init);                % transient
        options = odeset('RelTol', 1e-9, 'AbsTol',1e-9);
        [t0, x0] = ode45(odefun, [0 T0], x0_po(end,:)', options);    % steady state
        p0 = [omega0,epsilon,0];
    case 'linear'
        %% initial solution by solving linear equation M\ddot{x}+C\dot{x}+Kx = F cos(O*t)
        mass = obj.system.M;
        damp = obj.system.C;
        stif = obj.system.K;
        fext = 2*obj.system.fext.coeffs;
        fext = epsilon*fext;
        kapa = obj.system.fext.kappas(1);
        if kapa<0
            kapa = -kapa;
        end
        qcom = (-(kapa*omega0)^2*mass+1i*kapa*omega0*damp+stif)\fext(:,1);
        t0   = linspace(0,T0,100);
        solx = qcom*exp(1i * kapa * omega0 * t0);
        solv = 1i * kapa * omega0 *solx;
        x0   = real([solx; solv])';
        p0   = [omega0,epsilon,0];
    case 'given'
        assert(numel(varargin)>0,'initial guess is not provided');
        initSol = varargin{1};
        t0 = initSol.t0;
        x0 = initSol.x0;
        p0 = initSol.p0;
end
if strcmp(optnorm,'linf')
    bds = damped_backbone_linfnorm(obj,odedata,t0,x0,p0,optdof,parRange);
else
    bds = damped_backbone_l2norm(obj,odedata,t0,x0,p0,optdof,parRange);
end
% result visualization
% extract results
omega = [];
epsf  = [];
obj   = [];
stab  = logical([]);
for i=1:numel(bds)
    bd = bds{i};
    omegai = coco_bd_col(bd, 'omega');
    epsfi = coco_bd_col(bd, 'eps');
    stabi = coco_bd_col(bd, 'eigs')';
    stabi = abs(stabi);
    stabi = all(stabi<1, 2);
    obji  = coco_bd_col(bd, 'obj');
    if strcmp(optnorm,'linf')
        obji = abs(obji);
    else
        obji = sqrt(obji); 
    end
    if i==1
        figure;
        plot3(omegai(stabi), epsfi(stabi), obji(stabi), 'b*', 'MarkerSize', 6,'LineWidth',2,'DisplayName','FRC-COCO-stable');
        plot3(omegai(~stabi), epsfi(~stabi), obji(~stabi), 'r*', 'MarkerSize', 6,'LineWidth',2,'DisplayName','FRC-COCO-unstable');
    else
        omega = [omega, omegai];
        epsf  = [epsf, epsfi];
        stab  = [stab; stabi];
        obj   = [obj, obji];
    end
end
% plot results
figure(gcf); hold on
plot3(omega(stab), epsf(stab), obj(stab), 'g*', 'MarkerSize', 6,'LineWidth',2,'DisplayName','BC-COCO-stable');
plot3(omega(~stab), epsf(~stab), obj(~stab), 'c*', 'MarkerSize', 6,'LineWidth',2,'DisplayName','BC-COCO-unstable');
xlabel('$\Omega$','interpreter','latex','FontSize',14);
ylabel('$\epsilon$','interpreter','latex','FontSize',14);
zlabel('$||x_\mathrm{opt}||$','interpreter','latex','FontSize',14);
grid on; box on; 
legend('show');
end