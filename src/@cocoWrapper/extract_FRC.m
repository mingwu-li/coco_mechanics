function bds = extract_FRC(obj,parName,parRange, varargin)
%EXTRACTFRC This function extracts forced response curve with po toolbox of
% COCO
% BD = EXTRACT_FRC(OBJ, PARNAME, PARRANGE, VARARGIN)
%
% PARNAME:  'freq/amp' continuation with varied excitation frequency or
% amplitude
% PARRANGE: computational domain of parameters
% VARARGIN: initSol.(t0,x0,p0) initial solution guess

%% initial setup
dir_name = obj.Options.dir_name;
N = obj.system.N;
n = obj.system.n;
if strcmp(parName,'freq')
omega0 = parRange(1); % start from left end point if initial guess is not specified
end
T0     = 2*pi/omega0*obj.periodsRatio;
tf     = obj.nCycles*T0;
outdof = obj.outdof;
epsilon = obj.system.fext.epsilon;

assert(numel(obj.system.Omega)<=1, 'coco run assumes single freq component');
assert(obj.system.order == 2, 'fnl avaliable only for second-order systems')
for i=1:numel(obj.system.fnl)
    fnli = obj.system.fnl{i};
    if ~isempty(fnli)
        assert(size(fnli,2)==n, 'current implementation assumes f(x) instead of f(x,xd)');
    end
end

switch obj.initialGuess
    case 'forward'
        %% initial solution by forward simulation
        % ode45 is used here. Integration option may be added in future
        x0_init = zeros(N,1);
        obj.fnlTensor2Multi();
        odedata.fnl = obj.multiFnl;
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
        if numel(p0)==2
            p0 = [p0;0];
        end
end

%% continuation excitation frequency
% setup coco
prob = coco_prob();
if strcmp(obj.atlasAlg,'kd')
    coco_func_data.pointers('set', []); % only necessary for atlas_kd
    prob = coco_set(prob, 'all', 'CleanData', true);
end
prob = cocoSet(obj, prob);
prob = coco_set(prob, 'ode', 'autonomous', false);
prob = coco_set(prob, 'ode', 'vectorized', true);
obj.fnlTensor2Multi();
odedata.fnl = obj.multiFnl;
odefun = @(t,x,p) obj.ode_het(t,x,p,odedata);
odefun_dx = @(t,x,p) obj.ode_het_dx(t,x,p,odedata);
odefun_dp = @(t,x,p) obj.ode_het_dp(t,x,p,odedata);
odefun_dt = @(t,x,p) obj.ode_het_dt(t,x,p,odedata);
funcs = {odefun,odefun_dx,odefun_dp,odefun_dt};

coll_args = {funcs{:}, t0, x0, {'omega','eps','th'}, p0};   %#ok<CCAT>
prob = ode_isol2po(prob, '', coll_args{:});
[fdata, uidx] = coco_get_func_data(prob, 'po.orb.coll', 'data', 'uidx');
maps = fdata.coll_seg.maps;
omData = struct();
omData.periodsRatio = obj.periodsRatio;
prob = coco_add_func(prob, 'OmegaT', @OmegaT, @OmegaT_du, omData, 'zero',...
    'uidx', [uidx(maps.T_idx), uidx(maps.p_idx(1))]);

if strcmp(obj.atlasAlg,'kd')
    % integral functional
%     intdata = int_init_data(fdata, 'po.orb');
%     maps    = intdata.coll_seg.maps;
%     intuidx = uidx([maps.xbp_idx; maps.T_idx]);
%     prob = coco_add_func(prob, 'po.orb.int', @int, @int_du, intdata, ...
%       'inactive', 'po.L2norm', 'uidx', intuidx, 'remesh', @int_remesh);
    % initial state
    x0_idx = uidx(maps.xbp_idx);
    x0_idx = x0_idx(1:fdata.xdim);
    x0_names = cell(1,fdata.xdim);
    for k = 1:fdata.xdim
       x0_names{k} = strcat('x0_',num2str(k)); 
    end
    % monitor function for initial state
%     prob = coco_add_pars(prob,'initial_states',x0_idx, x0_names, 'active');
    % monitor function for normalized initial state: x0/numel(x0)
    prob = coco_add_func(prob, 'initial_states',@init,[],'active',x0_names,'uidx',x0_idx);
end

% track amplitude of outdof
ampdata.dof  = outdof;
ampdata.zdim = N;
numoutdof = numel(outdof);
ampNames = cell(1, numel(numoutdof));
for k = 1:numel(outdof)
   ampNames{k} = strcat('amp',num2str(outdof(k))); 
end
prob = coco_add_func(prob, 'amp', @amplitude, ampdata, 'regular', ampNames,...
    'uidx', uidx(maps.xbp_idx), 'remesh', @amplitude_remesh);
switch parName
    case 'freq'
        cont_args = {1, [{'omega'} {'po.period'} {'eps'} ampNames(:)'], [parRange(1) parRange(end)]};
        runid = coco_get_id(dir_name, 'freqFRC');
    case 'amp'
        cont_args = {1, [{'eps'} {'po.period'} {'omega'} ampNames(:)'], [parRange(1) parRange(end)]};
        runid = coco_get_id(dir_name, 'ampFRC');
    otherwise
        error('parName should be freq or amp');
end
    
fprintf('\n Run=''%s'': Continue primary family of periodic orbits.\n', ...
  runid);
bd0    = coco(prob, runid, [], cont_args{:});

% continuation along secondary branch if branch points are detected along
% the continuation above.

if obj.branchSwitch
    BPlab = coco_bd_labs(bd0, 'BP');
    numBP = numel(BPlab);
    bds = cell(numBP+1,1);
    runids = cell(numBP+1,1);
    bds{1} = bd0;
    runids{1} = runid;
    for k=1:numBP
        runidk = coco_get_id(runid, num2str(k));
        prob = coco_prob();
        prob = cocoSet(obj, prob);
        prob = coco_set(prob, 'ode', 'autonomous', false);
        prob = coco_set(prob, 'ode', 'vectorized', true);
        prob = ode_BP2po(prob, '', runid, BPlab(k));
        [fdata, uidx] = coco_get_func_data(prob, 'po.orb.coll', 'fdata', 'uidx');
        maps = fdata.coll_seg.maps;
        prob = coco_add_func(prob, 'OmegaT', @OmegaT, @OmegaT_du, omData, 'zero',...
            'uidx', [uidx(maps.T_idx), uidx(maps.p_idx(1))]);
        intdata = int_init_data(fdata, 'po.orb');
        maps    = intdata.coll_seg.maps;
        intuidx = uidx([maps.xbp_idx; maps.T_idx]);
        prob = coco_add_func(prob, 'po.orb.int', @int, @int_du, intdata, ...
          'inactive', 'po.L2norm', 'uidx', intuidx, 'remesh', @int_remesh);
        prob = coco_add_func(prob, 'amp', @amplitude, ampdata, 'regular', ampNames,...
            'uidx', uidx(maps.xbp_idx), 'remesh', @amplitude_remesh);
        fprintf('\n Run=''%s'': Continue equilibria along secondary branch.\n', ...
                  runidk);
        bdk = coco(prob, runidk, [], cont_args{:});
        bds{k+1} = bdk;
        runids{k+1} = runidk;
    end
else
    numBP = 0;
    bds = {bd0};
    runids = {runid};
end


% result visualization
% extract results
omega = [];
epsf  = [];
stab = logical([]);
amp = cell(numoutdof,1);
for i=1:numBP+1
    bd = bds{i};
    omegai = coco_bd_col(bd, 'omega');
    epsfi = coco_bd_col(bd, 'eps');
    stabi = coco_bd_col(bd, 'eigs')';
    stabi = abs(stabi);
    stabi = all(stabi<1, 2);
    omega = [omega, omegai];
    epsf  = [epsf, epsfi];
    stab  = [stab; stabi];
    for j=1:numoutdof
        ampj   = coco_bd_col(bd, ampNames{j});
        amp{j} = [amp{j}, ampj];
    end
end

% plot results
figure(gcf); hold on
if numoutdof>1
    for k=1:numoutdof
        subplot(numoutdof,1,k); hold on
        if strcmp(parName,'freq')
            plot(omega(stab), amp{k}(stab), 'g*', 'MarkerSize', 6,'LineWidth',2,'DisplayName','COCO-stable');
            plot(omega(~stab), amp{k}(~stab), 'c*', 'MarkerSize', 6,'LineWidth',2,'DisplayName','COCO-unstable');
        else
            plot(epsf(stab), amp{k}(stab), 'g*', 'MarkerSize', 6,'LineWidth',2,'DisplayName','COCO-stable');
            plot(epsf(~stab), amp{k}(~stab), 'c*', 'MarkerSize', 6,'LineWidth',2,'DisplayName','COCO-unstable');
        end
    end
else
    if strcmp(parName,'freq')
        plot(omega(stab), amp{1}(stab), 'g*', 'MarkerSize', 6,'LineWidth',2,'DisplayName','COCO-stable');
        plot(omega(~stab), amp{1}(~stab), 'c*', 'MarkerSize', 6,'LineWidth',2,'DisplayName','COCO-unstable');
    else
        plot(epsf(stab), amp{1}(stab), 'g*', 'MarkerSize', 6,'LineWidth',2,'DisplayName','COCO-stable');
        plot(epsf(~stab), amp{1}(~stab), 'c*', 'MarkerSize', 6,'LineWidth',2,'DisplayName','COCO-unstable');
    end        
end
legend('show');

% figure;
% thm = struct( 'special', {{'FP', 'PD', 'TR'}});
% thm.FP = {'LineStyle', 'none', 'LineWidth', 2, ...
%   'Color', 'cyan', 'Marker', 'v', 'MarkerSize', 10, 'MarkerEdgeColor', ...
%   'cyan', 'MarkerFaceColor', 'white'};
% thm.PD = {'LineStyle', 'none', 'LineWidth', 2, ...
%   'Color', 'black', 'Marker', 'o', 'MarkerSize', 8, 'MarkerEdgeColor', ...
%   'black', 'MarkerFaceColor', 'white'};
% thm.TR = {'LineStyle', 'none', 'LineWidth', 2, ...
%   'Color', 'red', 'Marker', 'd', 'MarkerSize', 8, 'MarkerEdgeColor', ...
%   'red', 'MarkerFaceColor', 'white'};
% for i=1:numBP+1
%     if isempty(varargin)
%         coco_plot_bd(thm, runids{i}, 'omega', ampNames{1})
%     else
%         coco_plot_bd(thm, runids{i}, 'eps', ampNames{1})
%     end
% end
% grid on; box on; 
% set(gca,'LineWidth',1.2);
% set(gca,'FontSize',14);
% if isempty(varargin)
% xlabel('$\Omega$','interpreter','latex','FontSize',16);
% else
% xlabel('$\epsilon$','interpreter','latex','FontSize',16);
% end    
% ylabel(strcat('$||z_{',num2str(outdof(1)),'}||_{\infty}$'),'interpreter','latex','FontSize',16);
% title('FRC by coco(solid/dashed - stable/unstable)');
end