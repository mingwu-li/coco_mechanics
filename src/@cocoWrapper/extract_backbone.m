function bd = extract_backbone(obj, omegaRange, varargin)
%EXTRACTBACKBONE This function extracts the backbone curve for a mode whose
%natural frequecy is closest to omega among all natural frequencies

warning('Backbone computation using coco assumes single symmetry');

% initial set up
dir_name = obj.Options.dir_name;
N = obj.system.N;
n = obj.system.n;
outdof = obj.outdof;

assert(obj.system.order == 2, 'fnl avaliable only for second-order systems')
for i=1:numel(obj.system.fnl)
    fnli = obj.system.fnl{i};
    if ~isempty(fnli)
        assert(size(fnli,2)==n, 'current implementation assumes f(x) instead of f(x,xd)');
    end
end

% find the natural frequency closest to omega
nEig = n;
[V, D, NOT_CONVERGED] = eigs(sparse(obj.system.K),sparse(obj.system.M),nEig,'smallestabs');
if NOT_CONVERGED
    error('The eigenvalue computation did not converge, please adjust the number of eigenvalues to be computed');
end
naturalFreq   = sqrt(diag(D));
[~, id_Freq] = min(abs(naturalFreq-omegaRange(1)));
priFreq = naturalFreq(id_Freq);
T = 2*pi/priFreq;
fprintf('Backbone curve around %d mode with freq %d will be calculated\n', id_Freq, priFreq);


% initial solution with zero damping and nonlinearity
if isempty(varargin)
    phi = V(:,id_Freq);
    phi = phi/sqrt(phi'*obj.system.M*phi);
    t0 = 0:T/200:T;
    x0 = phi*cos(priFreq*t0);
    v0 = -phi*priFreq*sin(priFreq*t0);
    z0 = [x0; v0];
else
    t0 = varargin{1};
    z0 = varargin{2};
end

% continuation in nonlinearity parameter
% setup coco
prob = coco_prob();
prob = cocoSet(obj, prob);
prob = coco_set(prob, 'ode', 'vectorized', true);
obj.fnlTensor2Multi();
odedata.fnl = obj.multiFnl;
odefun = @(z,p) obj.ode_aut(z,p,odedata);
funcs  = {odefun};

coll_args = {funcs{:}, t0', z0', {'epsilon', 'alpha'}, [0, 1]};   %#ok<CCAT>
prob = ode_isol2po(prob, '', coll_args{:});

% [data, uidx] = coco_get_func_data(prob, 'po.orb.coll', 'data', 'uidx');
% maps = data.coll_seg.maps;
% prob = coco_add_pars(prob, 'period', uidx(maps.x0_idx(1)), {'x0'});
% prob = coco_add_pars(prob, 'period', uidx(maps.T_idx), {'T'});

% prob = coco_add_event(prob, 'opt', 'alpha', '=', 1);
Trange = 2*pi./[omegaRange(2) omegaRange(1)];
cont_args = {1, { 'epsilon' 'po.period' 'alpha' }, {[],Trange}};

fprintf('\n Run=''%s'': Continue in nonlinear parameter.\n', ...
  dir_name);
runid = coco_get_id(dir_name, 'Backbone_nonlinear');
bd0 = coco(prob, runid, [], cont_args{:});


% continuation in damping parameter
lab = coco_bd_labs(bd0, 'opt');

prob = coco_prob();
prob = cocoSet(obj, prob);
prob = ode_po2po(prob, '', runid, lab);

[data, uidx] = coco_get_func_data(prob, 'po.orb.coll', 'data', 'uidx');
maps = data.coll_seg.maps;
% track amplitude of outdof
ampdata.dof  = outdof;
ampdata.zdim = N;
prob = coco_add_func(prob, 'amp', @amplitude, ampdata, 'regular', 'amp',...
    'uidx', uidx(maps.xbp_idx), 'remesh', @amplitude_remesh);


cont_args = {1, {'po.period' 'epsilon' 'amp'}, [0.8*T, 1.2*T]};
cont_args = {1, {'alpha' 'epsilon'},{[0 1]}};


fprintf('\n Run=''%s'': Continue in period.\n', ...
  dir_name);
runid = coco_get_id(dir_name, 'Backbone_omega');
bd = coco(prob, runid, [], cont_args{:});


% result visualization
hold on
omega = coco_bd_col(bd, 'po.period');
omega = 2*pi./omega;
amp   = coco_bd_col(bd, 'amp');
plot(omega, amp, 'rv', 'MarkerSize', 8,'LineWidth',2,'DisplayName','COCO');
legend('show');
% grid on; box on; 
% set(gca,'LineWidth',1.2);
% set(gca,'FontSize',14);
% xlabel('$\Omega$','interpreter','latex','FontSize',16);
% ylabel('$|z|_{\mathrm{out}}$','interpreter','latex','FontSize',16);
% title('Backbone curve by coco');

end

function [data, y] = amplitude(prob, data, u) %#ok<INUSL>

xbp = reshape(u, data.zdim, []);
y = xbp(data.dof,:);
y = max(abs(y));

end

function [prob, status, xtr] = amplitude_remesh(prob, data, chart, old_u, old_V) %#ok<INUSD>

[colldata, uidx] = coco_get_func_data(prob, 'po.orb.coll', 'data', 'uidx');
maps = colldata.coll_seg.maps;
xtr    = [];
prob   = coco_change_func(prob, data, 'uidx', uidx(maps.xbp_idx));
status = 'success';

end