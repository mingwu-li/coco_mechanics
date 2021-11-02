function bds = damped_backbone_linfnorm(obj,odedata,t0,x0,p0,optdof,parRange)

outdof   = obj.outdof;
dir_name = obj.Options.dir_name;
%% continuation in theta
% setup coco
prob = coco_prob();
prob = cocoSet(obj, prob);
prob = coco_set(prob, 'ode', 'autonomous', false);
prob = coco_set(prob, 'ode', 'vectorized', true);

%% zero problem
odefun    = @(t,x,p) obj.ode_het(t,x,p,odedata);
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

% objective functional
xoptidx = uidx(maps.x0_idx);
xoptidx = xoptidx(optdof);
prob    = coco_add_pars(prob, 'obj', xoptidx, 'obj');

% track amplitude of outdof
ampdata.dof  = outdof;
ampdata.zdim = obj.system.N;
numoutdof = numel(outdof);
ampNames = cell(1, numel(numoutdof));
for k = 1:numel(outdof)
   ampNames{k} = strcat('amp',num2str(outdof(k))); 
end
prob = coco_add_func(prob, 'amp', @amplitude, ampdata, 'regular', ampNames,...
    'uidx', uidx(maps.xbp_idx), 'remesh', @amplitude_remesh);


%% adjoints
prob = adjt_isol2po(prob, '');

[data, axidx] = coco_get_adjt_data(prob, 'po.orb.coll', 'data', 'axidx');
opt = data.coll_opt;
prob = coco_add_adjt(prob, 'OmegaT', 'aidx', ...
  axidx([opt.T_idx; opt.p_idx(1)]));

prob   = coco_add_adjt(prob, 'obj', 'd.obj', 'aidx', axidx(opt.x0_idx(optdof)));

cont_args = {1, [{'obj'} {'th'} {'po.period'} {'d.obj'} {'d.po.tinit'}...
    {'d.eps'} {'d.omega'} {'omega'} {'eps'} ampNames(:)']};
runid = coco_get_id(dir_name, 'phase');    
fprintf('\n Run=''%s'': Continue primary family of periodic orbits.\n', ...
  runid);
bd0    = coco(prob, runid, [], cont_args{:},{[], [-pi,pi]});

%% branch switch and continuation with d.obj=1
BPlab = coco_bd_labs(bd0, 'BP');
% find the one with maximum manginitude
objv  = coco_bd_col(bd0, 'obj');
BPidx = coco_bd_idxs(bd0, 'BP');
[~,idx] = max(abs(objv(BPidx)));
BPlab = BPlab(idx);

prob = coco_prob();
prob = cocoSet(obj, prob);
prob = coco_set(prob, 'ode', 'autonomous', false);
prob = coco_set(prob, 'ode', 'vectorized', true);
% zero problem
prob = ode_BP2po(prob, '', runid, BPlab);
[fdata, uidx] = coco_get_func_data(prob, 'po.orb.coll', 'data', 'uidx');
maps = fdata.coll_seg.maps;
prob = coco_add_func(prob, 'OmegaT', @OmegaT, @OmegaT_du, omData, 'zero',...
'uidx', [uidx(maps.T_idx), uidx(maps.p_idx(1))]);
xoptidx = uidx(maps.x0_idx);
xoptidx = xoptidx(optdof);
prob    = coco_add_pars(prob, 'obj', xoptidx, 'obj');
prob    = coco_add_func(prob, 'amp', @amplitude, ampdata, 'regular', ampNames,...
'uidx', uidx(maps.xbp_idx), 'remesh', @amplitude_remesh);  
% branch switch data
chart = coco_read_solution(runid, BPlab, 'chart');
cdata = coco_get_chart_data(chart, 'lsol');
% adjoint
prob = adjt_BP2po(prob, '', runid, BPlab);
[chart, aidx] = coco_read_adjoint('OmegaT', runid, BPlab, 'chart', 'lidx');
[data, axidx] = coco_get_adjt_data(prob, 'po.orb.coll', 'data', 'axidx');
opt = data.coll_opt;
prob = coco_add_adjt(prob, 'OmegaT', 'aidx', ...
axidx([opt.T_idx; opt.p_idx(1)]),'l0', chart.x, 'tl0', cdata.v(aidx));
[chart, aidx] = coco_read_adjoint('obj', runid, BPlab, 'chart', 'lidx');
prob = coco_add_adjt(prob, 'obj', 'd.obj', 'aidx', axidx(opt.x0_idx(optdof)), ...
    'l0', chart.x, 'tl0', cdata.v(aidx));
% continuation run
contBP_args = {1, [{'d.obj'} {'th'} {'po.period'} {'obj'} {'d.po.tinit'}...
    {'d.eps'} {'d.omega'} {'omega'} {'eps'} ampNames(:)']};
runidBP = coco_get_id(dir_name, 'phaseBP');    
bd1 = coco(prob, runidBP, [], contBP_args{:},[0,1]);

%% continuation in om with d.obj=1 
EPlab = coco_bd_labs(bd1,'EP');
EPlab = max(EPlab);
prob = coco_prob();
prob = cocoSet(obj, prob);
prob = coco_set(prob, 'ode', 'autonomous', false);
prob = coco_set(prob, 'ode', 'vectorized', true);
% zero problem
prob = ode_po2po(prob, '', runidBP, EPlab);
[fdata, uidx] = coco_get_func_data(prob, 'po.orb.coll', 'data', 'uidx');
maps = fdata.coll_seg.maps;
prob = coco_add_func(prob, 'OmegaT', @OmegaT, @OmegaT_du, omData, 'zero',...
'uidx', [uidx(maps.T_idx), uidx(maps.p_idx(1))]);
xoptidx = uidx(maps.x0_idx);
xoptidx = xoptidx(optdof);
prob    = coco_add_pars(prob, 'obj', xoptidx, 'obj');
prob = coco_add_func(prob, 'amp', @amplitude, ampdata, 'regular', ampNames,...
'uidx', uidx(maps.xbp_idx), 'remesh', @amplitude_remesh);  
% adjoint
prob = adjt_po2po(prob, '', runidBP, EPlab);
chart = coco_read_adjoint('OmegaT', runidBP, EPlab, 'chart');
[data, axidx] = coco_get_adjt_data(prob, 'po.orb.coll', 'data', 'axidx');
opt = data.coll_opt;
prob = coco_add_adjt(prob, 'OmegaT', 'aidx', ...
axidx([opt.T_idx; opt.p_idx(1)]),'l0', chart.x);
chart = coco_read_adjoint('obj', runidBP, EPlab, 'chart');
prob = coco_add_adjt(prob, 'obj', 'd.obj', 'aidx', axidx(opt.x0_idx(optdof)), ...
    'l0', chart.x);
prob = coco_add_event(prob, 'OPT', 'd.omega', '=', 0);
% continuation run
runidom = coco_get_id(dir_name, 'omega'); 
contom_args = {1, [{'d.omega'} {'omega'} {'th'} {'po.period'} {'obj'} {'d.po.tinit'}...
    {'d.eps'} {'d.obj'} {'eps'} ampNames(:)']};
bd2 = coco(prob, runidom, [], contom_args{:},{[], parRange{1}});
    
%% continuation in (epf,om) with d.obj=1 and d.om=0    
OPTlab = coco_bd_labs(bd2, 'OPT');
numOPT = numel(OPTlab);
bds    = cell(numOPT+1,1);
bds{1} = bd2;
contBC_args = {1, [{'d.eps'} {'eps'} {'omega'} {'th'} {'po.period'} {'obj'} {'d.po.tinit'}...
      {'d.omega'} {'d.obj'} ampNames(:)']};
runidBC = coco_get_id(dir_name,'epf');
for k=1:numOPT
    %% branch switch
    runidepfk = [runidBC, num2str(k)];
    prob = coco_prob();
    prob = cocoSet(obj, prob);
    prob = coco_set(prob, 'ode', 'autonomous', false);
    prob = coco_set(prob, 'ode', 'vectorized', true);
    % zero problem
    prob = ode_BP2po(prob, '', runidom, OPTlab(k));
    [fdata, uidx] = coco_get_func_data(prob, 'po.orb.coll', 'data', 'uidx');
    maps = fdata.coll_seg.maps;
    prob = coco_add_func(prob, 'OmegaT', @OmegaT, @OmegaT_du, omData, 'zero',...
        'uidx', [uidx(maps.T_idx), uidx(maps.p_idx(1))]);
    xoptidx = uidx(maps.x0_idx);
    xoptidx = xoptidx(optdof);
    prob    = coco_add_pars(prob, 'obj', xoptidx, 'obj');
    prob = coco_add_func(prob, 'amp', @amplitude, ampdata, 'regular', ampNames,...
        'uidx', uidx(maps.xbp_idx), 'remesh', @amplitude_remesh);  
    % adjoint
    prob = adjt_BP2po(prob, '', runidom, OPTlab(k));
    chart = coco_read_adjoint('OmegaT', runidom, OPTlab(k), 'chart');
    [data, axidx] = coco_get_adjt_data(prob, 'po.orb.coll', 'data', 'axidx');
    opt  = data.coll_opt;
    prob = coco_add_adjt(prob, 'OmegaT', 'aidx', ...
      axidx([opt.T_idx; opt.p_idx(1)]),'l0', chart.x);
    chart = coco_read_adjoint('obj', runidom, OPTlab(k), 'chart');
    prob  = coco_add_adjt(prob, 'obj', 'd.obj', 'aidx', axidx(opt.x0_idx(optdof)), ...
        'l0', chart.x); 
    prob = coco_add_event(prob, 'BC', 'BP', 'eps', '<', parRange{2}(1));
    % continuation run
    bdk = coco(prob, runidepfk, [], contBC_args{:},{[], parRange{2}, parRange{1}});
    bds{k+1} = bdk;
end

end