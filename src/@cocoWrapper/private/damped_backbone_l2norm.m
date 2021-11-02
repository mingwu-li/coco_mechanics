function bds = damped_backbone_l2norm(obj,odedata,t0,x0,p0,optdof,parRange)

outdof   = obj.outdof;
dir_name = obj.Options.dir_name;
%% continuation excitation frequency
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
data  = obj_init_data(fdata,optdof);
prob = coco_add_func(prob, 'obj', @objhan, @objhan_du, data, ...
  'inactive', 'obj', 'uidx', uidx, 'remesh', @obj_remesh);

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

data = adj_obj_init_data(data,optdof);
prob = coco_add_adjt(prob, 'obj', @adj_objhan, @adj_objhan_du, data, ...
  'd.obj', 'aidx', axidx([opt.xcn_idx; opt.T0_idx; opt.T_idx; opt.p_idx]), ...
  'remesh', @adj_obj_remesh);

cont_args = {1, [{'obj'} {'omega'} {'po.period'} {'d.obj'} {'d.po.tinit'}...
    {'d.eps'} {'d.th'} {'eps'} ampNames(:)']};
runid = coco_get_id(dir_name, 'freqFRC');    
fprintf('\n Run=''%s'': Continue primary family of periodic orbits.\n', ...
  runid);
bd0    = coco(prob, runid, [], cont_args{:},{[], [0.99*parRange{1}(1),parRange{1}(2)]});

%% branch switch and continuation with d.obj=1
BPlab = coco_bd_labs(bd0, 'BP');
numBP = numel(BPlab);
bds   = cell(numBP+1,1);
bds{1} = bd0;
contBP_args = {1, [{'d.obj'} {'omega'} {'po.period'} {'obj'} {'d.po.tinit'}...
    {'d.eps'} {'d.th'} {'eps'} ampNames(:)']};
contBC_args = {1, [{'d.eps'} {'eps'} {'omega'} {'po.period'} {'obj'} {'d.po.tinit'}...
      {'d.th'} {'d.obj'} ampNames(:)']};
runidBP = coco_get_id(runid,'BP');
runidBC = coco_get_id(runid,'BC');
for k=1:numBP
    %% branch switch
    runidBPk = [runidBP, num2str(k)];
    prob = coco_prob();
    prob = cocoSet(obj, prob);
    prob = coco_set(prob, 'ode', 'autonomous', false);
    prob = coco_set(prob, 'ode', 'vectorized', true);
    % zero problem
    prob = ode_BP2po(prob, '', runid, BPlab(k));
    [fdata, uidx] = coco_get_func_data(prob, 'po.orb.coll', 'data', 'uidx');
    maps = fdata.coll_seg.maps;
    prob = coco_add_func(prob, 'OmegaT', @OmegaT, @OmegaT_du, omData, 'zero',...
        'uidx', [uidx(maps.T_idx), uidx(maps.p_idx(1))]);
    data  = obj_init_data(fdata,optdof);
    prob = coco_add_func(prob, 'obj', @objhan, @objhan_du, data, ...
      'inactive', 'obj', 'uidx', uidx, 'remesh', @obj_remesh);
    prob = coco_add_func(prob, 'amp', @amplitude, ampdata, 'regular', ampNames,...
        'uidx', uidx(maps.xbp_idx), 'remesh', @amplitude_remesh);  
    % branch switch data
    chart = coco_read_solution(runid, BPlab(k), 'chart');
    cdata = coco_get_chart_data(chart, 'lsol');
    % adjoint
    prob = adjt_BP2po(prob, '', runid, BPlab(k));
    [chart, aidx] = coco_read_adjoint('OmegaT', runid, BPlab(k), 'chart', 'lidx');
    [data, axidx] = coco_get_adjt_data(prob, 'po.orb.coll', 'data', 'axidx');
    opt = data.coll_opt;
    prob = coco_add_adjt(prob, 'OmegaT', 'aidx', ...
      axidx([opt.T_idx; opt.p_idx(1)]),'l0', chart.x, 'tl0', cdata.v(aidx));
    [chart, aidx] = coco_read_adjoint('obj', runid, BPlab(k), 'chart', 'lidx');
    data = adj_obj_init_data(data,optdof);
    prob = coco_add_adjt(prob, 'obj', @adj_objhan, @adj_objhan_du, data, ...
      'd.obj', 'aidx', axidx([opt.xcn_idx; opt.T0_idx; opt.T_idx; opt.p_idx]), ...
      'l0', chart.x, 'tl0', cdata.v(aidx), 'remesh', @adj_obj_remesh);   
    % continuation run
    bdk = coco(prob, runidBPk, [], contBP_args{:},[0,1]);
    %% continuation in (epf,om) with d.obj=1 
    EPlab = coco_bd_labs(bdk,'EP');
    EPlab = max(EPlab);
    runidBCk = [runidBC, num2str(k)];
    prob = coco_prob();
    prob = cocoSet(obj, prob);
    prob = coco_set(prob, 'ode', 'autonomous', false);
    prob = coco_set(prob, 'ode', 'vectorized', true);
    % zero problem
    prob = ode_po2po(prob, '', runidBPk, EPlab);
    [fdata, uidx] = coco_get_func_data(prob, 'po.orb.coll', 'data', 'uidx');
    maps = fdata.coll_seg.maps;
    prob = coco_add_func(prob, 'OmegaT', @OmegaT, @OmegaT_du, omData, 'zero',...
        'uidx', [uidx(maps.T_idx), uidx(maps.p_idx(1))]);
    data  = obj_init_data(fdata,optdof);
    prob = coco_add_func(prob, 'obj', @objhan, @objhan_du, data, ...
      'inactive', 'obj', 'uidx', uidx, 'remesh', @obj_remesh);
    prob = coco_add_func(prob, 'amp', @amplitude, ampdata, 'regular', ampNames,...
        'uidx', uidx(maps.xbp_idx), 'remesh', @amplitude_remesh);  
    % adjoint
    prob = adjt_po2po(prob, '', runidBPk, EPlab);
    chart = coco_read_adjoint('OmegaT', runidBPk, EPlab, 'chart');
    [data, axidx] = coco_get_adjt_data(prob, 'po.orb.coll', 'data', 'axidx');
    opt = data.coll_opt;
    prob = coco_add_adjt(prob, 'OmegaT', 'aidx', ...
      axidx([opt.T_idx; opt.p_idx(1)]),'l0', chart.x);
    chart = coco_read_adjoint('obj', runidBPk, EPlab, 'chart');
    data = adj_obj_init_data(data,optdof);
    prob = coco_add_adjt(prob, 'obj', @adj_objhan, @adj_objhan_du, data, ...
      'd.obj', 'aidx', axidx([opt.xcn_idx; opt.T0_idx; opt.T_idx; opt.p_idx]), ...
      'l0', chart.x, 'remesh', @adj_obj_remesh);   
    prob = coco_add_event(prob, 'BC', 'BP', 'eps', '<', parRange{2}(1));
    % continuation run
    bdk = coco(prob, runidBCk, [], contBC_args{:},{[], parRange{2}, parRange{1}});
    bds{k+1} = bdk;
end

end