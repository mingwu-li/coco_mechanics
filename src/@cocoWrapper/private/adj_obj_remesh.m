function [prob, status, xtr, ftr] = adj_obj_remesh(prob, data, chart, lb, Vlb)  %#ok<INUSL>

[fdata, axidx] = coco_get_adjt_data(prob, 'po.orb.coll', 'data', 'axidx');
data = adj_obj_init_data(fdata,data.optdof);
fopt = fdata.coll_opt;

aidx = axidx([fopt.xcn_idx; fopt.T0_idx; fopt.T_idx; fopt.p_idx]);

xtr    = [];
ftr    = 1;
prob   = coco_change_adjt(prob, data, 'aidx', aidx, 'l0', lb, 'vecs', Vlb);
status = 'success';

end
