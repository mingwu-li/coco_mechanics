function [prob, status, xtr] = obj_remesh(prob, data, chart, old_u, old_V) %#ok<INUSD>

[fdata, uidx] = coco_get_func_data(prob, 'po.orb.coll', 'data', 'uidx');
data.coll_seg = fdata.coll_seg;

xtr    = [];
prob   = coco_change_func(prob, data, 'uidx', uidx);
status = 'success';

end