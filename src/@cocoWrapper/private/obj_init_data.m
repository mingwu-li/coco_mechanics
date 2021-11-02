function data = obj_init_data(fdata,optdof)

data.coll_seg = fdata.coll_seg;
data.ghan     = @ghan;
data.ghan_dx  = @ghan_dx;
data.optdof   = optdof;

data = coco_func_data(data);

end
