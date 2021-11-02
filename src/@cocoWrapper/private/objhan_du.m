function [data, J] = objhan_du(prob, data, u) %#ok<INUSL>

pr   = data.pr;
maps = pr.coll_seg.maps;
mesh = pr.coll_seg.mesh;

T = u(maps.T_idx);
x = u(maps.xbp_idx);
p = u(maps.p_idx);

xcn = reshape(maps.W*x, maps.x_shp);
% gcn = pr.ghan(xcn,data.optdof);
% gcn = mesh.gka.*gcn;

gdxcn = pr.ghan_dx(xcn,data.optdof);
gdxcn = mesh.gdxka.*gdxcn;
gdxcn = sparse(maps.gdxrows, maps.gdxcols, gdxcn(:));

J_xbp = (0.5/maps.NTST)*mesh.gwt*gdxcn*maps.W;
J_T0  = 0;
J_T   = 0;
J_p   = zeros(1,numel(p));

J = [J_xbp J_T0 J_T J_p];

% [data, Jd] = coco_ezDFDX('f(o,d,x)', prob, data, @objhan, u);
% max(max(abs(J-Jd)))
end