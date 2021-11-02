function y = dFextdp(obj,t,p)
% FEXT This function returns the derivative of external force F with
% respect to problem parameters p. EOM: M\ddot{u}+N(u,\dot{u})=F(t,p). 
% This function will be used in forward toolbox for finding periodic orbits
%
% Y = DFEXTDP(OBJ,T,P)
%
% obj: coco object with dynamical system included
% t:   time
% p:   problem parameter (omega,epsilon)
% 
% See also: FEXT

om = p(1);
ep = p(2);
fext_coeffs = obj.system.fext.coeffs(:,1);
fext_harm_ep = cos(obj.system.fext.kappas(1)*om*t);
fext_harm_om = -obj.system.fext.kappas(1)*t*ep*sin(obj.system.fext.kappas(1)*om*t);
y = 2*fext_coeffs*[fext_harm_om fext_harm_ep];

% Jp = coco_ezDFDP('f(x,p)', @obj.Fext, t, p);

end