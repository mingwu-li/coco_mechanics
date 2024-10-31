function [mass,damp,stiff,fnl,fnl_dx,fext] = build_model_func()

n   = 2;
om1 = 2;
om2 = 4.5;
D1  = 0.01;
D2  = 0.2;
f1  = 0.02;
mass = eye(n,n);
damp = [2*D1*om1, 0;
     0, 2*D2*om2];
stiff = [om1^2 0;0 om2^2];
fnl   = @(x) spring_fnl(x,om1,om2); % function handle for nonlinearity
fnl_dx = @(x) spring_dfnldx(x,om1,om2);     % function handle for jacobian of nonlinearity
fext = [f1;0];

end

function y = spring_fnl(x,om1,om2)
% setup
tmp1 = om1^2;
tmp2 = om2^2;
tmp3 = (tmp1+tmp2)/2;
q1   = x(1); 
q2   = x(2);
% fnl
y    = zeros(2,1);
y(1) = tmp1/2*(3*q1^2+q2^2)+tmp2*q1*q2+tmp3*q1*(q1^2+q2^2);
y(2) = tmp2/2*(3*q2^2+q1^2)+tmp1*q1*q2+tmp3*q2*(q1^2+q2^2);
end

function J = spring_dfnldx(x,om1,om2)
% setup
tmp1 = om1^2;
tmp2 = om2^2;
tmp3 = (tmp1+tmp2)/2;
q1   = x(1); 
q2   = x(2);
% fnl
J    = zeros(2,2);
J(1,1) = tmp1*3*q1+tmp2*q2+3*tmp3*q1^2;
J(1,2) = tmp1*q2+tmp2*q1+2*tmp3*q1*q2;
J(2,1) = tmp2*q1+tmp1*q2+2*tmp3*q1*q2;
J(2,2) = tmp2*3*q2+tmp1*q1+3*tmp3*q2^2;
end