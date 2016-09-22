% test for ComputeJacobianN for dim = 5;
% 
function J = jacob(a,b,c,d,e)
syms v w x y z
T = jacobian([cos(v), sin(v)*cos(w),sin(v)*sin(w)*cos(x),sin(v)*sin(w)*sin(x)*cos(y),sin(v)*sin(w)*sin(x)*sin(y)*cos(z),sin(v)*sin(w)*sin(x)*sin(y)*sin(z)],[v,w,x,y,z]);
J = double(subs(T,{v,w,x,y,z},{a,b,c,d,e}));
 