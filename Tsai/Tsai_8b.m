% f = Tsai_8b(params, params_const, xw, yw, zw, X, Y)
%
% **********************************************************************************************
% *******         Calibrating a Camera Using a Monoview Coplanar Set of Points           *******
% **********************************************************************************************
%                              6/2004   Simon Wan 
%                              //2006-03-04 如有疑问：simonwan1980@gmail.com (因为已从哈工大毕业，此地址已作废simonwan1980@hit.edu.cn)
%
% Note:    This is not called directly but as a function handle from the "fminsearch "
%
function f = Tsai_8b(params, params_const, xw, yw, zw, X, Y)
% unpack the params
    f  = params(1);
	Tz = params(2);
	k1 = params(3);
% unpack the params_const
	r4 = params_const(1);
	r5 = params_const(2);
	r6 = params_const(3);
	r7 = params_const(4);
	r8 = params_const(5);
	r9 = params_const(6);
	dx = params_const(7);
	dy = params_const(8);
	sx = params_const(9);
	Ty = params_const(10);
    
    rsq = (dx*X).^2 + (dy*Y).^2;
	res = (dy*Y).*(1+k1*rsq).*(r7*xw+r8*yw+r9*zw+Tz) - f*(r4*xw+r5*yw+r6*zw+Ty);
    f = norm(res, 2);