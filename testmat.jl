using MIRT
using MATLAB
if !@isdefined(irtdir)
	ENV["MATLAB_ROOT"] = "/Applications/matlab"

	irtdir = "/Users/jasonhu/Documents/MATLAB/phd/irt"
	tmp = "addpath('$irtdir')"
	eval_string(tmp)
	mat"setup"
end
mat"sg = sino_geom('ge1', 'units', 'mm', 'strip_width', 'd', 'down', 1);"

mat"ig = image_geom('nx', 420, 'dx', 500/512, 'down', 1);"

mat"ig.mask = ig.circ > 0;"
mat"A = Gtomo2_dscmex(sg, ig);"
