; Integrate random columns with NCL vibeta for comparison with arbudget.grid.vibeta
f = addfile("data/vibeta_in.nc", "r")
p = f->p
x = f->x
ps = f->ps
out = vibeta(p, x, 1, ps, 100000., 20000.)
system("rm -f data/vibeta_out.nc")
g = addfile("data/vibeta_out.nc", "c")
g->vint = out
