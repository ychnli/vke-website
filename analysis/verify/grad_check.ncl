f = addfile("data/grad_in.nc", "r")
z = f->z
lat = f->lat
lon = f->lon
p = f->p
gr = grad_latlon_cfd(z, lat, lon, True, False)
dzdp = center_finite_diff_n(z, p, False, 0, 0)
system("rm -f data/grad_out.nc")
g = addfile("data/grad_out.nc", "c")
g->gy = gr[0]
g->gx = gr[1]
g->dzdp = dzdp
