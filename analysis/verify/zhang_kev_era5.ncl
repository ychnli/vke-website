; Verification B: Zhang et al. (2026) KEV budget terms, run on global ERA5 pressure levels.
;
; The computational lines below are copied verbatim from
; reference/zhang_etal_2026/GRL2026_codes/ARanalysisMERRA2_IVKE1.ncl ("Time average steps"
; and the vibeta integrals). Adaptations, all to the INPUT side only:
;   * ERA5 is already on pressure levels, so the MERRA-2 model-level fields (Q0, U0, ...)
;     and the target-level fields (Q, U, ...) are the same arrays, PL is the pressure
;     level value broadcast to 4-D, and linint1_n_Wrap (model -> pressure levels) is skipped;
;   * ERA5 geopotential is already Phi (m2 s-2), so PHI = Z instead of PHI = H*g;
;   * latitude is supplied south-to-north, as in MERRA-2 (required by gradsf/uv2dvF).
; Inputs: f = addfile(global file with Q,U,V,W,Z,SP; level 1000..150 bottom-to-top, lat S->N)

g   =    9.81
pB  = 1000.
pT  =  200.

lev      = f->level
lat      = f->latitude
lon      = f->longitude
PS       = f->SP
Q0       = f->Q
U0       = f->U
V0       = f->V
OMEGA0   = f->W
PL       = conform(Q0, dble2flt(lev)*100., 1)

p        = dble2flt(lev)*100.
pB       = pB*100.
pT       = pT*100.
LAT      = dble2flt(lat)
LON      = dble2flt(lon)

; copy attributes (as in the original, which copies from DQVDTANA)
VKE1     = Q0 ; copy attributes
HAKE     = Q0 ; copy attributes
VAKE     = Q0 ; copy attributes
PEKE     = Q0 ; copy attributes
HAV      = Q0 ; copy attributes
VAV      = Q0 ; copy attributes
HCV      = Q0 ; copy attributes
VCV      = Q0 ; copy attributes

; ---- verbatim (Time average steps) ----
dUdp0    = center_finite_diff_n(U0       ,PL,False,0,1)
WdUdp0   = OMEGA0*dUdp0
dVdp0    = center_finite_diff_n(V0       ,PL,False,0,1)
WdVdp0   = OMEGA0*dVdp0
dQdp0    = center_finite_diff_n(Q0       ,PL,False,0,1)
WdQdp0   = OMEGA0*dQdp0
dQWdp0   = center_finite_diff_n(Q0*OMEGA0,PL,False,0,1)

gradKx0  = Q0
gradKy0  = Q0
gradQx0  = Q0
gradQy0  = Q0
gradsf((U0*U0+V0*V0)*0.5,gradKx0,gradKy0)
gradsf( Q0              ,gradQx0,gradQy0)
UgradQx0 = U0*gradQx0
VgradQy0 = V0*gradQy0

divQ0    = uv2dvF(Q0*U0,Q0*V0)

; (linint1_n_Wrap from model levels skipped: already on pressure levels)
Q        = Q0
U        = U0
V        = V0
K        = (U*U+V*V)*0.5
QU       = Q*U
QV       = Q*V

WdUdp    = WdUdp0
WdVdp    = WdVdp0
WdQdp    = WdQdp0
dQWdp    = dQWdp0
VAKE     = -QU*WdUdp-QV*WdVdp
VAV      = -K*WdQdp
VCV      = -K*dQWdp

PHI      = f->Z        ; adaptation: ERA5 geopotential is already Phi (was H*g)
gradPHI  = grad_latlon_cfd(PHI,LAT,LON,True,False)
gradKx   = gradKx0
gradKy   = gradKy0
UgradQx  = UgradQx0
VgradQy  = VgradQy0
divQ     = divQ0
PEKE     = -  QU*gradPHI[1]-  QV*gradPHI[0]
HAKE     = -  QU*gradKx    -  QV*gradKy
HAV      = - K*(UgradQx + VgradQy)
HCV      = - K*  divQ

n_lat    = dimsizes(lat)
IHAKE    = vibeta(p,HAKE (time|:,latitude|:,longitude|:,level|:),1,PS,pB,pT) /g
IVAKE    = vibeta(p,VAKE (time|:,latitude|:,longitude|:,level|:),1,PS,pB,pT) /g
IPEKE    = IHAKE
IPEKE(:,1:n_lat-2,:) = vibeta(p,PEKE (time|:,latitude|1:n_lat-2,longitude|:,level|:),1,PS(:,1:n_lat-2,:),pB,pT) /g
IHAV     = vibeta(p,HAV  (time|:,latitude|:,longitude|:,level|:),1,PS,pB,pT) /g
IVAV     = vibeta(p,VAV  (time|:,latitude|:,longitude|:,level|:),1,PS,pB,pT) /g
IHCV     = vibeta(p,HCV  (time|:,latitude|:,longitude|:,level|:),1,PS,pB,pT) /g
IVCV     = vibeta(p,VCV  (time|:,latitude|:,longitude|:,level|:),1,PS,pB,pT) /g
VKE1     = Q*K
IVKE1    = vibeta(p,VKE1(time|:,latitude|:,longitude|:,level|:),1,PS,pB,pT) /g

system("rm -f " + outfile)
o = addfile(outfile, "c")
o->IHAKE = IHAKE
o->IVAKE = IVAKE
o->IPEKE = IPEKE
o->IHAV  = IHAV
o->IVAV  = IVAV
o->IHCV  = IHCV
o->IVCV  = IVCV
o->IVKE1 = IVKE1
