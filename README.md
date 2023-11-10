X:
DT_1,   mu_r1-,           mu_v1-,             DV_1,           K_1, (1) (3) (3) (3) (18)
DT_2,   mu_r2-,           mu_v2-,             DV_2,           K_2,
...
0,      mu_r(nSeg + 1)-,  mu_v(nSeg + 1)-,    DV_(nSeg + 1),  K_(nSeg + 1).
XnLeg

G:
r_1,        v_1, (3) (3)
r_2,        v_2.
...
r_nSeg,     v_nSeg,
r_f,        v_f,

(sigma_rf,   sigma_vf,) (3) (3)

DV_1,
DV_2,
...
DV_(nSeg + 1),

(ToF_Leg)
XnLeg

norm(v_inf+_)
r_p>=r_min
XnFB

(ToF)
