run_id
test
istart   nwx      nwy      nwz     ibc     HALL_ON    gradP_ON  FLR_ON
0        2        2        1       4       1          1         1
xl       yl       zl       B00     angolo  U_0        L_U       Delta_N
150.0    15.0     1.0      1.0     0.02    0.5        3.0       0.25
dt       rk_ord   tmax     tt_w    tx_w    nx_lambda  grdlmbd
0.0025   3        3.01     0.005   100.0   100        10
omega_x  omega_y  omega_z  eta     ampl    alpha
0.5      0.6      0.66     1.e-4   1.e-3   1.0


