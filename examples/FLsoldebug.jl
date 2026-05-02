using Revise, FL2D_small

import FL2D_small.FLdata as FLdata

#----------------------------
s = 0.5

d = FL2D_small.disc(b = [2, 2, 2, 2, 2])

dp = FL2D_small.domprop(10, 0.1, 0.01, d)

f!, uex, _ = FLdata.makediscfuex(2, s);

b = FL2D_small.bvec(d, dp, s, f!);

FL2D_small.plotfunc(dp,d,b)

IV = FL2D_small.compress_vars(d, dp.N, s, 4, 1, 0.01; matrix_form=false);

vv = IV.IV1.Fsvec;

FL2D_small.Fsv!(IV, d, dp, s, 4);

FL2D_small.plotfunc(dp,d,IV.IV1.Fsvec)

#---------------------------------------

IntSFsG = FL2D_small.precompsL(d, dp, s, 4);

IntSFs = FL2D_small.precompsFs(d, dp, s, 4);

ERR = abs.(IntSFs .- vec(IntSFsG[1, :]));

maximum(ERR)

all(isfinite, IntSFsG)

any(isnan, IntSFsG)

all(isfinite, IntSFs)

any(isnan, IntSFs)

all(isfinite, ERR)

any(isnan, ERR)