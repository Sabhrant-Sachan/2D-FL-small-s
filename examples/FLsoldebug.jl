using Revise, FL2D_small

import FL2D_small.FLdata as FLdata

#----------------------------
s = 0.5

d = FL2D_small.disc(b = [1, 1, 1, 1, 1])

dp = FL2D_small.domprop(12, 0.1, 0.01, d)

f!, uex, _ = FLdata.makediscfuex(0, s);

# 8.263 s (109 allocations: 1.25 MiB)
b = FL2D_small.bvec(d, dp, f!);

FL2D_small.@btime FL2D_small.bvec($d, $dp, $f!);

FL2D_small.plotfunc(dp,d,b)

IV = FL2D_small.compress_vars(d, dp.N, s, 4, 1, 0.01; matrix_form=false);

vv = IV.IV1.Fsvec;

FL2D_small.Fsv!(IV, d, dp, s, 4);

#  95.616 s (270 allocations: 19.67 MiB)
FL2D_small.@btime FL2D_small.Fsv!($IV, $d, $dp, $s, 4);

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

#-----------------------------------------
# Convergence analysis of precomps 

s, p = 0.1, 4

δ, δclsbd = 0.1, 0.01

AN = 5

d = FL2D_small.disc(b=[1, 1, 1, 1, 1], L1=0.8, L2=0.8);

dp = FL2D_small.domprop(AN, δ, δclsbd, d);

Anₚᵣ = [32, 64, 128, 256, 512];

IntS_ex = FL2D_small.precomps(d, dp, s, p; n=1024);

for i in 1:5

    IntS = FL2D_small.precomps(d, dp, s, p; n=Anₚᵣ[i])

    Err = maximum(abs.(IntS .-  IntS_ex))

    display(Err)

end
