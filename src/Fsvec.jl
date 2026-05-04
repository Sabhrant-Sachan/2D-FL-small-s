# Compute  F_s[1](x) = ∫_Ω d(y)^s / |x-y|^(2s) dy
function Fsv!(IV::IVT, d::D, dp::domprop, s::Float64, p::Int; n::Int=256) where {D<:abstractdomain,IVT}

    (; IV1, IVr, IVbdth, IVt, IVbt1, IVbt2) = IV

    (; N, Np, M, Mbd, Fsvec) = IV1

    (; fwr, zx, zt, zy, Df) = IVr

    (; zx2, zy2) = IVbdth
    
    (; Zx, Zy, DJ, KIr) = IVt

    (; Zx₂, Zy₂, DJ₂, KIbd) = IVbt1

    (; mfw) = IVbt2

    # ----------------------------
    # main loop over target points
    # ----------------------------
    Lp = M * Np
    Lₚₙ = Lp + Mbd * N

    IntS = precompsFs(d, dp, s, p; n=n)

    v = Fsvec

    fill!(v, 0.0)

    @inbounds for k in 1:M
        isbdflag = (k in d.kd)

        if isbdflag
            # Zx₂, Zy₂: images of Chebyshev grid (zx2,zy2) on patch k
            mapxy_Dmap!(Zx₂, Zy₂, DJ₂, d, zx2, zy2, k) # nr×nr Zx₂, Zy₂, DJ₂

            hc = d.pths[k].ck1 - d.pths[k].ck0

            Dhc = hc^s
        else
            # Zx, Zy: images of Chebyshev grid (zx,zy) on patch k
            mapxy_Dmap!(Zx, Zy, DJ, d, zx, zy, k) # nr×nr Zx, Zy, DJ

            dfunc!(Df, d, k, zt, s)

            @. DJ = DJ * Df
        end

        @inbounds for row in 1:Lₚₙ

            if row <= Lp
                # Interior target row.
                ℓ = cld(row, Np)
                j = row - (ℓ - 1) * Np

                if k == ℓ
                    # Singular interior target on patch k.
                    col = dp.pthgo[k] + j - 1
                    v[row] += IntS[col]
                    continue
                end
            else
                # Boundary target row.
                k₀ = cld(row - Lp, N)
                ℓ = d.kd[k₀]

                if k == ℓ
                    # Singular boundary target on boundary patch k.
                    col = dp.pthgo[M+1] + row - Lp - 1 
                    v[row] += IntS[col]
                    continue
                end
            end

            Ikey = packkey(row, k)
            col = get(dp.hmap, Ikey, 0)
            if col != 0
                # ---------- Near-singular patch case ----------
                v[row] += IntS[col]
            else
                # ---------- Regular patch case ----------
                x1 = dp.tgtpts[1, row]
                x2 = dp.tgtpts[2, row]

                if isbdflag
                    @. KIbd = ((x1 - Zx₂)^2 + (x2 - Zy₂)^2)^(-s)
                    @. KIbd = KIbd * DJ₂
                    v[row] += Dhc * dot(mfw, KIbd, fwr)

                else
                    @. KIr = ((x1 - Zx)^2 + (x2 - Zy)^2)^(-s)
                    @. KIr = KIr * DJ
                    # Computes fwr' * KIr * fwr
                    v[row] += dot(fwr, KIr, fwr)
                end

            end

        end

    end

    return nothing
end

struct ThetaQuad
    nt::Int
    thet1::Vector{Float64}
    thet2::Vector{Float64}
    Sθ₁::Vector{Float64}
    Cθ₁::Vector{Float64}
    Sθ₂::Vector{Float64}
    Cθ₂::Vector{Float64}
    fwmt1::Vector{Float64}
    fwmt2::Vector{Float64}
    fwmt3::Vector{Float64}
    r₁::Matrix{Float64}
    r₂::Matrix{Float64}
    r₃::Matrix{Float64}
    rs₁::Matrix{Float64}
    rs₂::Matrix{Float64}
    rs₃::Matrix{Float64}
end

function makeThetaQuad(nt::Int, nr::Int, q::Int, s::Float64, wz1, wz2)
    zt  = Vector{Float64}(undef, nt)
    z2t = Vector{Float64}(undef, nt)

    # Fejér nodes and z2t
    @inbounds for i in 1:nt
        c₁ = π * (2 * i - 1) / (2 * nt)
        zt[i]  = cos(c₁)
        z2t[i] = sin(c₁ / 2)^2
    end

    # angular Fejér weights
    fwt = getF1W(nt)

    # polynomial change of variables
    thet1 = Vector{Float64}(undef, nt)
    thet2 = Vector{Float64}(undef, nt)
    wfunc!(thet1, q, z2t; α=-1.0, β=π/4)
    wfunc!(thet2, q, zt;  α=1.0,  β=π/8)

    Sθ₁ = similar(thet1); Cθ₁ = similar(thet1); Tθ₁ = similar(thet1)
    Sθ₂ = similar(thet1); Cθ₂ = similar(thet1); Tθ₂ = similar(thet1)

    @inbounds for i in 1:nt
        si₁, co₁ = sincos(thet1[i])
        si₂, co₂ = sincos(thet2[i])
        Sθ₁[i] = si₁;  Cθ₁[i] = co₁;  Tθ₁[i] = si₁ / co₁
        Sθ₂[i] = si₂;  Cθ₂[i] = co₂;  Tθ₂[i] = si₂ / co₂
    end

    # fwmt* (angular weights)
    fwmt1 = Vector{Float64}(undef, nt)
    fwmt2 = Vector{Float64}(undef, nt)
    fwmt3 = Vector{Float64}(undef, nt)

    tmp2  = Vector{Float64}(undef, nt)
    tmp3  = Vector{Float64}(undef, nt)
    dwfunc!(tmp2, q, z2t)
    dwfunc!(tmp3, q, zt)

    @inbounds for i in 1:nt
        fwmt1[i] = π * fwt[i] * tmp2[i] * (Cθ₁[i])^(2*(s - 1))   / 16
        fwmt2[i] = π * fwt[i] * tmp2[i] * (2*Cθ₁[i])^(2*(s - 1)) / 8
        fwmt3[i] = π * fwt[i] * tmp3[i] * (2*Cθ₂[i])^(2*(s - 1)) / 8
    end

    #Radial parts of singular integrals
    r₁ = Matrix{Float64}(undef, nt, nr) #sec(thet1) * w(-z2r)'
    r₂ = Matrix{Float64}(undef, nt, nr) #sec(thet1) * w(zr)' / 2
    r₃ = Matrix{Float64}(undef, nt, nr) #sec(thet2) * w(zr)' / 2
    rs₁= Matrix{Float64}(undef, nt, nr) #tan(thet1) * w(-z2r)'
    rs₂= Matrix{Float64}(undef, nt, nr) #tan(thet1) * w(zr)'
    rs₃= Matrix{Float64}(undef, nt, nr) #tan(thet2) * w(zr)'

    @inbounds for i in 1:nt
        c1 = Cθ₁[i]; t1 = Tθ₁[i]
        c2 = Cθ₂[i]; t2 = Tθ₂[i]
        @inbounds for j in 1:nr
            r₁[i,j]  = wz2[j] / c1
            r₂[i,j]  = wz1[j] / (2c1)
            r₃[i,j]  = wz1[j] / (2c2)
            rs₁[i,j] = t1 * wz2[j]
            rs₂[i,j] = t1 * wz1[j]
            rs₃[i,j] = t2 * wz1[j]
        end
    end

    return ThetaQuad(nt, thet1, thet2,
                     Sθ₁, Cθ₁, Sθ₂, Cθ₂,
                     fwmt1, fwmt2, fwmt3,
                     r₁, r₂, r₃, rs₁, rs₂, rs₃)
end

function precompsFs(d::abstractdomain, dp::domprop, s::Float64, p::Int
    ; n::Int=128)::Vector{Float64}

    #Bookkeeping
    M = d.Npat  # number of patches
    N = dp.N    # number of nodes per patch per Axis
    Np = N * N  # number of nodes per patch

    Mbd= length(d.kd) #Number of bonundary patches
    nbd = Mbd * N
    # Col index where the bd near–singular block starts in `prepts`
    bdnsp = dp.pthgo[M+1] + nbd

    # zp1 = zp+1, zp2 = xp-1 where zp = cos(pi*(2j-1)/(2N))
    zp  = Vector{Float64}(undef, N)
    zp1 = Vector{Float64}(undef, N)
    zp2 = Vector{Float64}(undef, N)

    @inbounds for j in 1:N
        zp[j]  = cospi((2 * j - 1) / (2N))
        zp1[j] = 2 * cospi((2 * j - 1) / (4N))^2
        zp2[j] =-2 * sinpi((2 * j - 1) / (4N))^2
    end

    #Parameter q for another polynomial change of variables is theta. 
    # It is usually taken as 4
    q = 4

    #Define nodes in radial and weights
    nr = n     
    fwr = getF1W(nr)

    zr = Vector{Float64}(undef, nr)
    z2r = Vector{Float64}(undef, nr)

    @inbounds for i in 1:nr
        c₁ = π * (2 * i - 1) / (2 * nr)
        zr[i] = cos(c₁)
        z2r[i] = sin(c₁ / 2)^2
    end

    #Precomputing modified weights
    fwmr1 = Vector{Float64}(undef, nr)
    fwmr2 = Vector{Float64}(undef, nr)
    fwmr3 = Vector{Float64}(undef, nr)

    #To store w(zr) and w(-z2r)
    wz1 = Vector{Float64}(undef, nr)
    wz2 = similar(wz1)
    wz3 = similar(wz1)

    wfunc!(wz1, p, zr)
    wfunc!(wz2, p, z2r; α=-1.0)
    wfunc!(wz3, p, zr; α=-1.0)

    qw1func!(fwmr1, p, z2r, s)
    qw1func!(fwmr3, p, zr, s; α=-1.0)

    @. fwmr1 = fwr * fwmr1
    @. fwmr2 = fwr * fwmr3
    @. fwmr3 = fwr * fwmr3

    # Storage of temporary variables inside the for loop
    # These are independent on nt

    yr  = Vector{Float64}(undef, nr)

    dfy = Vector{Float64}(undef, nr)

    #This is nt is the finest one
    nt = 2*n
    dfYc= Matrix{Float64}(undef, nt, nr)

    QuadT = makeThetaQuad(nt, nr, q, s, wz1, wz2)

    Sθ₁ = QuadT.Sθ₁
    Cθ₁ = QuadT.Cθ₁
    Sθ₂ = QuadT.Sθ₂
    Cθ₂ = QuadT.Cθ₂
    r₁  = QuadT.r₁
    r₂  = QuadT.r₂
    r₃  = QuadT.r₃
    rs₁ = QuadT.rs₁
    rs₂ = QuadT.rs₂
    rs₃ = QuadT.rs₃
    fwmt1 = QuadT.fwmt1
    fwmt2 = QuadT.fwmt2
    fwmt3 = QuadT.fwmt3

    d1   = Vector{Float64}(undef, nt)
    d2   = similar(d1)
    y1   = Matrix{Float64}(undef, nt, nr)
    y1tmp= similar(y1)
    y2   = similar(y1)
    dfY  = Matrix{Float64}(undef, nt, nr)
    J1   = Matrix{Float64}(undef, nt, nr)
    DIF  = Matrix{Float64}(undef, nt, nr)
    Zx   = Matrix{Float64}(undef, nt, nr)
    Zy   = Matrix{Float64}(undef, nt, nr)
    DJ   = Matrix{Float64}(undef, nt, nr)

    @inbounds for i in 1:nt
        cθ = sqrt(2) * (cos(QuadT.thet2[i] + π / 4) / QuadT.Cθ₂[i])
        @inbounds for j in 1:nr
            dfYc[i, j] = (wz3[j] + cθ * wz1[j])^s
        end
    end

    J2  = Vector{Float64}(undef, nr)
    G   = Vector{Float64}(undef, nr)

    # This is the case when s<0.5
    Lₚ = size(dp.prepts,2) # > dp.pthgo[M+1] - 1

    # First we will go over all the points in the interioir
    Dhc = Vector{Float64}(undef,M)

    @inbounds for k = 1:M
        Dhc[k] = d.pths[k].ck1 - d.pths[k].ck0
    end

    # knbd are patches which are not the boundary patches
    # and d.kd are patches are touching the boundary
    knbd = setdiff(collect(1:M), d.kd)

    IntS = zeros(Float64, Lₚ)

    #--------------------------------------------
    #------------Singular Integration------------
    #--Singular Integration of Boundary patches--
    @inbounds for j in 1:Np
        # j in the linear index of the Chebyshev
        # point. We will compute precomps of all 
        # points τₖ(x₁,x₂) at once by now varying
        # the patches k. This removes some 
        # repeated calculations and thereby saves 
        # time!
        qq, rr = divrem(j - 1, N)
        qq = qq + 1
        rr = rr + 1
        x1, x2 = zp[rr], zp[qq]
        x2p = zp1[qq]
        x1m = zp2[rr]
        x2m = zp2[qq]
        x̃₃ = -x1m * x2p
        x̃₄ = x1m * x2m

        # Singular and k is a boundary patch!
        #-----------------3rd part-----------------
        #------------1ˢᵗ half of part 3------------
        @. d1 = x1m * Cθ₁
        @. d2 = x2p * Sθ₁
        @. yr = x1 - x1m * wz1 / 2
        y1 .= yr'
        @. y2 = x2 - x2p * rs₂ / 2

        @. dfy = wz3 ^ s

        @inbounds for k in d.kd
            diff_rmap!(DIF, Zx, Zy, DJ, d,
                x1, x2, y1, y2, r₂, d1, d2, k)

            @. J1 = abs(DIF)^(-2 * s) * DJ

            c₁ = (Dhc[k] * (-x1m / 4.0))^s

            @. J2 = c₁ * dfy * fwmr2

            mul!(G, J1', fwmt2)

            IntS[dp.pthgo[k]+j-1] += x̃₃ * dot(G, J2)
        end

        # Singular and k is a boundary patch!
        #------------2ⁿᵈ half of part 3------------
        @. d1 = x1m * Sθ₂
        @. d2 = x2p * Cθ₂
        @. yr = x2 - x2p * wz1 / 2
        @. y1 = x1 - x1m * rs₃ / 2
        y2 .= yr'

        @inbounds for k in d.kd
            diff_rmap!(DIF, Zx, Zy, DJ, d,
                x1, x2, y1, y2, r₃, d1, d2, k)

            c₁ = (Dhc[k] * (-x1m / 4.0))^s

            @. J1 = c₁ * abs(DIF)^(-2 * s) * DJ * dfYc

            mul!(G, J1', fwmt3)

            IntS[dp.pthgo[k]+j-1] += x̃₃ * dot(G, fwmr3)
        end

        # Singular and k is a boundary patch!
        #-----------------4th part-----------------
        #------------1ˢᵗ half of part 4------------
        @. d1 = x1m * Cθ₁
        @. d2 = x2m * Sθ₁
        @. yr = x1 - x1m * wz1 / 2
        y1 .= yr'
        @. y2 = x2 - x2m * rs₂ / 2

        @inbounds for k in d.kd
            diff_rmap!(DIF, Zx, Zy, DJ, d,
                x1, x2, y1, y2, r₂, d1, d2, k)

            @. J1 = abs(DIF)^(-2 * s) * DJ

            c₁ = (Dhc[k] * (-x1m / 4.0))^s

            @. J2 = c₁ * dfy * fwmr2

            mul!(G, J1', fwmt2)
            
            IntS[dp.pthgo[k]+j-1] += x̃₄ * dot(G, J2)
        end

        # Singular and k is a boundary patch!
        #------------2ⁿᵈ half of part 4------------
        @. d1 = x1m * Sθ₂
        @. d2 = x2m * Cθ₂
        @. yr = x2 - x2m * wz1 / 2
        @. y1 = x1 - x1m * rs₃ / 2
        y2 .= yr'

        @inbounds for k in d.kd
            diff_rmap!(DIF, Zx, Zy, DJ, d,
                x1, x2, y1, y2, r₃, d1, d2, k)

            c₁ = (Dhc[k] * (-x1m / 4.0))^s

            @. J1 = c₁ * abs(DIF)^(-2 * s) * DJ * dfYc

            mul!(G, J1', fwmt3)

            IntS[dp.pthgo[k]+j-1] += x̃₄ * dot(G, fwmr3)

        end

    end

    Lₚₛ = dp.pthgo[M+1] - 1
    #-Singular Integration For pts on bonudary-
    #These points are ofcourse always on bd patches
    @inbounds for j in 1:N
        # j in the linear index of the Chebyshev
        # point. We will compute precomps of all 
        # points τₖ(x₁,x₂) at once by now varying
        # the patches k.

        x2 = zp[j]
        x2p = zp1[j]
        x2m = zp2[j]
        x̃₁ = 2.0 * x2p
        x̃₂ = -2.0 * x2m

        #Singular integration is sum of four parts
        #I1,I2,I3,I4, they are all N*N matrices
        #-----------------1st part-----------------
        #------------1ˢᵗ half of part 1------------
        @. d1 = 2.0 * Cθ₁
        @. d2 = x2p * Sθ₁
        @. yr = 1.0 - 2.0 * wz2
        y1 .= yr'
        @. y2 = x2 - x2p * rs₁

        @. yr = 2.0 * wz2

        @inbounds for ℓ in 1:Mbd

            k = d.kd[ℓ]

            diff_rmap!(DIF, Zx, Zy, DJ, d,
                1.0, x2, y1, y2, r₁, d1, d2, k)

            @. J1 = abs(DIF)^(-2 * s) * DJ

            dfunc!(dfy, d, k, yr, s)

            @. J2 = dfy * fwmr1

            mul!(G, J1', fwmt1)

            IntS[Lₚₛ+(ℓ-1)*N+j] += x̃₁ * dot(G, J2)
        end

        #------------2ⁿᵈ half of part 1------------

        @. d1 = 2.0 * Sθ₁
        @. d2 = x2p * Cθ₁
        @. yr = x2 - x2p * wz2
        @. y1 = 1.0 - 2.0 * rs₁
        y2 .= yr'
        @. y1tmp = 2.0 * rs₁

        @inbounds for ℓ in 1:Mbd

            k = d.kd[ℓ]

            diff_rmap!(DIF, Zx, Zy, DJ, d,
                1.0, x2, y1, y2, r₁, d1, d2, k)

            dfunc!(dfY, d, k, y1tmp, s)

            @. J1 = abs(DIF)^(-2 * s) * DJ * dfY

            mul!(G, J1', fwmt1)

            IntS[Lₚₛ+(ℓ-1)*N+j] += x̃₁ * dot(G, fwmr1)
        end

        #-----------------2nd part-----------------
        #------------1ˢᵗ half of part 2------------
        @. d1 = 2.0 * Cθ₁
        @. d2 = x2m * Sθ₁
        @. yr = 1.0 - 2.0 * wz2
        y1 .= yr'
        @. y2 = x2 - x2m * rs₁
        @. yr = 2.0 * wz2

        @inbounds for ℓ in 1:Mbd

            k = d.kd[ℓ]

            diff_rmap!(DIF, Zx, Zy, DJ, d,
                1.0, x2, y1, y2, r₁, d1, d2, k)

            @. J1 = abs(DIF)^(-2 * s) * DJ

            dfunc!(dfy, d, k, yr, s)

            @. J2 = dfy * fwmr1

            mul!(G, J1', fwmt1)

            IntS[Lₚₛ+(ℓ-1)*N+j] += x̃₂ * dot(G, J2)
        end

        #------------2ⁿᵈ half of part 2------------

        @. d1 = 2.0 * Sθ₁
        @. d2 = x2m * Cθ₁
        @. yr = x2 - x2m * wz2
        @. y1 = 1.0 - 2.0 * rs₁
        y2 .= yr'
        @. y1tmp = 2.0 * rs₁

        @inbounds for ℓ in 1:Mbd

            k = d.kd[ℓ]

            diff_rmap!(DIF, Zx, Zy, DJ, d,
                1.0, x2, y1, y2, r₁, d1, d2, k)

            dfunc!(dfY, d, k, y1tmp, s)

            @. J1 = abs(DIF)^(-2 * s) * DJ * dfY

            mul!(G, J1', fwmt1)

            IntS[Lₚₛ+(ℓ-1)*N+j] += x̃₂ * dot(G, fwmr1)
        end

    end


    #-Singular Integration of Interior patches-
    @inbounds for j in 1:Np
        # j in the linear index of the Chebyshev
        # point. We will compute precomps of all 
        # points τₖ(x₁,x₂) at once by now varying
        # the patches k. This removes some 
        # repeated calculations and thereby saves 
        # time!
        qq, rr = divrem(j - 1, N) .+ 1

        x1, x2 = zp[rr], zp[qq]
        x1p, x2p = zp1[rr], zp1[qq]
        x1m, x2m = zp2[rr], zp2[qq]
        x̃₁ = x1p * x2p
        x̃₂ = -x1p * x2m
        x̃₃ = -x1m * x2p
        x̃₄ = x1m * x2m
       
        #Singular integration is sum of four parts
        #I1,I2,I3,I4, they are all N*N matrices
        #-----------------1st part-----------------
        #------------1ˢᵗ half of part 1------------
        @. d1 = x1p * Cθ₁
        @. d2 = x2p * Sθ₁
        @. yr = x1 - x1p * wz2
        y1 .= yr'
        @. y2 = x2 - x2p * rs₁
        @. yr = 1 - yr

        @inbounds for k in 1:M

            diff_rmap!(DIF, Zx, Zy, DJ, d,
                x1, x2, y1, y2, r₁, d1, d2, k)

            @. J1 = abs(DIF)^(-2 * s) * DJ

            dfunc!(dfy, d, k, yr, s)

            @. J2 = dfy * fwmr1

            mul!(G, J1', fwmt1)

            IntS[dp.pthgo[k]+j-1] += x̃₁ * dot(G, J2)
        end

        #------------2ⁿᵈ half of part 1------------

        @. d1 = x1p * Sθ₁
        @. d2 = x2p * Cθ₁
        @. yr = x2 - x2p * wz2
        @. y1 = x1 - x1p * rs₁
        y2 .= yr'
        @. y1tmp = 1 - y1

        @inbounds for k in 1:M
            diff_rmap!(DIF, Zx, Zy, DJ, d,
                x1, x2, y1, y2, r₁, d1, d2, k)

            dfunc!(dfY, d, k, y1tmp, s)

            @. J1 = abs(DIF)^(-2 * s) * DJ * dfY

            mul!(G, J1', fwmt1)

            IntS[dp.pthgo[k]+j-1] += x̃₁ * dot(G, fwmr1)
        end

        #-----------------2nd part-----------------
        #------------1ˢᵗ half of part 2------------
        @. d1 = x1p * Cθ₁
        @. d2 = x2m * Sθ₁
        @. yr = x1 - x1p * wz2
        y1 .= yr'
        @. y2 = x2 - x2m * rs₁
        @. yr = 1 - yr

        @inbounds for k in 1:M
            diff_rmap!(DIF, Zx, Zy, DJ, d,
                x1, x2, y1, y2, r₁, d1, d2, k)

            @. J1 = abs(DIF)^(-2 * s) * DJ

            dfunc!(dfy, d, k, yr, s)

            @. J2 = dfy * fwmr1

            mul!(G, J1', fwmt1)

            IntS[dp.pthgo[k]+j-1] += x̃₂ * dot(G, J2)
        end

        #------------2ⁿᵈ half of part 2------------

        @. d1 = x1p * Sθ₁
        @. d2 = x2m * Cθ₁
        @. yr = x2 - x2m * wz2
        @. y1 = x1 - x1p * rs₁
        y2 .= yr'
        @. y1tmp = 1 - y1

        @inbounds for k in 1:M
            diff_rmap!(DIF, Zx, Zy, DJ, d,
                x1, x2, y1, y2, r₁, d1, d2, k)

            dfunc!(dfY, d, k, y1tmp, s)

            @. J1 = abs(DIF)^(-2 * s) * DJ * dfY

            mul!(G, J1', fwmt1)

            IntS[dp.pthgo[k]+j-1] += x̃₂ * dot(G, fwmr1)
        end

        #-----------------3rd part-----------------
        #------------1ˢᵗ half of part 3------------
        @. d1 = x1m * Cθ₁
        @. d2 = x2p * Sθ₁
        @. yr = x1 - x1m * wz2
        y1 .= yr'
        @. y2 = x2 - x2p * rs₁
        @. yr = 1 - yr

        @inbounds for k in knbd

            diff_rmap!(DIF, Zx, Zy, DJ, d,
                x1, x2, y1, y2, r₁, d1, d2, k)

            @. J1 = abs(DIF)^(-2 * s) * DJ

            dfunc!(dfy, d, k, yr, s)

            @. J2 = dfy * fwmr1

            mul!(G, J1', fwmt1)

            IntS[dp.pthgo[k]+j-1] += x̃₃ * dot(G, J2)
        end

        #------------2ⁿᵈ half of part 3------------

        @. d1 = x1m * Sθ₁
        @. d2 = x2p * Cθ₁
        @. yr = x2 - x2p * wz2
        @. y1 = x1 - x1m * rs₁
        y2 .= yr'

        @. y1tmp = 1 - y1

        @inbounds for k in knbd
            diff_rmap!(DIF, Zx, Zy, DJ, d,
                x1, x2, y1, y2, r₁, d1, d2, k)

            dfunc!(dfY, d, k, y1tmp, s)

            @. J1 = abs(DIF)^(-2 * s) * DJ * dfY

            mul!(G, J1', fwmt1)

            IntS[dp.pthgo[k]+j-1] += x̃₃ * dot(G, fwmr1)
        end

        #-----------------4th part-----------------
        #------------1ˢᵗ half of part 4------------
        @. d1 = x1m * Cθ₁
        @. d2 = x2m * Sθ₁
        @. yr = x1 - x1m * wz2
        y1 .= yr'
        @. y2 = x2 - x2m * rs₁
        @. yr = 1 - yr

        @inbounds for k in knbd

            diff_rmap!(DIF, Zx, Zy, DJ, d,
                x1, x2, y1, y2, r₁, d1, d2, k)

            @. J1 = abs(DIF)^(-2 * s) * DJ

            dfunc!(dfy, d, k, yr, s)

            @. J2 = dfy * fwmr1

            mul!(G, J1', fwmt1)

            IntS[dp.pthgo[k]+j-1] += x̃₄ * dot(G, J2)
        end

        #------------2ⁿᵈ half of part 4------------
        @. d1 = x1m * Sθ₁
        @. d2 = x2m * Cθ₁
        @. yr = x2 - x2m * wz2
        y2 .= yr'
        @. y1 = x1 - x1m * rs₁
        @. y1tmp = 1 - y1

        @inbounds for k in knbd
            diff_rmap!(DIF, Zx, Zy, DJ, d,
                x1, x2, y1, y2, r₁, d1, d2, k)

            dfunc!(dfY, d, k, y1tmp, s)

            @. J1 = abs(DIF)^(-2 * s) * DJ * dfY

            mul!(G, J1', fwmt1)

            IntS[dp.pthgo[k]+j-1] += x̃₄ * dot(G, fwmr1)
        end

    end

    #--------------------------------------------
    #----------Near Singular Integration---------
    # z1 = (1+z)/2, z2 = (1-z)/2 where z = cos(pi*(2j-1)/(2n)) 
    # are Chebyshev nodes in the interval [-1,1]. This is for 
    # near singular integrals 
    #n = n #Maybe n here can be fixed to 128 ? (only for small s)
    z = Vector{Float64}(undef, n)
    z1 = Vector{Float64}(undef, n)
    z2 = Vector{Float64}(undef, n)

    #Store Fejer 1st quadrature weights for near singular integrals
    fw = getF1W(n)  

    d1 = Vector{Float64}(undef, n)
    d2 = Vector{Float64}(undef, n)

    @inbounds for i in 1:n
        c₁ = π * (2 * i - 1) / (2 * n)
        z[i] = cos(c₁)
        z1[i] = cos(c₁ / 2)^2
        z2[i] = sin(c₁ / 2)^2
    end

    y1  = Vector{Float64}(undef, n)
    y1tmp=Vector{Float64}(undef, n)
    y2  = Vector{Float64}(undef, n)
    t1  = Matrix{Float64}(undef, n, n)  # meshgrid of y1/y2 
    t2  = Matrix{Float64}(undef, n, n)  # (column = y2[j], row = y1[i])
    DJ  = Matrix{Float64}(undef, n, n)  # To store the Jacobian
    DIF = Matrix{Float64}(undef, n, n)
    Zx  = Matrix{Float64}(undef, n, n)
    Zy  = Matrix{Float64}(undef, n, n)


    TNL = Vector{Float64}(undef, n)
    TNR = Vector{Float64}(undef, n)
    A   = Matrix{Float64}(undef, n, n)  
    dw  = Vector{Float64}(undef, n)
    dfy = Vector{Float64}(undef, n)

    #-------------------------------------------
    # small helpers (no allocations)
    @inline function fill_meshgrid!(T1, T2, y1, y2)
        @inbounds for j in eachindex(y2)
             for i in eachindex(y1)
                T1[i, j] = y1[i]
                T2[i, j] = y2[j]
            end
        end
        return nothing
    end
    #-------------------------------------------
    #A vector of Bool, initialized to true for all
    #indices from 1:Lₚ. They will be updated as 
    #False  for singular points. (and ofcourse the
    #points left are near singular, which are true)
    NSI = trues(Lₚ)

    @inbounds for j in 1:Np
        @inbounds for k in 1:M
            NSI[dp.pthgo[k]+j-1] = false
        end
    end

    @inbounds for j in dp.pthgo[M+1]:bdnsp-1
        NSI[j] = false
    end

    A1 = B1 = A2 = B2 = 0.0

    #Now the big for loop for near singular points
    @inbounds for i in 1:Lₚ
        if NSI[i] == true

            k = dp.prepts[2,i]

            # (prepts index) -> `ll` (column in invpts)
            ll = if i < bdnsp
                # interior near–singular to patch k
                i - k * Np
            else
                # boundary near–singulars
                i - M * Np - nbd
            end

            α₁ = dp.invpts[1, ll]
            α₂ = dp.invpts[2, ll]

            if 1 ≤ α₁
                B1 = -winv(p, (α₁ - 1) / (α₁ + 1))
            elseif α₁ ≤ -1
                A1 = -winv(p, (-1 - α₁) / (1 - α₁))
            end

            if 1 ≤ α₂
                B2 = -winv(p, (α₂ - 1) / (α₂ + 1))
            elseif α₂ ≤ -1
                A2 = -winv(p, (-1 - α₂) / (1 - α₂))
            end

            #Not a boundary patch
            if k in knbd
                if 1 ≤ α₁
                    wfunc!(d1, p, z1; α=-B1, β=α₁ + 1.0)
                    @. y1 = α₁ - d1
                    @. y1tmp = 1 - y1
                    if 1 ≤ α₂
                        wfunc!(d2, p, z1; α=-B2, β=α₂ + 1.0)
                        @. y2 = α₂ - d2

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, α₁, α₂, t1, t2, d1, d2, k)

                        dfunc!(dfy, d, k, y1tmp, s)   # dfy = dfac(k, 1-y1)

                        #Left weights: fw .* dw(B1*z1) .* dfy 
                        dwfunc!(dw, p, z1; α=B1)
                        @. TNL = fw * dw * dfy

                        # Right weights: dw(B2*z1).*fw
                        dwfunc!(dw, p, z1; α=B2)
                        @. TNR = dw * fw   # n×N scaled row-wise
                        # Middle terms together
                        @. A = abs(DIF)^(-2*s) * DJ   

                        IntS[i] = dot(TNL, A, TNR) * B1 * B2 * (α₁ + 1.0) * (α₂ + 1.0) / 4.0

                    elseif -1 < α₂ && α₂ < 1
                        # ---------- 1st half in α₂ ----------
                        wfunc!(d2, p, z1; α=-1.0, β=α₂ + 1.0)
                        @. y2 = α₂ - d2

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, α₁, α₂, t1, t2, d1, d2, k)

                        dfunc!(dfy, d, k, y1tmp, s)   # dfy = dfac(k, 1-y1)

                        # Left weights: fw .* dw(B1*z1) .* dfy
                        dwfunc!(dw, p, z1; α=B1)
                        @. TNL = fw * dw * dfy

                        # Right weights: fw .* dw(z1)
                        dwfunc!(dw, p, z1)
                        @. TNR = fw * dw

                        @. A = abs(DIF)^(-2 * s) * DJ

                        I₁ = dot(TNL, A, TNR)

                        # ---------- 2nd half in α₂ ----------
                        wfunc!(d2, p, z2; α=-1.0, β=α₂ - 1.0)
                        @. y2 = α₂ - d2

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, α₁, α₂, t1, t2, d1, d2, k)

                        # Left weights are unchanged.
                        # Right weights: fw .* dw(z2)
                        dwfunc!(dw, p, z2)
                        @. TNR = fw * dw

                        @. A = abs(DIF)^(-2 * s) * DJ

                        I₂ = dot(TNL, A, TNR)

                        IntS[i] = B1 * (α₁ + 1.0) * ((α₂ + 1.0) * I₁ + (1.0 - α₂) * I₂) / 4.0

                    elseif α₂ ≤ -1
                        wfunc!(d2, p, z2; α=-A2, β=α₂ - 1.0)
                        @. y2 = α₂ - d2

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, α₁, α₂, t1, t2, d1, d2, k)

                        dfunc!(dfy, d, k, y1tmp, s)   # dfy = dfac(k, 1-y1)

                        # Left weights: fw .* dw(B1*z1) .* dfy
                        dwfunc!(dw, p, z1; α=B1)
                        @. TNL = fw * dw * dfy

                        # Right weights: fw .* dw(A2*z2)
                        dwfunc!(dw, p, z2; α=A2)
                        @. TNR = fw * dw

                        @. A = abs(DIF)^(-2 * s) * DJ

                        IntS[i] = B1 * A2 * (α₁ + 1.0) * (1.0 - α₂) * dot(TNL, A, TNR) / 4.0
                    end
                elseif -1 < α₁ && α₁ < 1
                    if 1 ≤ α₂
                        # ------------ 1st half in α₁ ------------
                        wfunc!(d1, p, z1; α=-1.0, β=α₁ + 1.0)
                        wfunc!(d2, p, z1; α=-B2, β=α₂ + 1.0)

                        @. y1 = α₁ - d1
                        @. y2 = α₂ - d2
                        @. y1tmp = 1 - y1

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, α₁, α₂, t1, t2, d1, d2, k)

                        dfunc!(dfy, d, k, y1tmp, s)

                        # Left weights: fw .* dw(z1) .* dfy
                        dwfunc!(dw, p, z1)
                        @. TNL = fw * dw * dfy

                        # Right weights: fw .* dw(B2*z1)
                        dwfunc!(dw, p, z1; α=B2)
                        @. TNR = fw * dw

                        @. A = abs(DIF)^(-2 * s) * DJ

                        I₁ = dot(TNL, A, TNR)

                        # ------------ 2nd half in α₁ ------------
                        wfunc!(d1, p, z2; α=-1.0, β=α₁ - 1.0)
                        wfunc!(d2, p, z1; α=-B2, β=α₂ + 1.0)

                        @. y1 = α₁ - d1
                        @. y2 = α₂ - d2
                        @. y1tmp = 1 - y1

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, α₁, α₂, t1, t2, d1, d2, k)

                        dfunc!(dfy, d, k, y1tmp, s)

                        # Left weights: fw .* dw(z2) .* dfy
                        dwfunc!(dw, p, z2)
                        @. TNL = fw * dw * dfy

                        # Right weights: fw .* dw(B2*z1)
                        dwfunc!(dw, p, z1; α=B2)
                        @. TNR = fw * dw

                        @. A = abs(DIF)^(-2 * s) * DJ

                        I₂ = dot(TNL, A, TNR)

                        IntS[i] = B2 * (α₂ + 1.0) * ((α₁ + 1.0) * I₁ + (1.0 - α₁) * I₂) / 4.0

                    elseif α₂ ≤ -1
                        # ------------ 1st half in α₁ ------------
                        wfunc!(d1, p, z1; α=-1.0, β=α₁ + 1.0)
                        wfunc!(d2, p, z2; α=-A2, β=α₂ - 1.0)

                        @. y1 = α₁ - d1
                        @. y2 = α₂ - d2
                        @. y1tmp = 1 - y1

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, α₁, α₂, t1, t2, d1, d2, k)

                        dfunc!(dfy, d, k, y1tmp, s)

                        # Left weights: fw .* dw(z1) .* dfy
                        dwfunc!(dw, p, z1)
                        @. TNL = fw * dw * dfy

                        # Right weights: fw .* dw(A2*z2)
                        dwfunc!(dw, p, z2; α=A2)
                        @. TNR = fw * dw

                        @. A = abs(DIF)^(-2 * s) * DJ

                        I₁ = dot(TNL, A, TNR)

                        # ------------ 2nd half in α₁ ------------
                        wfunc!(d1, p, z2; α=-1.0, β=α₁ - 1.0)
                        wfunc!(d2, p, z2; α=-A2, β=α₂ - 1.0)

                        @. y1 = α₁ - d1
                        @. y2 = α₂ - d2
                        @. y1tmp = 1 - y1

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, α₁, α₂, t1, t2, d1, d2, k)

                        dfunc!(dfy, d, k, y1tmp, s)

                        # Left weights: fw .* dw(z2) .* dfy
                        dwfunc!(dw, p, z2)
                        @. TNL = fw * dw * dfy

                        # Right weights: fw .* dw(A2*z2)
                        dwfunc!(dw, p, z2; α=A2)
                        @. TNR = fw * dw

                        @. A = abs(DIF)^(-2 * s) * DJ

                        I₂ = dot(TNL, A, TNR)

                        IntS[i] = A2 * (1.0 - α₂) *((α₁ + 1.0) * I₁ + (1.0 - α₁) * I₂) / 4.0

                    else
                        error("Fatal error! mapinv not working correctly.")
                    end
                elseif α₁ ≤ -1
                    wfunc!(d1, p, z2; α=-A1, β=α₁ - 1.0)
                    @. y1 = α₁ - d1
                    @. y1tmp = 1 - y1

                    if 1 ≤ α₂
                        # ============================================================
                        # α₁ is below/left, α₂ is above/right
                        # ============================================================

                        wfunc!(d2, p, z1; α=-B2, β=α₂ + 1.0)
                        @. y2 = α₂ - d2

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, α₁, α₂, t1, t2, d1, d2, k)

                        dfunc!(dfy, d, k, y1tmp, s)

                        # Left weights: fw .* dw(A1*z2) .* dfy
                        dwfunc!(dw, p, z2; α=A1)
                        @. TNL = fw * dw * dfy

                        # Right weights: fw .* dw(B2*z1)
                        dwfunc!(dw, p, z1; α=B2)
                        @. TNR = fw * dw

                        @. A = abs(DIF)^(-2 * s) * DJ

                        IntS[i] = A1 * B2 * (1.0 - α₁) * (α₂ + 1.0) * dot(TNL, A, TNR) / 4.0

                    elseif -1 < α₂ && α₂ < 1
                        # ============================================================
                        # α₁ is below/left, α₂ is inside: split α₂ into z1 and z2
                        # ============================================================

                        # ------------ 1st half in α₂ ------------
                        wfunc!(d2, p, z1; α=-1.0, β=α₂ + 1.0)
                        @. y2 = α₂ - d2

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, α₁, α₂, t1, t2, d1, d2, k)

                        dfunc!(dfy, d, k, y1tmp, s)

                        # Left weights: fw .* dw(A1*z2) .* dfy
                        dwfunc!(dw, p, z2; α=A1)
                        @. TNL = fw * dw * dfy

                        # Right weights: fw .* dw(z1)
                        dwfunc!(dw, p, z1)
                        @. TNR = fw * dw

                        @. A = abs(DIF)^(-2 * s) * DJ

                        I₁ = dot(TNL, A, TNR)

                        # ------------ 2nd half in α₂ ------------
                        wfunc!(d2, p, z2; α=-1.0, β=α₂ - 1.0)
                        @. y2 = α₂ - d2

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, α₁, α₂, t1, t2, d1, d2, k)

                        # Left weights are unchanged.
                        # Right weights: fw .* dw(z2)
                        dwfunc!(dw, p, z2)
                        @. TNR = fw * dw

                        @. A = abs(DIF)^(-2 * s) * DJ

                        I₂ = dot(TNL, A, TNR)

                        IntS[i] = A1 * (1.0 - α₁) * ((α₂ + 1.0) * I₁ + (1.0 - α₂) * I₂) / 4.0

                    elseif α₂ ≤ -1
                        # ============================================================
                        # α₁ is below/left, α₂ is below/left
                        # ============================================================

                        wfunc!(d2, p, z2; α=-A2, β=α₂ - 1.0)
                        @. y2 = α₂ - d2

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, α₁, α₂, t1, t2, d1, d2, k)

                        dfunc!(dfy, d, k, y1tmp, s)

                        # Left weights: fw .* dw(A1*z2) .* dfy
                        dwfunc!(dw, p, z2; α=A1)
                        @. TNL = fw * dw * dfy

                        # Right weights: fw .* dw(A2*z2)
                        dwfunc!(dw, p, z2; α=A2)
                        @. TNR = fw * dw

                        @. A = abs(DIF)^(-2 * s) * DJ

                        IntS[i] = A1 * A2 * (1.0 - α₁) * (1.0 - α₂) * dot(TNL, A, TNR) / 4.0
                    end
                end
            else
                # Near singular for boundary patches.
                if 1 ≤ α₁
                    # α₁ = 1
                    wfunc!(d1, p, z1; α=-1.0, β=2.0)
                    @. y1 = 1.0 - d1

                    if 1 ≤ α₂
                        # ============================================================
                        # Boundary patch: α₁ = 1, α₂ above/right
                        # ============================================================

                        wfunc!(d2, p, z1; α=-B2, β=α₂ + 1.0)
                        @. y2 = α₂ - d2

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, 1.0, α₂, t1, t2, d1, d2, k)

                        dfunc!(dfy, d, k, d1, s)   # dfy = dfac(k, 1-y1) = dfac(k, d1)

                        # Left weights: fw .* dw(z1) .* dfy
                        dwfunc!(dw, p, z1; α=1.0)
                        @. TNL = fw * dw * dfy

                        # Right weights: fw .* dw(B2*z1)
                        dwfunc!(dw, p, z1; α=B2)
                        @. TNR = fw * dw

                        @. A = abs(DIF)^(-2 * s) * DJ

                        IntS[i] = B2 * (α₂ + 1.0) * dot(TNL, A, TNR) / 2.0

                    elseif -1 < α₂ && α₂ < 1
                        # ------------ 1st half in α₂ ------------
                        wfunc!(d2, p, z1; α=-1.0, β=α₂ + 1.0)
                        @. y2 = α₂ - d2

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, 1.0, α₂, t1, t2, d1, d2, k)

                        dfunc!(dfy, d, k, d1, s)

                        # Left weights: fw .* dw(z1) .* dfy
                        dwfunc!(dw, p, z1)
                        @. TNL = fw * dw * dfy

                        # Right weights: fw .* dw(z1)
                        # dw already contains dw(z1) from the previous dwfunc! call.
                        @. TNR = fw * dw

                        @. A = abs(DIF)^(-2 * s) * DJ

                        I₁ = dot(TNL, A, TNR)

                        # ------------ 2nd half in α₂ ------------
                        wfunc!(d2, p, z2; α=-1.0, β=α₂ - 1.0)
                        @. y2 = α₂ - d2

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, 1.0, α₂, t1, t2, d1, d2, k)

                        # Left weights are unchanged.
                        # Right weights: fw .* dw(z2)
                        dwfunc!(dw, p, z2)
                        @. TNR = fw * dw

                        @. A = abs(DIF)^(-2 * s) * DJ

                        I₂ = dot(TNL, A, TNR)

                        IntS[i] = ((α₂ + 1.0) * I₁ + (1.0 - α₂) * I₂) / 2.0

                    elseif α₂ ≤ -1
                        wfunc!(d2, p, z2; α=-A2, β=α₂ - 1.0)
                        @. y2 = α₂ - d2

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, 1.0, α₂, t1, t2, d1, d2, k)

                        dfunc!(dfy, d, k, d1, s)

                        # Left weights: fw .* dw(z1) .* dfy
                        dwfunc!(dw, p, z1)
                        @. TNL = fw * dw * dfy

                        # Right weights: fw .* dw(A2*z2)
                        dwfunc!(dw, p, z2; α=A2)
                        @. TNR = fw * dw

                        @. A = abs(DIF)^(-2 * s) * DJ

                        IntS[i] = A2 * (1.0 - α₂) * dot(TNL, A, TNR) / 2.0
                    end
                elseif -1 < α₁ && α₁ < 1
                    if 1 ≤ α₂
                        # ------------ 1st half in α₁ ------------
                        wfunc!(d1, p, z1; α=-1.0, β=α₁ + 1.0)
                        wfunc!(d2, p, z1; α=-B2, β=α₂ + 1.0)

                        @. y1 = α₁ - d1
                        @. y2 = α₂ - d2
                        @. y1tmp = 1 - y1

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, α₁, α₂, t1, t2, d1, d2, k)

                        dfunc!(dfy, d, k, y1tmp, s)

                        # Left weights: fw .* dw(z1) .* dfy
                        dwfunc!(dw, p, z1)
                        @. TNL = fw * dw * dfy

                        # Right weights: fw .* dw(B2*z1)
                        dwfunc!(dw, p, z1; α=B2)
                        @. TNR = fw * dw

                        @. A = abs(DIF)^(-2 * s) * DJ

                        I₁ = dot(TNL, A, TNR)

                        # ------------ 2nd half in α₁ ------------
                        wfunc!(d1, p, z; α=1.0, β=(α₁ - 1.0) / 2.0)
                        wfunc!(d2, p, z1; α=-B2, β=α₂ + 1.0)

                        @. y1 = α₁ - d1
                        @. y2 = α₂ - d2

                        # y1tmp = 1 - y1, using the same changed variable as the original code
                        wfunc!(y1tmp, p, z; α=-1.0, β=(1.0 - α₁) / 2.0)

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, α₁, α₂, t1, t2, d1, d2, k)

                        dfunc!(dfy, d, k, y1tmp, s)

                        # Left weights: fw .* dw(z) .* dfy
                        dwfunc!(dw, p, z)
                        @. TNL = fw * dw * dfy

                        # Right weights: fw .* dw(B2*z1)
                        dwfunc!(dw, p, z1; α=B2)
                        @. TNR = fw * dw

                        @. A = abs(DIF)^(-2 * s) * DJ

                        I₂ = dot(TNL, A, TNR)

                        IntS[i] = B2 * (α₂ + 1.0) * ((α₁ + 1.0) * I₁ + (1.0 - α₁) * I₂) / 4.0

                    elseif α₂ ≤ -1

                        # ------------ 1st half in α₁ ------------
                        wfunc!(d1, p, z1; α=-1.0, β=α₁ + 1.0)
                        wfunc!(d2, p, z2; α=-A2, β=α₂ - 1.0)

                        @. y1 = α₁ - d1
                        @. y2 = α₂ - d2
                        @. y1tmp = 1 - y1

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, α₁, α₂, t1, t2, d1, d2, k)

                        dfunc!(dfy, d, k, y1tmp, s)

                        # Left weights: fw .* dw(z1) .* dfy
                        dwfunc!(dw, p, z1)
                        @. TNL = fw * dw * dfy

                        # Right weights: fw .* dw(A2*z2)
                        dwfunc!(dw, p, z2; α=A2)
                        @. TNR = fw * dw

                        @. A = abs(DIF)^(-2 * s) * DJ

                        I₁ = dot(TNL, A, TNR)

                        # ------------ 2nd half in α₁ ------------
                        wfunc!(d1, p, z; α=1.0, β=(α₁ - 1.0) / 2.0)
                        wfunc!(d2, p, z2; α=-A2, β=α₂ - 1.0)

                        @. y1 = α₁ - d1
                        @. y2 = α₂ - d2

                        # y1tmp = 1 - y1, using the same changed variable as the original code
                        wfunc!(y1tmp, p, z; α=-1.0, β=(1.0 - α₁) / 2.0)

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, α₁, α₂, t1, t2, d1, d2, k)

                        dfunc!(dfy, d, k, y1tmp, s)

                        # Left weights: fw .* dw(z) .* dfy
                        dwfunc!(dw, p, z)
                        @. TNL = fw * dw * dfy

                        # Right weights: fw .* dw(A2*z2)
                        dwfunc!(dw, p, z2; α=A2)
                        @. TNR = fw * dw

                        @. A = abs(DIF)^(-2 * s) * DJ

                        I₂ = dot(TNL, A, TNR)

                        IntS[i] = A2 * (1.0 - α₂) * ((α₁ + 1.0) * I₁ + (1.0 - α₁) * I₂) / 4.0

                    else
                        error("Fatal error! mapinv not working correctly.")
                    end

                elseif α₁ ≤ -1
                    A1 = 1.0 - winv(p, 2.0 * (-1.0 - α₁) / (1.0 - α₁))

                    wfunc!(d1, p, z2; α=-A1, β=(α₁ - 1.0) / 2.0, γ=1.0)
                    @. y1 = α₁ - d1

                    if 1 ≤ α₂

                        wfunc!(d2, p, z1; α=-B2, β=α₂ + 1.0)
                        @. y2 = α₂ - d2

                        wfunc!(y1tmp, p, z2; α=A1, β=(1.0 - α₁) / 2.0, γ=-1.0)

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, α₁, α₂, t1, t2, d1, d2, k)

                        dfunc!(dfy, d, k, y1tmp, s)

                        # Left weights: fw .* dw(special α₁ transform) .* dfy
                        dwfunc!(dw, p, z2; α=-A1, β=1.0, γ=1.0)
                        @. TNL = fw * dw * dfy

                        # Right weights: fw .* dw(B2*z1)
                        dwfunc!(dw, p, z1; α=B2)
                        @. TNR = fw * dw

                        @. A = abs(DIF)^(-2 * s) * DJ

                        IntS[i] = A1 * B2 * (1.0 - α₁) * (α₂ + 1.0) * dot(TNL, A, TNR) / 8.0

                    elseif -1 < α₂ && α₂ < 1

                        # ------------ 1st half in α₂ ------------
                        wfunc!(d2, p, z1; α=-1.0, β=α₂ + 1.0)
                        @. y2 = α₂ - d2

                        wfunc!(y1tmp, p, z2; α=A1, β=(1.0 - α₁) / 2.0, γ=-1.0)

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, α₁, α₂, t1, t2, d1, d2, k)

                        dfunc!(dfy, d, k, y1tmp, s)

                        # Left weights: fw .* dw(special α₁ transform) .* dfy
                        dwfunc!(dw, p, z2; α=-A1, β=1.0, γ=1.0)
                        @. TNL = fw * dw * dfy

                        # Right weights: fw .* dw(z1)
                        dwfunc!(dw, p, z1)
                        @. TNR = fw * dw

                        @. A = abs(DIF)^(-2 * s) * DJ

                        I₁ = dot(TNL, A, TNR)

                        # ------------ 2nd half in α₂ ------------
                        wfunc!(d2, p, z2; α=-1.0, β=α₂ - 1.0)
                        @. y2 = α₂ - d2

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, α₁, α₂, t1, t2, d1, d2, k)

                        # Left weights are unchanged.
                        # Right weights: fw .* dw(z2)
                        dwfunc!(dw, p, z2)
                        @. TNR = fw * dw

                        @. A = abs(DIF)^(-2 * s) * DJ

                        I₂ = dot(TNL, A, TNR)

                        IntS[i] = A1 * (1.0 - α₁) * ((α₂ + 1.0) * I₁ + (1.0 - α₂) * I₂) / 8.0

                    elseif α₂ ≤ -1
                        wfunc!(d2, p, z2; α=-A2, β=α₂ - 1.0)
                        @. y2 = α₂ - d2

                        wfunc!(y1tmp, p, z2; α=A1, β=(1.0 - α₁) / 2.0, γ=-1.0)

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, α₁, α₂, t1, t2, d1, d2, k)

                        dfunc!(dfy, d, k, y1tmp, s)

                        # Left weights: fw .* dw(special α₁ transform) .* dfy
                        dwfunc!(dw, p, z2; α=-A1, β=1.0, γ=1.0)
                        @. TNL = fw * dw * dfy

                        # Right weights: fw .* dw(A2*z2)
                        dwfunc!(dw, p, z2; α=A2)
                        @. TNR = fw * dw

                        @. A = abs(DIF)^(-2 * s) * DJ

                        IntS[i] = A1 * A2 * (1.0 - α₁) * (1.0 - α₂) * dot(TNL, A, TNR) / 8.0
                    end
                end
            end

        end
    end

    return IntS
end