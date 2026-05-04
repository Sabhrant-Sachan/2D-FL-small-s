function precompsLs(d::abstractdomain, dp::domprop, s::Float64, p::Int
    ; n::Int=128)::Matrix{Float64}

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

    # z1 = (1+z)/2, z2 = (1-z)/2 where z = cos(pi*(2j-1)/(2n)) 
    # are Chebyshev nodes in the interval [-1,1]. This is for 
    # near singular integrals 
    #n = n #Maybe n here can be fixed to 128 ? (only for small s)
    z  = Vector{Float64}(undef, n)
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

    wmz₂ = Vector{Float64}(undef, n)
    wz   = Vector{Float64}(undef, n)
    wmz  = Vector{Float64}(undef, n)

    wfunc!(wz, p, z) #wz = w(z)
    wfunc!(wmz, p, z; α=-1.0) #wmz = w(z)
    wfunc!(wmz₂, p, z2; α=-1.0) #wmz₂ = w(-z2)

    fwm = similar(fw)
    fwl = similar(fw)
    dwfunc!(fwl, p, z)
    dwfunc!(fwm, p, z2)   # fwm := dw(z2)
    @. fwm = fw * fwm
    @. fwl = fw * fwl

    y1  = Vector{Float64}(undef, n)
    wdf = similar(y1)
    y1tmp=similar(y1)
    y2  = similar(y1)
    t1  = Matrix{Float64}(undef, n, n)  # meshgrid of y1/y2 
    t2  = Matrix{Float64}(undef, n, n)  # (column = y2[j], row = y1[i])
    DJ  = Matrix{Float64}(undef, n, n)  # To store the Jacobian
    DIF = similar(DJ)
    Zx  = similar(DJ)
    Zy  = similar(DJ)

    TN_y1 = Matrix{Float64}(undef, n, N)
    TN_y2 = Matrix{Float64}(undef, n, N)

    TNL = Matrix{Float64}(undef, N, n)
    TNR = Matrix{Float64}(undef, n, N)
    A   = Matrix{Float64}(undef, n, n)    
    Tmp = Matrix{Float64}(undef, N, n)
    dw  = Vector{Float64}(undef, n)
    dfy = Vector{Float64}(undef, n)

    I = Matrix{Float64}(undef, N, N)

    Lₚ = size(dp.prepts,2) # > dp.pthgo[M+1] - 1

    # First we will go over all the points in the interioir

    Dhc = Vector{Float64}(undef,M)

    @inbounds for k = 1:M
        Dhc[k] = d.pths[k].ck1 - d.pths[k].ck0
    end

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

    # knbd are patches which are not the boundary patches
    # and d.kd are patches are touching the boundary
    knbd = setdiff(collect(1:M), d.kd)

    IntS = zeros(Float64, Np, Lₚ)

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
        x1m, x2m = zp2[rr], zp2[qq]
        x̃₃ = -x1m * x2p / 4.0
        x̃₄ = x1m * x2m / 4.0

        # Singular and k is a boundary patch!
        #-----------------3rd part-----------------

        @. d1 = x1m * wz / 2.0
        @. d2 = x2p * wmz₂
        @. y1 = x1 - d1
        @. y2 = x2 - d2
        @. y1tmp = -x1m * wmz / 2.0

        fill_meshgrid!(t1, t2, y1, y2)

        ChebyTN!(TN_y1, N, y1)
        ChebyTN!(TN_y2, N, y2)

        @inbounds for k in d.kd

            diff_map!(DIF, Zx, Zy, DJ, d, x1, x2, t1, t2, d1, d2, k)

            dfunc!(dfy, d, k, y1tmp, s)   # dfy = dfac(k, 1-y1)

            #Left weights: fwl .* dfy .* TN(y1)'
            @. wdf = fwl * dfy
            @. TNL = TN_y1' * wdf'

            # Right weights: TN(y2).*fwm'
            @. TNR = TN_y2 * fwm   # n×N scaled row-wise
            # Middle terms together
            @. A = expm1(-2 * s * log(DIF)) * DJ

            mul!(Tmp, TNL, A)   # Tmp = TNL * A   (N×n)
            mul!(I, Tmp, TNR)   # I   = Tmp * TNR (N×N)
            
            @inbounds for m in 1:Np
                IntS[m, dp.pthgo[k]+j-1] += x̃₃ * I[m]
            end
        end

        #-----------------4th part-----------------
        #(reuse d1, y1, TN_y1, y1tmp; update d2,y2,t1,t2)
        @. d2 = x2m * wmz₂
        @. y2 = x2 - d2

        fill_meshgrid!(t1, t2, y1, y2)

        ChebyTN!(TN_y2, N, y2)

        @inbounds for k in d.kd

            diff_map!(DIF, Zx, Zy, DJ, d, x1, x2, t1, t2, d1, d2, k)

            dfunc!(dfy, d, k, y1tmp, s)   # dfy = dfac(k, 1-y1)

            #Left weights: fwl .* dfy .* TN(y1)'
            @. wdf = fwl * dfy
            @. TNL = TN_y1' * wdf'

            # Right weights: TN(y2).*fwm'
            @. TNR = TN_y2 * fwm   # n×N scaled row-wise
            # Middle terms together
            @. A = expm1(-2 * s * log(DIF)) * DJ

            mul!(Tmp, TNL, A)   # Tmp = TNL * A   (N×n)
            mul!(I, Tmp, TNR)   # I   = Tmp * TNR (N×N)
            
            @inbounds for m in 1:Np
                IntS[m, dp.pthgo[k]+j-1] += x̃₄ * I[m]
            end
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
        x̃₁ = x2p / 2.0
        x̃₂ = -x2m / 2.0

        #Singular integration is sum of four parts
        #I1,I2,I3,I4, they are all N*N matrices
        #-----------------1st part-----------------

        @. d1 = 2.0 * wmz₂
        @. d2 = x2p * wmz₂
        @. y1 = 1.0 - d1
        @. y2 = x2 - d2

        fill_meshgrid!(t1, t2, y1, y2)

        ChebyTN!(TN_y1, N, y1)
        ChebyTN!(TN_y2, N, y2)

        @inbounds for ℓ in 1:Mbd

            k = d.kd[ℓ]

            diff_map!(DIF, Zx, Zy, DJ, d, 1.0, x2, t1, t2, d1, d2, k)

            dfunc!(dfy, d, k, d1, s)   # dfy = dfac(k, 1-y1=d1)

            #Left weights: fwm .* dfy .* TN(y1)'
            @. wdf = fwm * dfy
            @. TNL = TN_y1' * wdf'

            # Right weights: TN(y2).*fwm'
            @. TNR = TN_y2 * fwm   # n×N scaled row-wise
            # Middle terms together
            @. A = expm1(-2 * s * log(DIF)) * DJ

            mul!(Tmp, TNL, A)   # Tmp = TNL * A   (N×n)
            mul!(I, Tmp, TNR)   # I   = Tmp * TNR (N×N)
            
            @inbounds for m in 1:Np
                IntS[m, Lₚₛ+(ℓ-1)*N+j] += x̃₁ * I[m]
            end
        end

        #-----------------2nd part-----------------
        #(reuse d1, y1, TN_y1; update d2,y2,t1,t2)
        @. d2 = x2m * wmz₂
        @. y2 = x2 - d2

        fill_meshgrid!(t1, t2, y1, y2)

        ChebyTN!(TN_y2, N, y2)

        @inbounds for ℓ in 1:Mbd

            k = d.kd[ℓ]

            diff_map!(DIF, Zx, Zy, DJ, d, 1.0, x2, t1, t2, d1, d2, k)

            dfunc!(dfy, d, k, d1, s)   # dfy = dfac(k, 1-y1)

            #Left weights: fwm .* dfy .* TN(y1)'
            @. wdf = fwm * dfy
            @. TNL = TN_y1' * wdf'

            # Right weights: TN(y2).*fwm'
            @. TNR = TN_y2 * fwm   # n×N scaled row-wise
            # Middle terms together
            @. A = expm1(-2 * s * log(DIF)) * DJ

            mul!(Tmp, TNL, A)   # Tmp = TNL * A   (N×n)
            mul!(I, Tmp, TNR)   # I   = Tmp * TNR (N×N)
            
            @inbounds for m in 1:Np
                IntS[m, Lₚₛ+(ℓ-1)*N+j] += x̃₂ * I[m]
            end
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
        qq, rr = divrem(j - 1, N)
        qq = qq + 1
        rr = rr + 1

        x1, x2 = zp[rr], zp[qq]
        x1p, x2p = zp1[rr], zp1[qq]
        x1m, x2m = zp2[rr], zp2[qq]
        x̃₁ = x1p * x2p / 4.0
        x̃₂ = -x1p * x2m / 4.0
        x̃₃ = -x1m * x2p / 4.0
        x̃₄ = x1m * x2m / 4.0
       
        #Singular integration is sum of four parts
        #I1,I2,I3,I4, they are all N*N matrices
        #-----------------1st part-----------------

        @. d1 = x1p * wmz₂
        @. d2 = x2p * wmz₂
        @. y1 = x1 - d1
        @. y2 = x2 - d2
        @. y1tmp = 1 - y1

        fill_meshgrid!(t1, t2, y1, y2)

        ChebyTN!(TN_y1, N, y1)
        ChebyTN!(TN_y2, N, y2)

        @inbounds for k in 1:M

            diff_map!(DIF, Zx, Zy, DJ, d, x1, x2, t1, t2, d1, d2, k)

            dfunc!(dfy, d, k, y1tmp, s)   # dfy = dfac(k, 1-y1)

            #Left weights: fwm .* dfy .* TN(y1)'
            @. wdf = fwm * dfy
            @. TNL = TN_y1' * wdf'

            # Right weights: TN(y2).*fwm'
            @. TNR = TN_y2 * fwm   # n×N scaled row-wise
            # Middle terms together
            @. A = expm1(-2 * s * log(DIF)) * DJ

            mul!(Tmp, TNL, A)   # Tmp = TNL * A   (N×n)
            mul!(I, Tmp, TNR)   # I   = Tmp * TNR (N×N)
            
            @inbounds for m in 1:Np
                IntS[m, dp.pthgo[k]+j-1] += x̃₁ * I[m]
            end
        end

        #-----------------2nd part-----------------
        #(reuse d1, y1, TN_y1, y1tmp; update d2,y2,t1,t2)
        @. d2 = x2m * wmz₂
        @. y2 = x2 - d2

        fill_meshgrid!(t1, t2, y1, y2)

        ChebyTN!(TN_y2, N, y2)

        @inbounds for k in 1:M

            diff_map!(DIF, Zx, Zy, DJ, d, x1, x2, t1, t2, d1, d2, k)

            dfunc!(dfy, d, k, y1tmp, s)   # dfy = dfac(k, 1-y1)

            #Left weights: fwm .* dfy .* TN(y1)'
            @. wdf = fwm * dfy
            @. TNL = TN_y1' * wdf'

            # Right weights: TN(y2).*fwm'
            @. TNR = TN_y2 * fwm   # n×N scaled row-wise
            # Middle terms together
            @. A = expm1(-2 * s * log(DIF)) * DJ

            mul!(Tmp, TNL, A)   # Tmp = TNL * A   (N×n)
            mul!(I, Tmp, TNR)   # I   = Tmp * TNR (N×N)
            
            @inbounds for m in 1:Np
                IntS[m, dp.pthgo[k]+j-1] += x̃₂ * I[m]
            end
        end

        #-----------------3rd part-----------------

        @. d1 = x1m * wmz₂
        @. d2 = x2p * wmz₂
        @. y1 = x1 - d1
        @. y2 = x2 - d2
        @. y1tmp = 1 - y1

        fill_meshgrid!(t1, t2, y1, y2)

        ChebyTN!(TN_y1, N, y1)
        ChebyTN!(TN_y2, N, y2)

        @inbounds for k in knbd

            diff_map!(DIF, Zx, Zy, DJ, d, x1, x2, t1, t2, d1, d2, k)

            dfunc!(dfy, d, k, y1tmp, s)   # dfy = dfac(k, 1-y1)

            #Left weights: fwm .* dfy .* TN(y1)'
            @. wdf = fwm * dfy
            @. TNL = TN_y1' * wdf'

            # Right weights: TN(y2).*fwm'
            @. TNR = TN_y2 * fwm   # n×N scaled row-wise
            # Middle terms together
            @. A = expm1(-2 * s * log(DIF)) * DJ

            mul!(Tmp, TNL, A)   # Tmp = TNL * A   (N×n)
            mul!(I, Tmp, TNR)   # I   = Tmp * TNR (N×N)
            
            @inbounds for m in 1:Np
                IntS[m, dp.pthgo[k]+j-1] += x̃₃ * I[m]
            end
        end

        #-----------------4th part-----------------
        #(reuse d1, y1, TN_y1, y1tmp; update d2,y2,t1,t2)
        @. d2 = x2m * wmz₂
        @. y2 = x2 - d2

        fill_meshgrid!(t1, t2, y1, y2)

        ChebyTN!(TN_y2, N, y2)

        @inbounds for k in knbd

            diff_map!(DIF, Zx, Zy, DJ, d, x1, x2, t1, t2, d1, d2, k)

            dfunc!(dfy, d, k, y1tmp, s)   # dfy = dfac(k, 1-y1)

            #Left weights: fwm .* dfy .* TN(y1)'
            @. wdf = fwm * dfy
            @. TNL = TN_y1' * wdf'

            # Right weights: TN(y2).*fwm'
            @. TNR = TN_y2 * fwm   # n×N scaled row-wise
            # Middle terms together
            @. A = expm1(-2 * s * log(DIF)) * DJ

            mul!(Tmp, TNL, A)   # Tmp = TNL * A   (N×n)
            mul!(I, Tmp, TNR)   # I   = Tmp * TNR (N×N)
            
            @inbounds for m in 1:Np
                IntS[m, dp.pthgo[k]+j-1] += x̃₄ * I[m]
            end
        end

    end

    #--------------------------------------------
    #----------Near Singular Integration---------
    #-------------------------------------------
    I₁ = Matrix{Float64}(undef, N, N)
    I₂ = Matrix{Float64}(undef, N, N)

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

                        ChebyTN!(TN_y1, N, y1)
                        ChebyTN!(TN_y2, N, y2)

                        #Left weights: fw .* dw(B1*z1) .* dfy .* TN(y1)'
                        dwfunc!(dw, p, z1; α=B1)
                        @. wdf = fw * dw * dfy
                        @. TNL = TN_y1' * wdf'

                        # Right weights: TN(y2).*dw(B2*z1)'.*fw'
                        dwfunc!(dw, p, z1; α=B2)
                        @. TNR = TN_y2 * dw * fw   # n×N scaled row-wise
                        # Middle terms together
                        @. A = expm1(-2*s*log(DIF)) * DJ   

                        mul!(Tmp, TNL, A)   # Tmp = TNL * A   (N×n)
                        mul!(I, Tmp, TNR)   # I   = Tmp * TNR (N×N)

                        I .*= B1 * B2 * (α₁ + 1) * (α₂ + 1) / 4

                    elseif -1 < α₂ && α₂ < 1
                        wfunc!(d2, p, z1; α=-1.0, β=α₂ + 1.0)
                        @. y2 = α₂ - d2

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, α₁, α₂, t1, t2, d1, d2, k)

                        dfunc!(dfy, d, k, y1tmp, s)   # dfy = dfac(k, 1-y1)

                        ChebyTN!(TN_y1, N, y1)
                        ChebyTN!(TN_y2, N, y2)

                        # Left weights: fw .* dw(B1*z1) .* dfy .* TN(y1)'
                        dwfunc!(dw, p, z1; α=B1)
                        @. wdf = fw * dw * dfy
                        @. TNL = TN_y1' * wdf'

                        # Right weights: TN(y2).*dw(z1)'.*fw'
                        dwfunc!(dw, p, z1)
                        @. TNR = TN_y2 * dw * fw   # n×N scaled row-wise
                        # Middle terms together
                        @. A = expm1(-2*s*log(DIF)) * DJ   

                        mul!(Tmp, TNL, A)   # Tmp = TNL * A   (N×n)
                        mul!(I₁, Tmp, TNR)   # I   = Tmp * TNR (N×N)

                        #------------2ⁿᵈ half ------------
                        wfunc!(d2, p, z2; α=-1.0, β=α₂ - 1.0)
                        @. y2 = α₂ - d2

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, α₁, α₂, t1, t2, d1, d2, k)

                        #Left weights same as before
                        #Right weights: TN(y2).*dw(z2)'.*fw'
                        ChebyTN!(TN_y2, N, y2)
                        dwfunc!(dw, p, z2)
                        @. TNR = TN_y2 * dw * fw   # n×N scaled row-wise
                        # Middle terms together
                        @. A = expm1(-2*s*log(DIF)) * DJ   

                        mul!(Tmp, TNL, A)   # Tmp = TNL * A   (N×n)
                        mul!(I₂, Tmp, TNR)  # I   = Tmp * TNR (N×N)

                        @. I = B1*(α₁ + 1)*( (α₂ + 1)*I₁ + (1 - α₂)*I₂ )/4

                    elseif α₂ ≤ -1
                        wfunc!(d2, p, z2; α=-A2, β=α₂ - 1.0)
                        @. y2 = α₂ - d2

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, α₁, α₂, t1, t2, d1, d2, k)

                        dfunc!(dfy, d, k, y1tmp, s)   # dfy = dfac(k, 1-y1)

                        ChebyTN!(TN_y1, N, y1)
                        ChebyTN!(TN_y2, N, y2)

                        # Left weights: fw .* dw(B1*z1) .* dfy .* TN(y1)'
                        dwfunc!(dw, p, z1; α=B1)
                        @. wdf = fw * dw * dfy
                        @. TNL = TN_y1' * wdf'

                        # Right weights: TN(y2).*dw(A2*z2)'.*fw'
                        dwfunc!(dw, p, z2; α=A2)
                        @. TNR = TN_y2 * dw * fw   # n×N scaled row-wise
                        # Middle terms together
                        @. A = expm1(-2*s*log(DIF)) * DJ   

                        mul!(Tmp, TNL, A)   # Tmp = TNL * A   (N×n)
                        mul!(I, Tmp, TNR)   # I   = Tmp * TNR (N×N)

                        I .*= B1 * A2 * (α₁ + 1) * (1 - α₂) / 4
                    end
                elseif -1 < α₁ && α₁ < 1
                    if 1 ≤ α₂
                        wfunc!(d1, p, z1; α=-1.0, β=α₁ + 1.0)
                        wfunc!(d2, p, z1; α=-B2, β=α₂ + 1.0)
                        @. y1 = α₁ - d1
                        @. y2 = α₂ - d2 
                        @. y1tmp = 1 - y1

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, α₁, α₂, t1, t2, d1, d2, k)

                        dfunc!(dfy, d, k, y1tmp, s)   # dfy = dfac(k, 1-y1)

                        ChebyTN!(TN_y1, N, y1)
                        ChebyTN!(TN_y2, N, y2)

                        # Left weights: fw .* dw(z1) .* dfy .* TN(y1)'
                        dwfunc!(dw, p, z1)
                        @. wdf = fw * dw * dfy
                        @. TNL = TN_y1' * wdf'

                        # Right weights: TN(y2).*dw(B2*z1)'.*fw'
                        dwfunc!(dw, p, z1; α=B2)
                        @. TNR = TN_y2 * dw * fw   # n×N scaled row-wise
                        # Middle terms together
                        @. A = expm1(-2*s*log(DIF)) * DJ   

                        mul!(Tmp, TNL, A)   # Tmp = TNL * A   (N×n)
                        mul!(I₁, Tmp, TNR)  # I   = Tmp * TNR (N×N)

                        #------------2ⁿᵈ half ------------
                        wfunc!(d1, p, z2; α=-1.0, β=α₁ - 1.0)
                        wfunc!(d2, p, z1; α=-B2, β=α₂ + 1.0)
                        @. y1 = α₁ - d1
                        @. y2 = α₂ - d2 
                        @. y1tmp = 1 - y1

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, α₁, α₂, t1, t2, d1, d2, k)

                        dfunc!(dfy, d, k, y1tmp, s)   # dfy = dfac(k, 1-y1)

                        ChebyTN!(TN_y1, N, y1)
                        ChebyTN!(TN_y2, N, y2)

                        # Left weights: fw .* dw(z2) .* dfy .* TN(y1)'
                        dwfunc!(dw, p, z2)
                        @. wdf = fw * dw * dfy
                        @. TNL = TN_y1' * wdf'

                        # Right weights: TN(y2).*dw(B2*z1)'.*fw'
                        dwfunc!(dw, p, z1; α=B2)
                        @. TNR = TN_y2 * dw * fw   # n×N scaled row-wise
                        # Middle terms together
                        @. A = expm1(-2*s*log(DIF)) * DJ   

                        mul!(Tmp, TNL, A)   # Tmp = TNL * A   (N×n)
                        mul!(I₂, Tmp, TNR)  # I   = Tmp * TNR (N×N)

                        @. I = B2*(α₂ + 1)*( (α₁ + 1)*I₁ + (1 - α₁)*I₂ )/4

                    elseif α₂ ≤ -1
                        wfunc!(d1, p, z1; α=-1.0, β=α₁ + 1.0)
                        wfunc!(d2, p, z2; α=-A2, β=α₂ - 1.0)
                        @. y1 = α₁ - d1
                        @. y2 = α₂ - d2 
                        @. y1tmp = 1 - y1

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, α₁, α₂, t1, t2, d1, d2, k)

                        dfunc!(dfy, d, k, y1tmp, s)   # dfy = dfac(k, 1-y1)

                        ChebyTN!(TN_y1, N, y1)
                        ChebyTN!(TN_y2, N, y2)

                        # Left weights: fw .* dw(z1) .* dfy .* TN(y1)'
                        dwfunc!(dw, p, z1)
                        @. wdf = fw * dw * dfy
                        @. TNL = TN_y1' * wdf'

                        # Right weights: TN(y2).*dw(A2*z2)'.*fw'
                        dwfunc!(dw, p, z2; α=A2)
                        @. TNR = TN_y2 * dw * fw   # n×N scaled row-wise
                        # Middle terms together
                        @. A = expm1(-2*s*log(DIF)) * DJ   

                        mul!(Tmp, TNL, A)   # Tmp = TNL * A   (N×n)
                        mul!(I₁, Tmp, TNR)  # I   = Tmp * TNR (N×N)

                        #------------2ⁿᵈ half ------------
                        wfunc!(d1, p, z2; α=-1.0, β=α₁ - 1.0)
                        wfunc!(d2, p, z2; α=-A2, β=α₂ - 1.0)
                        @. y1 = α₁ - d1
                        @. y2 = α₂ - d2 
                        @. y1tmp = 1 - y1

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, α₁, α₂, t1, t2, d1, d2, k)

                        dfunc!(dfy, d, k, y1tmp, s)   # dfy = dfac(k, 1-y1)

                        ChebyTN!(TN_y1, N, y1)
                        ChebyTN!(TN_y2, N, y2)

                        # Left weights: fw .* dw(z2) .* dfy .* TN(y1)'
                        dwfunc!(dw, p, z2)
                        @. wdf = fw * dw * dfy
                        @. TNL = TN_y1' * wdf'

                        # Right weights: TN(y2).*dw(A2*z2)'.*fw'
                        dwfunc!(dw, p, z2; α=A2)
                        @. TNR = TN_y2 * dw * fw   # n×N scaled row-wise
                        # Middle terms together
                        @. A = expm1(-2*s*log(DIF)) * DJ   

                        mul!(Tmp, TNL, A)   # Tmp = TNL * A   (N×n)
                        mul!(I₂, Tmp, TNR)  # I   = Tmp * TNR (N×N)

                        @. I = A2*(1 - α₂)*( (α₁ + 1)*I₁ + (1 - α₁)*I₂ )/4
                    else
                        error("Fatal error! mapinv not working correctly.")
                    end

                elseif α₁ ≤ -1
                    wfunc!(d1, p, z2; α=-A1, β=α₁ - 1.0)
                    @. y1 = α₁ - d1
                    @. y1tmp = 1 - y1
                    if 1 ≤ α₂
                        wfunc!(d2, p, z1; α=-B2, β=α₂ + 1.0)
                        @. y2 = α₂ - d2

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, α₁, α₂, t1, t2, d1, d2, k)

                        dfunc!(dfy, d, k, y1tmp, s)   # dfy = dfac(k, 1-y1)

                        ChebyTN!(TN_y1, N, y1)
                        ChebyTN!(TN_y2, N, y2)

                        #Left weights: fw .* dw(A1*z2) .* dfy .* TN(y1)'
                        dwfunc!(dw, p, z2; α=A1)
                        @. wdf = fw * dw * dfy
                        @. TNL = TN_y1' * wdf'

                        # Right weights: TN(y2).*dw(B2*z1)'.*fw'
                        dwfunc!(dw, p, z1; α=B2)
                        @. TNR = TN_y2 * dw * fw   # n×N scaled row-wise
                        # Middle terms together
                        @. A = expm1(-2*s*log(DIF)) * DJ   

                        mul!(Tmp, TNL, A)   # Tmp = TNL * A   (N×n)
                        mul!(I, Tmp, TNR)   # I   = Tmp * TNR (N×N)

                        I .*= A1 * B2 * (1 - α₁) * (α₂ + 1) / 4

                    elseif -1 < α₂ && α₂ < 1
                        wfunc!(d2, p, z1; α=-1.0, β=α₂ + 1.0)
                        @. y2 = α₂ - d2

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, α₁, α₂, t1, t2, d1, d2, k)

                        dfunc!(dfy, d, k, y1tmp, s)   # dfy = dfac(k, 1-y1)

                        ChebyTN!(TN_y1, N, y1)
                        ChebyTN!(TN_y2, N, y2)

                        # Left weights: fw .* dw(A1*z2) .* dfy .* TN(y1)'
                        dwfunc!(dw, p, z2; α=A1)
                        @. wdf = fw * dw * dfy
                        @. TNL = TN_y1' * wdf'

                        # Right weights: TN(y2).*dw(z1)'.*fw'
                        dwfunc!(dw, p, z1)
                        @. TNR = TN_y2 * dw * fw   # n×N scaled row-wise
                        # Middle terms together
                        @. A = expm1(-2*s*log(DIF)) * DJ   

                        mul!(Tmp, TNL, A)   # Tmp = TNL * A   (N×n)
                        mul!(I₁, Tmp, TNR)   # I   = Tmp * TNR (N×N)

                        #------------2ⁿᵈ half ------------
                        wfunc!(d2, p, z2; α=-1.0, β=α₂ - 1.0)
                        @. y2 = α₂ - d2

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, α₁, α₂, t1, t2, d1, d2, k)

                        #Left weights same as before
                        #Right weights: TN(y2).*dw(z2)'.*fw'
                        ChebyTN!(TN_y2, N, y2)
                        dwfunc!(dw, p, z2)
                        @. TNR = TN_y2 * dw * fw   # n×N scaled row-wise
                        # Middle terms together
                        @. A = expm1(-2*s*log(DIF)) * DJ   

                        mul!(Tmp, TNL, A)   # Tmp = TNL * A   (N×n)
                        mul!(I₂, Tmp, TNR)  # I   = Tmp * TNR (N×N)

                        @. I = A1*(1 - α₁)*( (α₂ + 1)*I₁ + (1 - α₂)*I₂ )/4
                    elseif α₂ ≤ -1
                        wfunc!(d2, p, z2; α=-A2, β=α₂ - 1.0)
                        @. y2 = α₂ - d2;

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, α₁, α₂, t1, t2, d1, d2, k)

                        dfunc!(dfy, d, k, y1tmp, s)   # dfy = dfac(k, 1-y1)

                        ChebyTN!(TN_y1, N, y1)
                        ChebyTN!(TN_y2, N, y2)

                        # Left weights: fw .* dw(A1*z2) .* dfy .* TN(y1)'
                        dwfunc!(dw, p, z2; α=A1)
                        @. wdf = fw * dw * dfy
                        @. TNL = TN_y1' * wdf'

                        # Right weights: TN(y2).*dw(A2*z2)'.*fw'
                        dwfunc!(dw, p, z2; α=A2)
                        @. TNR = TN_y2 * dw * fw   # n×N scaled row-wise
                        # Middle terms together
                        @. A = expm1(-2*s*log(DIF)) * DJ   

                        mul!(Tmp, TNL, A)   # Tmp = TNL * A   (N×n)
                        mul!(I, Tmp, TNR)   # I   = Tmp * TNR (N×N)

                        I .*= A1 * A2 * (1 - α₁) * (1 - α₂) / 4
                    end
                end
            else
                #Near singular for boundary patches.
                if 1 ≤ α₁ 
                    #α₁ = 1
                    wfunc!(d1, p, z1; α=-1.0, β=2.0)
                    @. y1 = 1 - d1
                    if 1 ≤ α₂
                        wfunc!(d2, p, z1; α=-B2, β=α₂ + 1.0)
                        @. y2 = α₂ - d2

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, 1.0, α₂, t1, t2, d1, d2, k)

                        dfunc!(dfy, d, k, d1, s)   # dfy = dfac(k, 1-y1) = dfac(k, d1)

                        ChebyTN!(TN_y1, N, y1)
                        ChebyTN!(TN_y2, N, y2)

                        #Left weights: fw .* dw(z1) .* dfy .* TN(y1)'
                        dwfunc!(dw, p, z1; α=1.0)
                        @. wdf = fw * dw * dfy
                        @. TNL = TN_y1' * wdf'

                        # Right weights: TN(y2).*dw(B2*z1)'.*fw'
                        dwfunc!(dw, p, z1; α=B2)
                        @. TNR = TN_y2 * dw * fw   # n×N scaled row-wise
                        # Middle terms together
                        @. A = expm1(-2*s*log(DIF)) * DJ   

                        mul!(Tmp, TNL, A)   # Tmp = TNL * A   (N×n)
                        mul!(I, Tmp, TNR)   # I   = Tmp * TNR (N×N)

                        I .*= B2 * (α₂ + 1) / 2

                    elseif -1 < α₂ && α₂ < 1
                        wfunc!(d2, p, z1; α=-1.0, β=α₂ + 1.0)
                        @. y2 = α₂ - d2

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, 1.0, α₂, t1, t2, d1, d2, k)

                        dfunc!(dfy, d, k, d1, s)   # dfy = dfac(k, 1-y1)

                        ChebyTN!(TN_y1, N, y1)
                        ChebyTN!(TN_y2, N, y2)

                        # Left weights: fw .* dw(z1) .* dfy .* TN(y1)'
                        dwfunc!(dw, p, z1)
                        @. wdf = fw * dw * dfy
                        @. TNL = TN_y1' * wdf'

                        # Right weights: TN(y2).*dw(z1)'.*fw'
                        @. TNR = TN_y2 * dw * fw   # n×N scaled row-wise
                        # Middle terms together
                        @. A = expm1(-2*s*log(DIF)) * DJ

                        mul!(Tmp, TNL, A)   # Tmp = TNL * A   (N×n)
                        mul!(I₁, Tmp, TNR)   # I   = Tmp * TNR (N×N)

                        #------------2ⁿᵈ half ------------
                        wfunc!(d2, p, z2; α=-1.0, β=α₂ - 1.0)
                        @. y2 = α₂ - d2

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, 1.0, α₂, t1, t2, d1, d2, k)

                        #Left weights same as before
                        #Right weights: TN(y2).*dw(z2)'.*fw'
                        ChebyTN!(TN_y2, N, y2)
                        dwfunc!(dw, p, z2)
                        @. TNR = TN_y2 * dw * fw   # n×N scaled row-wise
                        # Middle terms together
                        @. A = expm1(-2*s*log(DIF)) * DJ

                        mul!(Tmp, TNL, A)   # Tmp = TNL * A   (N×n)
                        mul!(I₂, Tmp, TNR)  # I   = Tmp * TNR (N×N)

                        @. I = ((α₂ + 1) * I₁ + (1 - α₂) * I₂) / 2

                    elseif α₂ ≤ -1
                        wfunc!(d2, p, z2; α=-A2, β=α₂ - 1.0)
                        @. y2 = α₂ - d2

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, 1.0, α₂, t1, t2, d1, d2, k)

                        dfunc!(dfy, d, k, d1, s)   # dfy = dfac(k, 1-y1)

                        ChebyTN!(TN_y1, N, y1)
                        ChebyTN!(TN_y2, N, y2)

                        # Left weights: fw .* dw(z1) .* dfy .* TN(y1)'
                        dwfunc!(dw, p, z1)
                        @. wdf = fw * dw * dfy
                        @. TNL = TN_y1' * wdf'

                        # Right weights: TN(y2).*dw(A2*z2)'.*fw'
                        dwfunc!(dw, p, z2; α=A2)
                        @. TNR = TN_y2 * dw * fw   # n×N scaled row-wise
                        # Middle terms together
                        @. A = expm1(-2*s*log(DIF)) * DJ   

                        mul!(Tmp, TNL, A)   # Tmp = TNL * A   (N×n)
                        mul!(I, Tmp, TNR)   # I   = Tmp * TNR (N×N)

                        I .*=  A2 * (1 - α₂) / 2
                    end
                elseif -1 < α₁ && α₁ < 1
                    if 1 ≤ α₂
                        wfunc!(d1, p, z1; α=-1.0, β=α₁ + 1.0)
                        wfunc!(d2, p, z1; α=-B2, β=α₂ + 1.0)
                        @. y1 = α₁ - d1
                        @. y2 = α₂ - d2 
                        @. y1tmp = 1 - y1

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, α₁, α₂, t1, t2, d1, d2, k)

                        dfunc!(dfy, d, k, y1tmp, s)   # dfy = dfac(k, 1-y1)

                        ChebyTN!(TN_y1, N, y1)
                        ChebyTN!(TN_y2, N, y2)

                        # Left weights: fw .* dw .* dfy .* TN(y1)'
                        dwfunc!(dw, p, z1)
                        @. wdf = fw * dw * dfy
                        @. TNL = TN_y1' * wdf'

                        # Right weights: TN(y2).*dw(B2*z1)'.*fw'
                        dwfunc!(dw, p, z1; α=B2)
                        @. TNR = TN_y2 * dw * fw   # n×N scaled row-wise
                        # Middle terms together
                        @. A = expm1(-2*s*log(DIF)) * DJ   

                        mul!(Tmp, TNL, A)   # Tmp = TNL * A   (N×n)
                        mul!(I₁, Tmp, TNR)  # I   = Tmp * TNR (N×N)

                        #------------2ⁿᵈ half ------------
                        wfunc!(d1, p, z; α=1.0, β=(α₁ - 1.0) / 2)
                        wfunc!(d2, p, z1; α=-B2, β=α₂ + 1.0)
                        @. y1 = α₁ - d1
                        @. y2 = α₂ - d2 
                        wfunc!(y1tmp, p, z; α=-1.0, β=(1.0 - α₁) / 2)

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, α₁, α₂, t1, t2, d1, d2, k)

                        dfunc!(dfy, d, k, y1tmp, s)   # dfy = dfac(k, 1-y1)

                        ChebyTN!(TN_y1, N, y1)
                        ChebyTN!(TN_y2, N, y2)

                        # Left weights: fw .* dw .* dfy .* TN(y1)'
                        dwfunc!(dw, p, z)
                        @. wdf = fw * dw * dfy
                        @. TNL = TN_y1' * wdf'

                        # Right weights: TN(y2).*dw(B2*z1)'.*fw'
                        dwfunc!(dw, p, z1; α=B2)
                        @. TNR = TN_y2 * dw * fw   # n×N scaled row-wise
                        # Middle terms together
                        @. A = expm1(-2*s*log(DIF)) * DJ   

                        mul!(Tmp, TNL, A)   # Tmp = TNL * A   (N×n)
                        mul!(I₂, Tmp, TNR)  # I   = Tmp * TNR (N×N)

                        @. I = B2*(α₂ + 1)*( (α₁ + 1)*I₁ + (1 - α₁) * I₂ )/4

                    elseif α₂ ≤ -1
                        wfunc!(d1, p, z1; α=-1.0, β=α₁ + 1.0)
                        wfunc!(d2, p, z2; α=-A2, β=α₂ - 1.0)
                        @. y1 = α₁ - d1
                        @. y2 = α₂ - d2 
                        @. y1tmp = 1 - y1

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, α₁, α₂, t1, t2, d1, d2, k)

                        dfunc!(dfy, d, k, y1tmp, s)   # dfy = dfac(k, 1-y1)

                        ChebyTN!(TN_y1, N, y1)
                        ChebyTN!(TN_y2, N, y2)

                        # Left weights: fw .* dw(z1) .* dfy .* TN(y1)'
                        dwfunc!(dw, p, z1)
                        @. wdf = fw * dw * dfy
                        @. TNL = TN_y1' * wdf'

                        # Right weights: TN(y2).*dw(A2*z2)'.*fw'
                        dwfunc!(dw, p, z2; α=A2)
                        @. TNR = TN_y2 * dw * fw   # n×N scaled row-wise
                        # Middle terms together
                        @. A = expm1(-2*s*log(DIF)) * DJ   

                        mul!(Tmp, TNL, A)   # Tmp = TNL * A   (N×n)
                        mul!(I₁, Tmp, TNR)  # I   = Tmp * TNR (N×N)

                        #------------2ⁿᵈ half ------------
                        wfunc!(d1, p, z; α=1.0, β=(α₁ - 1.0) / 2)
                        wfunc!(d2, p, z2; α=-A2, β=α₂ - 1.0)
                        @. y1 = α₁ - d1
                        @. y2 = α₂ - d2 
                        wfunc!(y1tmp, p, z; α=-1.0, β=(1.0 - α₁) / 2)

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, α₁, α₂, t1, t2, d1, d2, k)

                        dfunc!(dfy, d, k, y1tmp, s)   # dfy = dfac(k, 1-y1)

                        ChebyTN!(TN_y1, N, y1)
                        ChebyTN!(TN_y2, N, y2)

                        # Left weights: fw .* dw .* dfy .* TN(y1)'
                        dwfunc!(dw, p, z)
                        @. wdf = fw * dw * dfy
                        @. TNL = TN_y1' * wdf'

                        # Right weights: TN(y2).*dw(A2*z2)'.*fw'
                        dwfunc!(dw, p, z2; α=A2)
                        @. TNR = TN_y2 * dw * fw   # n×N scaled row-wise
                        # Middle terms together
                        @. A = expm1(-2*s*log(DIF)) * DJ   

                        mul!(Tmp, TNL, A)   # Tmp = TNL * A   (N×n)
                        mul!(I₂, Tmp, TNR)  # I   = Tmp * TNR (N×N)

                        @. I = A2*(1 - α₂)*( (α₁ + 1) * I₁ + (1 - α₁) * I₂ )/4
                    else
                        error("Fatal error! mapinv not working correctly.")
                    end

                elseif α₁ ≤ -1
                    A1 = 1 - winv(p, 2 * (-1 - α₁) / (1 - α₁))
                    wfunc!(d1, p, z2; α=-A1, β=(α₁ - 1.0) / 2, γ=1.0)
                    @. y1 = α₁ - d1
                    if 1 ≤ α₂
                        wfunc!(d2, p, z1; α=-B2, β=α₂ + 1.0)
                        @. y2 = α₂ - d2

                        wfunc!(y1tmp, p, z2; α=A1, β=(1.0 - α₁) / 2, γ=-1.0)

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, α₁, α₂, t1, t2, d1, d2, k)

                        dfunc!(dfy, d, k, y1tmp, s)   # dfy = dfac(k, 1-y1)

                        ChebyTN!(TN_y1, N, y1)
                        ChebyTN!(TN_y2, N, y2)

                        #Left weights: fw .* dw .* dfy .* TN(y1)'
                        dwfunc!(dw, p, z2; α=-A1, β=1.0, γ=1.0)
                        @. wdf = fw * dw * dfy
                        @. TNL = TN_y1' * wdf'

                        # Right weights: TN(y2).*dw(B2*z1)'.*fw'
                        dwfunc!(dw, p, z1; α=B2)
                        @. TNR = TN_y2 * dw * fw   # n×N scaled row-wise
                        # Middle terms together
                        @. A = expm1(-2*s*log(DIF)) * DJ   

                        mul!(Tmp, TNL, A)   # Tmp = TNL * A   (N×n)
                        mul!(I, Tmp, TNR)   # I   = Tmp * TNR (N×N)

                        I .*=  A1 * B2 * (1-α₁)*(α₂+1)/8

                    elseif -1 < α₂ && α₂ < 1
                        wfunc!(d2, p, z1; α=-1.0, β=α₂ + 1.0)
                        @. y2 = α₂ - d2

                        wfunc!(y1tmp, p, z2; α=A1, β=(1.0 - α₁) / 2, γ=-1.0)

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, α₁, α₂, t1, t2, d1, d2, k)

                        dfunc!(dfy, d, k, y1tmp, s)   # dfy = dfac(k, 1-y1)

                        ChebyTN!(TN_y1, N, y1)
                        ChebyTN!(TN_y2, N, y2)

                        #Left weights: fw .* dw .* dfy .* TN(y1)'
                        dwfunc!(dw, p, z2; α=-A1, β=1.0, γ=1.0)
                        @. wdf = fw * dw * dfy
                        @. TNL = TN_y1' * wdf'

                        # Right weights: TN(y2).*dw(z1)'.*fw'
                        dwfunc!(dw, p, z1)
                        @. TNR = TN_y2 * dw * fw   # n×N scaled row-wise
                        # Middle terms together
                        @. A = expm1(-2*s*log(DIF)) * DJ   

                        mul!(Tmp, TNL, A)   # Tmp = TNL * A   (N×n)
                        mul!(I₁, Tmp, TNR)   # I   = Tmp * TNR (N×N)

                        #------------2ⁿᵈ half ------------
                        wfunc!(d2, p, z2; α=-1.0, β=α₂ - 1.0)
                        @. y2 = α₂ - d2

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, α₁, α₂, t1, t2, d1, d2, k)

                        #Left weights same as before
                        #Right weights: TN(y2).*dw(z2)'.*fw'
                        ChebyTN!(TN_y2, N, y2)
                        dwfunc!(dw, p, z2)
                        @. TNR = TN_y2 * dw * fw   # n×N scaled row-wise
                        # Middle terms together
                        @. A = expm1(-2*s*log(DIF)) * DJ   

                        mul!(Tmp, TNL, A)   # Tmp = TNL * A   (N×n)
                        mul!(I₂, Tmp, TNR)  # I   = Tmp * TNR (N×N)

                        @. I = A1*(1-α₁)*( (α₂+1)*I₁ + (1-α₂)*I₂ )/8
                        
                    elseif α₂ ≤ -1
                        wfunc!(d2, p, z2; α=-A2, β=α₂ - 1.0)
                        @. y2 = α₂ - d2

                        wfunc!(y1tmp, p, z2; α=A1, β=(1.0 - α₁) / 2, γ=-1.0)

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, α₁, α₂, t1, t2, d1, d2, k)

                        dfunc!(dfy, d, k, y1tmp, s)   # dfy = dfac(k, 1-y1)

                        ChebyTN!(TN_y1, N, y1)
                        ChebyTN!(TN_y2, N, y2)

                        #Left weights: fw .* dw .* dfy .* TN(y1)'
                        dwfunc!(dw, p, z2; α=-A1, β=1.0, γ=1.0)
                        @. wdf = fw * dw * dfy
                        @. TNL = TN_y1' * wdf'

                        # Right weights: TN(y2).*dw(A2*z2)'.*fw'
                        dwfunc!(dw, p, z2; α=A2)
                        @. TNR = TN_y2 * dw * fw   # n×N scaled row-wise
                        # Middle terms together
                        @. A = expm1(-2*s*log(DIF)) * DJ   
                        
                        mul!(Tmp, TNL, A)   # Tmp = TNL * A   (N×n)
                        mul!(I, Tmp, TNR)   # I   = Tmp * TNR (N×N)

                        I .*= A1*A2*(1-α₁)*(1-α₂)/8

                    end
                end
            end

            #Allocation free version of IntS[:,i] = reshape(I,(Np,1));
            @inbounds for m in 1:Np
                IntS[m, i] = I[m]
            end

        end
    end
  
    return IntS
end

