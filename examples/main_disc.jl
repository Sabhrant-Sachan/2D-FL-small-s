using Revise, Dates, FL2D_small

using FL2D_small.FLdata

dobenchmark, docondnum = false, true

# =============================================
#Direct solver for many value of s

open("solve_outputs0500_direct.txt", "w") do io
    # =============================================
    s, p = 0.5, 4

    println(io, "Run started: ", Dates.now())
    println(io, "s = ", s)
    
    δ, δclsbd = 0.1, 0.01

    dom = FL2D_small.disc(b=[1, 1, 1, 1, 1], L1=0.8, L2=0.8)

    Anₚᵣ = [64, 64, 64, 128, 128, 128, 128, 128, 128]

    AN = [3, 4, 5, 6, 7, 8, 9, 10, 12]

    f!, uex, fv = makediscfuex(5, s)

    for i in 1:9

        nₚᵣ, N = Anₚᵣ[i], AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end

    # =============================================
    dom = FL2D_small.disc(b=[2, 2, 2, 2, 2], L1=0.8, L2=0.8)

    for i in 1:9

        nₚᵣ, N = 128, AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end
    # =============================================
    dom = FL2D_small.disc(b=[3, 3, 3, 3, 3], L1=0.8, L2=0.8)

    for i in 1:9

        nₚᵣ, N = 128, AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end
    # =============================================
    dom = FL2D_small.disc(b=[5, 5, 5, 5, 5], a=[3, 3, 3, 3, 4], L1=0.8, L2=0.8)

    for i in 1:9

        nₚᵣ, N = 128, AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end

    # =============================================
    dom = FL2D_small.disc(b=[6, 6, 6, 6, 6], a=[3, 3, 3, 3, 5], L1=0.8, L2=0.8)

    for i in 1:9

        nₚᵣ, N = 256, AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=2, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end

    println(io, "==End of Disc==")

end

open("solve_outputs0250_direct.txt", "w") do io
    # =============================================
    s, p = 0.25, 4

    println(io, "Run started: ", Dates.now())
    println(io, "s = ", s)
    
    δ, δclsbd = 0.1, 0.01

    dom = FL2D_small.disc(b=[1, 1, 1, 1, 1], L1=0.8, L2=0.8)

    Anₚᵣ = [64, 64, 64, 128, 128, 128, 128, 128, 128]

    AN = [3, 4, 5, 6, 7, 8, 9, 10, 12]

    f!, uex, fv = makediscfuex(5, s)

    for i in 1:9

        nₚᵣ, N = Anₚᵣ[i], AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end

    # =============================================
    dom = FL2D_small.disc(b=[2, 2, 2, 2, 2], L1=0.8, L2=0.8)

    for i in 1:9

        nₚᵣ, N = 128, AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end
    # =============================================
    dom = FL2D_small.disc(b=[3, 3, 3, 3, 3], L1=0.8, L2=0.8)

    for i in 1:9

        nₚᵣ, N = 128, AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end
    # =============================================
    dom = FL2D_small.disc(b=[5, 5, 5, 5, 5], a=[3, 3, 3, 3, 4], L1=0.8, L2=0.8)

    for i in 1:9

        nₚᵣ, N = 128, AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end

    # =============================================
    dom = FL2D_small.disc(b=[6, 6, 6, 6, 6], a=[3, 3, 3, 3, 5], L1=0.8, L2=0.8)

    for i in 1:9

        nₚᵣ, N = 256, AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=2, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end

    println(io, "==End of Disc==")

end
# Start from here for overnight
open("solve_outputs0100_direct.txt", "w") do io
    # =============================================
    s, p = 0.1, 4

    println(io, "Run started: ", Dates.now())
    println(io, "s = ", s)
    
    δ, δclsbd = 0.1, 0.01

    dom = FL2D_small.disc(b=[1, 1, 1, 1, 1], L1=0.8, L2=0.8)

    Anₚᵣ = [64, 64, 64, 128, 128, 128, 128, 128, 128]

    AN = [3, 4, 5, 6, 7, 8, 9, 10, 12]

    f!, uex, fv = makediscfuex(5, s)

    for i in 1:9

        nₚᵣ, N = Anₚᵣ[i], AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end

    # =============================================
    dom = FL2D_small.disc(b=[2, 2, 2, 2, 2], L1=0.8, L2=0.8)

    for i in 1:9

        nₚᵣ, N = 128, AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end
    # =============================================
    dom = FL2D_small.disc(b=[3, 3, 3, 3, 3], L1=0.8, L2=0.8)

    for i in 1:9

        nₚᵣ, N = 128, AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end
    # =============================================
    dom = FL2D_small.disc(b=[5, 5, 5, 5, 5], a=[3, 3, 3, 3, 4], L1=0.8, L2=0.8)

    for i in 1:9

        nₚᵣ, N = 128, AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end

    nₚᵣ, N = 256, 12

    prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
        dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

    opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

    core_res = solveFL(prob; opts=opts)

    println(io, SolveView(prob, opts, core_res))

    flush(io)

    # =============================================
    dom = FL2D_small.disc(b=[6, 6, 6, 6, 6], a=[3, 3, 3, 3, 5], L1=0.8, L2=0.8)

    for i in 1:9

        nₚᵣ, N = 256, AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=2, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end

    println(io, "==End of Disc==")

end

open("solve_outputs0010_direct.txt", "w") do io
    # =============================================
    s, p = 0.01, 4

    println(io, "Run started: ", Dates.now())
    println(io, "s = ", s)
    
    δ, δclsbd = 0.1, 0.01

    dom = FL2D_small.disc(b=[1, 1, 1, 1, 1], L1=0.8, L2=0.8)

    Anₚᵣ = [64, 64, 64, 128, 128, 128, 128, 128, 128]

    AN = [3, 4, 5, 6, 7, 8, 9, 10, 12]

    f!, uex, fv = makediscfuex(5, s)

    for i in 1:9

        nₚᵣ, N = Anₚᵣ[i], AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end

    # =============================================
    dom = FL2D_small.disc(b=[2, 2, 2, 2, 2], L1=0.8, L2=0.8)

    for i in 1:9

        nₚᵣ, N = 128, AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end
    # =============================================
    dom = FL2D_small.disc(b=[3, 3, 3, 3, 3], L1=0.8, L2=0.8)

    for i in 1:9

        nₚᵣ, N = 128, AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end
    # =============================================
    dom = FL2D_small.disc(b=[5, 5, 5, 5, 5], a=[3, 3, 3, 3, 4], L1=0.8, L2=0.8)

    for i in 1:9

        nₚᵣ, N = 128, AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end

    nₚᵣ, N = 256, 12

    prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
        dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

    opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

    core_res = solveFL(prob; opts=opts)

    println(io, SolveView(prob, opts, core_res))

    flush(io)

    # =============================================
    dom = FL2D_small.disc(b=[6, 6, 6, 6, 6], a=[3, 3, 3, 3, 5], L1=0.8, L2=0.8)

    for i in 1:9

        nₚᵣ, N = 256, AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=2, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end

    println(io, "==End of Disc==")

end

open("solve_outputs0001_direct.txt", "w") do io
    # =============================================
    s, p = 0.001, 4

    println(io, "Run started: ", Dates.now())
    println(io, "s = ", s)
    
    δ, δclsbd = 0.1, 0.01

    dom = FL2D_small.disc(b=[1, 1, 1, 1, 1], L1=0.8, L2=0.8)

    Anₚᵣ = [64, 64, 64, 128, 128, 128, 128, 128, 128]

    AN = [3, 4, 5, 6, 7, 8, 9, 10, 12]

    f!, uex, fv = makediscfuex(5, s)

    for i in 1:9

        nₚᵣ, N = Anₚᵣ[i], AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end

    # =============================================
    dom = FL2D_small.disc(b=[2, 2, 2, 2, 2], L1=0.8, L2=0.8)

    for i in 1:9

        nₚᵣ, N = 128, AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end
    # =============================================
    dom = FL2D_small.disc(b=[3, 3, 3, 3, 3], L1=0.8, L2=0.8)

    for i in 1:9

        nₚᵣ, N = 128, AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end
    # =============================================
    dom = FL2D_small.disc(b=[5, 5, 5, 5, 5], a=[3, 3, 3, 3, 4], L1=0.8, L2=0.8)

    for i in 1:9

        nₚᵣ, N = 128, AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end

    nₚᵣ, N = 256, 12

    prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
        dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

    opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

    core_res = solveFL(prob; opts=opts)

    println(io, SolveView(prob, opts, core_res))

    flush(io)

    # =============================================
    dom = FL2D_small.disc(b=[6, 6, 6, 6, 6], a=[3, 3, 3, 3, 5], L1=0.8, L2=0.8)

    for i in 1:9

        nₚᵣ, N = 256, AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=2, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end

    println(io, "==End of Disc==")

end

open("solve_outputs1e_4_direct.txt", "w") do io
    # =============================================
    s, p = 1e-4, 4

    println(io, "Run started: ", Dates.now())
    println(io, "s = ", s)
    
    δ, δclsbd = 0.1, 0.01

    dom = FL2D_small.disc(b=[1, 1, 1, 1, 1], L1=0.8, L2=0.8)

    Anₚᵣ = [64, 64, 64, 128, 128, 128, 128, 128, 128]

    AN = [3, 4, 5, 6, 7, 8, 9, 10, 12]

    f!, uex, fv = makediscfuex(5, s)

    for i in 1:9

        nₚᵣ, N = Anₚᵣ[i], AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end

    # =============================================
    dom = FL2D_small.disc(b=[2, 2, 2, 2, 2], L1=0.8, L2=0.8)

    for i in 1:9

        nₚᵣ, N = 128, AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end
    # =============================================
    dom = FL2D_small.disc(b=[3, 3, 3, 3, 3], L1=0.8, L2=0.8)

    for i in 1:9

        nₚᵣ, N = 128, AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end
    # =============================================
    dom = FL2D_small.disc(b=[5, 5, 5, 5, 5], a=[3, 3, 3, 3, 4], L1=0.8, L2=0.8)

    for i in 1:9

        nₚᵣ, N = 128, AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end

    nₚᵣ, N = 256, 12

    prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
        dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

    opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

    core_res = solveFL(prob; opts=opts)

    println(io, SolveView(prob, opts, core_res))

    flush(io)

    # =============================================
    dom = FL2D_small.disc(b=[6, 6, 6, 6, 6], a=[3, 3, 3, 3, 5], L1=0.8, L2=0.8)

    for i in 1:9

        nₚᵣ, N = 256, AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=2, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end

    println(io, "==End of Disc==")

end

open("solve_outputs1e_5_direct.txt", "w") do io
    # =============================================
    s, p = 1e-5, 4

    println(io, "Run started: ", Dates.now())
    println(io, "s = ", s)
    
    δ, δclsbd = 0.1, 0.01

    dom = FL2D_small.disc(b=[1, 1, 1, 1, 1], L1=0.8, L2=0.8)

    Anₚᵣ = [64, 64, 64, 128, 128, 128, 128, 128, 128]

    AN = [3, 4, 5, 6, 7, 8, 9, 10, 12]

    f!, uex, fv = makediscfuex(5, s)

    for i in 1:9

        nₚᵣ, N = Anₚᵣ[i], AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end

    # =============================================
    dom = FL2D_small.disc(b=[2, 2, 2, 2, 2], L1=0.8, L2=0.8)

    for i in 1:9

        nₚᵣ, N = 128, AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end
    # =============================================
    dom = FL2D_small.disc(b=[3, 3, 3, 3, 3], L1=0.8, L2=0.8)

    for i in 1:9

        nₚᵣ, N = 128, AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end
    # =============================================
    dom = FL2D_small.disc(b=[5, 5, 5, 5, 5], a=[3, 3, 3, 3, 4], L1=0.8, L2=0.8)

    for i in 1:9

        nₚᵣ, N = 128, AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end

    # =============================================
    dom = FL2D_small.disc(b=[6, 6, 6, 6, 6], a=[3, 3, 3, 3, 5], L1=0.8, L2=0.8)

    for i in 1:9

        nₚᵣ, N = 256, AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=2, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end

    println(io, "==End of Disc==")

end

open("solve_outputs1e_6_direct.txt", "w") do io
    # =============================================
    s, p = 1e-6, 4

    println(io, "Run started: ", Dates.now())
    println(io, "s = ", s)
    
    δ, δclsbd = 0.1, 0.01

    dom = FL2D_small.disc(b=[1, 1, 1, 1, 1], L1=0.8, L2=0.8)

    Anₚᵣ = [64, 64, 64, 128, 128, 128, 128, 128, 128]

    AN = [3, 4, 5, 6, 7, 8, 9, 10, 12]

    f!, uex, fv = makediscfuex(5, s)

    for i in 1:9

        nₚᵣ, N = Anₚᵣ[i], AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end

    # =============================================
    dom = FL2D_small.disc(b=[2, 2, 2, 2, 2], L1=0.8, L2=0.8)

    for i in 1:9

        nₚᵣ, N = 128, AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end
    # =============================================
    dom = FL2D_small.disc(b=[3, 3, 3, 3, 3], L1=0.8, L2=0.8)

    for i in 1:9

        nₚᵣ, N = 128, AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end
    # =============================================
    dom = FL2D_small.disc(b=[5, 5, 5, 5, 5], a=[3, 3, 3, 3, 4], L1=0.8, L2=0.8)

    for i in 1:9

        nₚᵣ, N = 128, AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end

    # =============================================
    dom = FL2D_small.disc(b=[6, 6, 6, 6, 6], a=[3, 3, 3, 3, 5], L1=0.8, L2=0.8)

    for i in 1:9

        nₚᵣ, N = 256, AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=2, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end

    println(io, "==End of Disc==")

end

open("solve_outputs1e_7_direct.txt", "w") do io
    # =============================================
    s, p = 1e-7, 4

    println(io, "Run started: ", Dates.now())
    println(io, "s = ", s)
    
    δ, δclsbd = 0.1, 0.01

    dom = FL2D_small.disc(b=[1, 1, 1, 1, 1], L1=0.8, L2=0.8)

    Anₚᵣ = [64, 64, 64, 128, 128, 128, 128, 128, 128]

    AN = [3, 4, 5, 6, 7, 8, 9, 10, 12]

    f!, uex, fv = makediscfuex(5, s)

    for i in 1:9

        nₚᵣ, N = Anₚᵣ[i], AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end

    # =============================================
    dom = FL2D_small.disc(b=[2, 2, 2, 2, 2], L1=0.8, L2=0.8)

    for i in 1:9

        nₚᵣ, N = 128, AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end
    # =============================================
    dom = FL2D_small.disc(b=[3, 3, 3, 3, 3], L1=0.8, L2=0.8)

    for i in 1:9

        nₚᵣ, N = 128, AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end
    # =============================================
    dom = FL2D_small.disc(b=[5, 5, 5, 5, 5], a=[3, 3, 3, 3, 4], L1=0.8, L2=0.8)

    for i in 1:9

        nₚᵣ, N = 128, AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end

    # =============================================
    dom = FL2D_small.disc(b=[6, 6, 6, 6, 6], a=[3, 3, 3, 3, 5], L1=0.8, L2=0.8)

    for i in 1:9

        nₚᵣ, N = 256, AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=2, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end

    println(io, "==End of Disc==")

end

open("solve_outputs1e_8_direct.txt", "w") do io
    # =============================================
    s, p = 1e-8, 4

    println(io, "Run started: ", Dates.now())
    println(io, "s = ", s)
    
    δ, δclsbd = 0.1, 0.01

    dom = FL2D_small.disc(b=[1, 1, 1, 1, 1], L1=0.8, L2=0.8)

    Anₚᵣ = [64, 64, 64, 128, 128, 128, 128, 128, 128]

    AN = [3, 4, 5, 6, 7, 8, 9, 10, 12]

    f!, uex, fv = makediscfuex(5, s)

    for i in 1:9

        nₚᵣ, N = Anₚᵣ[i], AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end

    # =============================================
    dom = FL2D_small.disc(b=[2, 2, 2, 2, 2], L1=0.8, L2=0.8)

    for i in 1:9

        nₚᵣ, N = 128, AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end
    # =============================================
    dom = FL2D_small.disc(b=[3, 3, 3, 3, 3], L1=0.8, L2=0.8)

    for i in 1:9

        nₚᵣ, N = 128, AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end
    # =============================================
    dom = FL2D_small.disc(b=[5, 5, 5, 5, 5], a=[3, 3, 3, 3, 4], L1=0.8, L2=0.8)

    for i in 1:9

        nₚᵣ, N = 128, AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end

    # =============================================
    dom = FL2D_small.disc(b=[6, 6, 6, 6, 6], a=[3, 3, 3, 3, 5], L1=0.8, L2=0.8)

    for i in 1:9

        nₚᵣ, N = 256, AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=2, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end

    println(io, "==End of Disc==")

end

open("solve_outputs1e_9_direct.txt", "w") do io
    # =============================================
    s, p = 1e-9, 4

    println(io, "Run started: ", Dates.now())
    println(io, "s = ", s)
    
    δ, δclsbd = 0.1, 0.01

    dom = FL2D_small.disc(b=[1, 1, 1, 1, 1], L1=0.8, L2=0.8)

    Anₚᵣ = [64, 64, 64, 128, 128, 128, 128, 128, 128]

    AN = [3, 4, 5, 6, 7, 8, 9, 10, 12]

    f!, uex, fv = makediscfuex(5, s)

    for i in 1:9

        nₚᵣ, N = Anₚᵣ[i], AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end

    # =============================================
    dom = FL2D_small.disc(b=[2, 2, 2, 2, 2], L1=0.8, L2=0.8)

    for i in 1:9

        nₚᵣ, N = 128, AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end
    # =============================================
    dom = FL2D_small.disc(b=[3, 3, 3, 3, 3], L1=0.8, L2=0.8)

    for i in 1:9

        nₚᵣ, N = 128, AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end
    # =============================================
    dom = FL2D_small.disc(b=[5, 5, 5, 5, 5], a=[3, 3, 3, 3, 4], L1=0.8, L2=0.8)

    for i in 1:9

        nₚᵣ, N = 128, AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end

    # =============================================
    dom = FL2D_small.disc(b=[6, 6, 6, 6, 6], a=[3, 3, 3, 3, 5], L1=0.8, L2=0.8)

    for i in 1:9

        nₚᵣ, N = 256, AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=2, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end

    println(io, "==End of Disc==")

end

open("solve_outputs1e_10_direct.txt", "w") do io
    # =============================================
    s, p = 1e-10, 4

    println(io, "Run started: ", Dates.now())
    println(io, "s = ", s)
    
    δ, δclsbd = 0.1, 0.01

    dom = FL2D_small.disc(b=[1, 1, 1, 1, 1], L1=0.8, L2=0.8)

    Anₚᵣ = [64, 64, 64, 128, 128, 128, 128, 128, 128]

    AN = [3, 4, 5, 6, 7, 8, 9, 10, 12]

    f!, uex, fv = makediscfuex(5, s)

    for i in 1:9

        nₚᵣ, N = Anₚᵣ[i], AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end

    # =============================================
    dom = FL2D_small.disc(b=[2, 2, 2, 2, 2], L1=0.8, L2=0.8)

    for i in 1:9

        nₚᵣ, N = 128, AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end
    # =============================================
    dom = FL2D_small.disc(b=[3, 3, 3, 3, 3], L1=0.8, L2=0.8)

    for i in 1:9

        nₚᵣ, N = 128, AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end
    # =============================================
    dom = FL2D_small.disc(b=[5, 5, 5, 5, 5], a=[3, 3, 3, 3, 4], L1=0.8, L2=0.8)

    for i in 1:9

        nₚᵣ, N = 128, AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end

    # =============================================
    dom = FL2D_small.disc(b=[6, 6, 6, 6, 6], a=[3, 3, 3, 3, 5], L1=0.8, L2=0.8)

    for i in 1:9

        nₚᵣ, N = 256, AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=2, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end

    println(io, "==End of Disc==")

end

open("solve_outputs1e_11_direct.txt", "w") do io
    # =============================================
    s, p = 1e-11, 4

    println(io, "Run started: ", Dates.now())
    println(io, "s = ", s)
    
    δ, δclsbd = 0.1, 0.01

    dom = FL2D_small.disc(b=[1, 1, 1, 1, 1], L1=0.8, L2=0.8)

    Anₚᵣ = [64, 64, 64, 128, 128, 128, 128, 128, 128]

    AN = [3, 4, 5, 6, 7, 8, 9, 10, 12]

    f!, uex, fv = makediscfuex(5, s)

    for i in 1:9

        nₚᵣ, N = Anₚᵣ[i], AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end

    # =============================================
    dom = FL2D_small.disc(b=[2, 2, 2, 2, 2], L1=0.8, L2=0.8)

    for i in 1:9

        nₚᵣ, N = 128, AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end
    # =============================================
    dom = FL2D_small.disc(b=[3, 3, 3, 3, 3], L1=0.8, L2=0.8)

    for i in 1:9

        nₚᵣ, N = 128, AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end
    # =============================================
    dom = FL2D_small.disc(b=[5, 5, 5, 5, 5], a=[3, 3, 3, 3, 4], L1=0.8, L2=0.8)

    for i in 1:9

        nₚᵣ, N = 128, AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end

    # =============================================
    dom = FL2D_small.disc(b=[6, 6, 6, 6, 6], a=[3, 3, 3, 3, 5], L1=0.8, L2=0.8)

    for i in 1:9

        nₚᵣ, N = 256, AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=2, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end

    println(io, "==End of Disc==")

end

open("solve_outputs1e_12_direct.txt", "w") do io
    # =============================================
    s, p = 1e-12, 4

    println(io, "Run started: ", Dates.now())
    println(io, "s = ", s)
    
    δ, δclsbd = 0.1, 0.01

    dom = FL2D_small.disc(b=[1, 1, 1, 1, 1], L1=0.8, L2=0.8)

    Anₚᵣ = [64, 64, 64, 128, 128, 128, 128, 128, 128]

    AN = [3, 4, 5, 6, 7, 8, 9, 10, 12]

    f!, uex, fv = makediscfuex(5, s)

    for i in 1:9

        nₚᵣ, N = Anₚᵣ[i], AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end

    # =============================================
    dom = FL2D_small.disc(b=[2, 2, 2, 2, 2], L1=0.8, L2=0.8)

    for i in 1:9

        nₚᵣ, N = 128, AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end
    # =============================================
    dom = FL2D_small.disc(b=[3, 3, 3, 3, 3], L1=0.8, L2=0.8)

    for i in 1:9

        nₚᵣ, N = 128, AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end
    # =============================================
    dom = FL2D_small.disc(b=[5, 5, 5, 5, 5], a=[3, 3, 3, 3, 4], L1=0.8, L2=0.8)

    for i in 1:9

        nₚᵣ, N = 128, AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end

    # =============================================
    dom = FL2D_small.disc(b=[6, 6, 6, 6, 6], a=[3, 3, 3, 3, 5], L1=0.8, L2=0.8)

    for i in 1:9

        nₚᵣ, N = 256, AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=2, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end

    println(io, "==End of Disc==")

end

open("solve_outputs1e_13_direct.txt", "w") do io
    # =============================================
    s, p = 1e-13, 4

    println(io, "Run started: ", Dates.now())
    println(io, "s = ", s)
    
    δ, δclsbd = 0.1, 0.01

    dom = FL2D_small.disc(b=[1, 1, 1, 1, 1], L1=0.8, L2=0.8)

    Anₚᵣ = [64, 64, 64, 128, 128, 128, 128, 128, 128]

    AN = [3, 4, 5, 6, 7, 8, 9, 10, 12]

    f!, uex, fv = makediscfuex(5, s)

    for i in 1:9

        nₚᵣ, N = Anₚᵣ[i], AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end

    # =============================================
    dom = FL2D_small.disc(b=[2, 2, 2, 2, 2], L1=0.8, L2=0.8)

    for i in 1:9

        nₚᵣ, N = 128, AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end
    # =============================================
    dom = FL2D_small.disc(b=[3, 3, 3, 3, 3], L1=0.8, L2=0.8)

    for i in 1:9

        nₚᵣ, N = 128, AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end
    # =============================================
    dom = FL2D_small.disc(b=[5, 5, 5, 5, 5], a=[3, 3, 3, 3, 4], L1=0.8, L2=0.8)

    for i in 1:9

        nₚᵣ, N = 128, AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end

    # =============================================
    dom = FL2D_small.disc(b=[6, 6, 6, 6, 6], a=[3, 3, 3, 3, 5], L1=0.8, L2=0.8)

    for i in 1:9

        nₚᵣ, N = 256, AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=2, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end

    println(io, "==End of Disc==")

end

open("solve_outputs1e_14_direct.txt", "w") do io
    # =============================================
    s, p = 1e-14, 4

    println(io, "Run started: ", Dates.now())
    println(io, "s = ", s)
    
    δ, δclsbd = 0.1, 0.01

    dom = FL2D_small.disc(b=[1, 1, 1, 1, 1], L1=0.8, L2=0.8)

    Anₚᵣ = [64, 64, 64, 128, 128, 128, 128, 128, 128]

    AN = [3, 4, 5, 6, 7, 8, 9, 10, 12]

    f!, uex, fv = makediscfuex(5, s)

    for i in 1:9

        nₚᵣ, N = Anₚᵣ[i], AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end

    # =============================================
    dom = FL2D_small.disc(b=[2, 2, 2, 2, 2], L1=0.8, L2=0.8)

    for i in 1:9

        nₚᵣ, N = 128, AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end
    # =============================================
    dom = FL2D_small.disc(b=[3, 3, 3, 3, 3], L1=0.8, L2=0.8)

    for i in 1:9

        nₚᵣ, N = 128, AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end
    # =============================================
    dom = FL2D_small.disc(b=[5, 5, 5, 5, 5], a=[3, 3, 3, 3, 4], L1=0.8, L2=0.8)

    for i in 1:9

        nₚᵣ, N = 128, AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end

    # =============================================
    dom = FL2D_small.disc(b=[6, 6, 6, 6, 6], a=[3, 3, 3, 3, 5], L1=0.8, L2=0.8)

    for i in 1:9

        nₚᵣ, N = 256, AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=2, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end

    println(io, "==End of Disc==")

end

open("solve_outputs1e_15_direct.txt", "w") do io
    # =============================================
    s, p = 1e-15, 4

    println(io, "Run started: ", Dates.now())
    println(io, "s = ", s)
    
    δ, δclsbd = 0.1, 0.01

    dom = FL2D_small.disc(b=[1, 1, 1, 1, 1], L1=0.8, L2=0.8)

    Anₚᵣ = [64, 64, 64, 128, 128, 128, 128, 128, 128]

    AN = [3, 4, 5, 6, 7, 8, 9, 10, 12]

    f!, uex, fv = makediscfuex(5, s)

    for i in 1:9

        nₚᵣ, N = Anₚᵣ[i], AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end

    # =============================================
    dom = FL2D_small.disc(b=[2, 2, 2, 2, 2], L1=0.8, L2=0.8)

    for i in 1:9

        nₚᵣ, N = 128, AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end
    # =============================================
    dom = FL2D_small.disc(b=[3, 3, 3, 3, 3], L1=0.8, L2=0.8)

    for i in 1:9

        nₚᵣ, N = 128, AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end
    # =============================================
    dom = FL2D_small.disc(b=[5, 5, 5, 5, 5], a=[3, 3, 3, 3, 4], L1=0.8, L2=0.8)

    for i in 1:9

        nₚᵣ, N = 128, AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end

    # =============================================
    dom = FL2D_small.disc(b=[6, 6, 6, 6, 6], a=[3, 3, 3, 3, 5], L1=0.8, L2=0.8)

    for i in 1:9

        nₚᵣ, N = 256, AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=2, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end

    println(io, "==End of Disc==")

end
#End here
open("solve_outputs1e_16_direct.txt", "w") do io
    # =============================================
    s, p = 1e-16, 4

    println(io, "Run started: ", Dates.now())
    println(io, "s = ", s)
    
    δ, δclsbd = 0.1, 0.01

    dom = FL2D_small.disc(b=[1, 1, 1, 1, 1], L1=0.8, L2=0.8)

    Anₚᵣ = [64, 64, 64, 128, 128, 128, 128, 128, 128]

    AN = [3, 4, 5, 6, 7, 8, 9, 10, 12]

    f!, uex, fv = makediscfuex(5, s)

    for i in 1:9

        nₚᵣ, N = Anₚᵣ[i], AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end

    # =============================================
    dom = FL2D_small.disc(b=[2, 2, 2, 2, 2], L1=0.8, L2=0.8)

    for i in 1:9

        nₚᵣ, N = 128, AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end
    # =============================================
    dom = FL2D_small.disc(b=[3, 3, 3, 3, 3], L1=0.8, L2=0.8)

    for i in 1:9

        nₚᵣ, N = 128, AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end
    # =============================================
    dom = FL2D_small.disc(b=[5, 5, 5, 5, 5], a=[3, 3, 3, 3, 4], L1=0.8, L2=0.8)

    for i in 1:9

        nₚᵣ, N = 128, AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end

    # =============================================
    dom = FL2D_small.disc(b=[6, 6, 6, 6, 6], a=[3, 3, 3, 3, 5], L1=0.8, L2=0.8)

    for i in 1:9

        nₚᵣ, N = 256, AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=2, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end

    println(io, "==End of Disc==")

end

#flush(io) forces Julia to write any buffered output to the actual file immediately.