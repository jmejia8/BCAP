function addProcesses(n::Int = Sys.CPU_THREADS)
    n = min(Sys.CPU_THREADS, n)

    while nprocs() > n
        t = rmprocs(nprocs(), waitfor=0)
        wait(t)
    end

    while nprocs() < n
        addprocs(1)
    end
end # function
