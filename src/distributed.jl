function addProcesses(n::Int = Sys.CPU_THREADS)
    n = min(Sys.CPU_THREADS, n)

    while nprocs() < n
        addprocs(1)
    end
end # function
