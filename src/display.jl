import Base.Multimedia.display


function display(ll::LLSolution)
    println("")
    ll.isfeasible ? println("Feasible solution") : println("Infeasible solution")
    println("instance_values: ")
    display(ll.instance_values)
    @show ll.solved_instances
    @show ll.evaluated_instances
    @show ll.seed
end
