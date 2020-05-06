function hill_climbing(pre_hill, k_level, num_inter, step)
    #step is used to keep track of all the rows we have checked
    #  so we don't waste resources rechecking them
    #Iterate through every row in the array
    for i in step:length(pre_hill)
        #if it is found that a row can be removed...
        # println("aaa")
        # @show i
        tmp = pre_hill[i]
        if check_interactions(deleteat!(pre_hill, i), k_level, num_inter)
            # println("bbb")
            #print("We have to go deeper..." + str(i))
            #...recurse the program on itself with said row removed
            @show length(pre_hill)
            return hill_climbing(pre_hill, k_level, num_inter, i)
        end
        insert!(pre_hill, i, tmp)
        # println("ccc")
    end
    return pre_hill
end

function check_interactions(test, k_level, num_inter)
    #get all existing combinations (interactions) from the data
    for p in IterTools.subsets(1:length(k_level), num_inter)
        #Mul represents the number of unique iterations we expect for the exerpt
        mul = 1
        for i in p
            mul *= k_level[i]
        end
        #check if every possible interaction is represented
        vv = length( unique(map( t -> t[p], test)) )
        if vv < mul
            #if not, return false
            return false
        end
    end
    #if we made it thro
    return true
end

function create_cover(k_level, num_inter)
    #get the full array that we can take random elements out of
    full = create_full_array(k_level)

    #get the predicted covering array size
    cov_size = min(find_size(num_inter, length(k_level), maximum(k_level)), length(full))

    running = 100
    test = []
    working_test = []

    while (true)
        #Create a new test array
        test = create_test_array(full, cov_size)
        #check if all interactions are represented in this test array
        if (check_interactions(test, k_level, num_inter))
            #save the working test
            #cov_size = cov_size - 1
            working_test = test
            break
        else
            #else, iterate down
            running = running - 1
            #if there have been enough failed attempts
            #  (usually happens on larger arrays)
            #increase the size of the covariance matrix and try again
            if (running == 0)
                cov_size = min(mul_array(k_level), cov_size + 1)
                running = 100
            end
        end
    end

    return working_test
end

function mul_array(k_level)
    mul = 1
    for i = 1:length(k_level)
        mul = mul * k_level[i]
    end
    return mul
end

function create_full_array(input_array)
    #Create an array in a format that can be passed to itertools
    give = []
    for i in input_array
        push!(give, 1:i)
    end
    #Pass said array to itertools

    # return Iterators.product(give...)
    arr = Tuple[]
    I = Iterators.product(give...)
    # display(I)
    for i in I
        push!(arr, i)
    end
    return arr
end



function create_test_array(array, size_)
    #Sorted for readability
    #random.sample to prevent duplicate elements
    #min to make sure that we don't run into a sample error
    ii = sort(sample(1:length(array), size_, replace = false))
    return array[ii]

end



function find_size(t, k, v)
    top = log(choose(k, t)) + (t * log(v))
    bottom = log((v^t) / (v^t - 1))
    return floor(Int, top / bottom) + 1
end


function choose(n, k)
    """
    A fast way to calculate binomial coefficients by Andrew Dalke (contrib).
    """
    if 0 <= k <= n
        ntok = 1
        ktok = 1
        for t = 1:min(k, n - k)+1
            ntok *= n
            ktok *= t
            n -= 1
        end
        return ntok ÷ ktok
    end

    return 0
end


function main()
    k_level = repeat([2], 10)
    # k_level[end] = 20
    num_inter = 2


    k_level = sort(k_level)

    CA = create_cover(k_level, num_inter)
    @show length(CA)
    hill_climbing(CA, k_level, num_inter, 1)
end

@time r = main()

println("-------------------")
display(r)
