using JLD2
using StatsBase
using FileIO
using ArgParse
using DataFrames
using CSV
using Distributed

addprocs(40) 
@everywhere include("ONS_Fixed_Links.jl")

@everywhere function run_worker(inputs, results)
    include("ONS_Fixed_Links.jl")
    while true
        pard = take!(inputs)
        println(pard["ben"], " ", pard["cl"], " in pard")
        coopFreq = runSimsReturn(; B=pard["ben"], C=0.5, D=0.0, CL=pard["cl"], gen=100000, dbOrder=deathbirth, findMom=anyMom, pn=0.5, distInherit=false, pnd=false, pr=0.0001, prd=false, muP=0.01, delta=0.5, sigmapn=0.01, sigmapr=0.01, reps=10)
        Keys = ["pnc_end","pnd_end","prc_end","prd_end","degree","assortment","distance","inclusion","coopFreq","fitness","shortestPath","connComponents","meanConnCompSize","largestConnComp","meanConnDistance"]
        temp = Dict(zip(Keys, coopFreq))
        temp = merge(pard, temp)
        println(temp["ben"], " ", temp["cl"], " CF: ", round(temp["coopFreq"]; digits = 3), " PNC: ", round(temp["pnc_end"]; digits = 3), " PND: ", round(temp["pnd_end"]; digits = 3), " PRC: ", round(temp["prc_end"]; digits = 3), " PRD: ", round(temp["prd_end"]; digits = 3) )
        put!(results, temp)
    end
end

function fill_inputs(range,pars, nruns)
    for(i) in 1:range #PNC/PND 
        for(j) in 1:range #PNR 
            vals = (10*i/range, 0.4*j/range)
            temp = copy(pars)
            temp["ben"] = vals[1]
            temp["cl"] = vals[2]
            println(temp["ben"], " ", temp["cl"], "into inputs")
            nruns+=1
            put!(inputs, temp)
            #vals_arr.push(vals)
        end
    end
    return nruns
end

#notebook for running below
@time begin
    range = 20        
    inputs  = RemoteChannel(()->Channel{Dict}(range*range)) #2*nsets*maximum(pars["num_crossings"])
    results = RemoteChannel(()->Channel{Dict}(range*range))
            
    vals_arr = Array{Tuple}
    pars = Dict([
            "ben" => 0.0,
            "cl" => 0.0,
            "pnc_end" => 0.0,
            "pnd_end" => 0.0,
            "prc_end" => 0.0,
            "prd_end" => 0.0,
            "degree" => 0.0,
            "assortment" => 0.0,
            "distance" => 0.0,
            "inclusion" => 0.0,
            "coopFreq" => 0.0,
            "fitness" => 0.0,
            "shortestPath" => 0.0,
            "connComponents" => 0.0,
            "meanConnCompSize" => 0.0,
            "largestConnComp" => 0.0,
            "meanConnDistance" => 0.0,
        ])
    nruns = fill_inputs(range,pars, 0)

    for w in workers() # start tasks on the workers to process requests in parallel
        remote_do(run_worker, w, inputs, results)
    end


    file = "anyMom_no_pnd_new.csv"
        cols = push!(sort(collect(keys(pars))),
                    ["ben", "cl"]...)
        dat = DataFrame(Dict([(c, Any[]) for c in cols]))

    for sim in 1:nruns
        # get results from parallel jobs
        flush(stdout)
        resd = take!(results)
        # add to table (must convert dict keys to symbols) and save
        push!(dat, Dict([(Symbol(k), resd[k]) for k in keys(resd)]))
        CSV.write(file, dat)
    end

end
println("Data Collected")
