using JLD2
using StatsBase
using FileIO
using ArgParse
using Plots
using DataFrames
using CSV
using Distributed
backend(:plotly)

addprocs(40) 
@everywhere include("ONS_Fixed_Links.jl")

@everywhere function run_worker(inputs, results)
    include("ONS_Fixed_Links.jl")
    while true
        pard = take!(inputs)
        println(pard["ben"], " ", pard["cl"], " in pard")
        coopFreq = runSimsReturn(; B=pard["ben"], C=0.5, D=0.0, CL=pard["cl"], gen=100000, pnc=0.5, pnd=0.5, pr=0.0001, muP=0.001, delta=0.1, sigmapn=0.01, sigmapr=0.01, reps=10)
        #println(pard["pn"], " ", pard["pr"], " CF: ", coopFreq[8])
        Keys = ["pnc_end","pnd_end","pr_end","degree","assortment","distance","inclusion","coopFreq"]
        temp = Dict(zip(Keys, coopFreq))
        temp = merge(pard, temp)
        println(temp["ben"], " ", temp["cl"], " CF: ", temp["coopFreq"])
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
range = 20        
inputs  = RemoteChannel(()->Channel{Dict}(range*range)) #2*nsets*maximum(pars["num_crossings"])
results = RemoteChannel(()->Channel{Dict}(range*range))
        
vals_arr = Array{Tuple}
pars = Dict([
        "ben" => 0.0,
        "cl" => 0.0,
        "pnc_end" => 0.0,
        "pnd_end" => 0.0,
        "pr_end" => 0.0,
        "degree" => 0.0,
        "assortment" => 0.0,
        "distance" => 0.0,
        "inclusion" => 0.0,
        "coopFreq" => 0.0,
    ])
nruns = fill_inputs(range,pars, 0)

for w in workers() # start tasks on the workers to process requests in parallel
    remote_do(run_worker, w, inputs, results)
end

file = "fig_4_small_sigmapn.csv"
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
end#hmap(results, 8, "Cooperation Frequencies")

#save("parker-hmap1.jld", "matr", data)
#currentDict = load("akcay-hmap1.jld")
#data2 = currentDict["matr"]
#diff = data[:, :, 8] - data2[:, :, 1]
#hmap(diff, 10, 1)

