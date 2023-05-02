using JLD2
using StatsBase
using FileIO
using ArgParse
using Plots
using Distributed
backend(:plotly)

addprocs(10) 
@everywhere include("ONS_Fixed_Links.jl")


@everywhere function run_worker(inputs, results)
    include("ONS_Fixed_Links.jl")
    while true
        pard = take!(inputs)
        println(pard["pn"], " ", pard["pr"], "in pard")
        coopFreq = runSimsReturn(; B=2.0, C=0.5, D=0.0, CL=0.0, gen=500, pnc=pard["pn"], pnd=pard["pn"], pr=pard["pr"], muP=0.001, reps=100)
        println(pard["pn"], " ", pard["pr"], " CF: ", coopFreq[8])
        pard["data"] = coopFreq
        put!(results, pard)
    end
end

function fill_inputs(range)
    pard = Dict(zip(keys(pars), zeros(8)))
    for(i) in 1:range #PNC/PND 
        for(j) in 1:range #PNR 
            vals = (i/(range), j/(20*range))
            temp = copy(pard)
            temp["pn"] = vals[1]
            temp["pr"] = vals[2]
            println(temp["pn"], " ", temp["pr"], "into inputs")
            put!(inputs, temp)
            #vals_arr.push(vals)
        end
    end
end

function hmap(results, index::Int, chartTitle::String) #index 8 is coop frequency
    x_axis = String[]
    y_axis = String[]
    for vec in results #PNC/PND times 10
        push!(x_axis, (string(round(vec["pr"]; digits = 3))))
        push!(y_axis, (string(round(vec["pn"]; digits = 3))))
    end
    p = heatmap(x_axis, y_axis, results["data"][index]; title = chartTitle,)
    gui(p)
end

#notebook for running below
        
inputs  = RemoteChannel(()->Channel{Dict}(4000)) #2*nsets*maximum(pars["num_crossings"])
results = RemoteChannel(()->Channel{Dict}(4000))
        
vals_arr = Array{Tuple}
pars = Dict{String,Any}([
        "pn"     => Dict("value" => 0, "type" => Float64),
        "pr" => Dict("value" => 0, "type" => Float64),
        "data" => Dict("value" => zeros(8), "type" => Vector{Float64}),
    ])
        
fill_inputs(10)

for w in workers() # start tasks on the workers to process requests in parallel
    remote_do(run_worker, w, inputs, results)
end

#hmap(results, 8, "Cooperation Frequencies")

#save("parker-hmap1.jld", "matr", data)
#currentDict = load("akcay-hmap1.jld")
#data2 = currentDict["matr"]
#diff = data[:, :, 8] - data2[:, :, 1]
#hmap(diff, 10, 1)

#dataArray[1] += network.meanProbNeighborCoop
#dataArray[2] += network.meanProbNeighborDef
#dataArray[3] += network.meanProbRandom
#dataArray[4] += network.meanDegree
#dataArray[5] += network.meanAssortment
#dataArray[6] += network.meanCoopDefDistance
#dataArray[7] += network.meanDistInclusion
#dataArray[8] += network.meanCoopFreq
