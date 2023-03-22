using JLD2
using StatsBase
using FileIO
using ArgParse
using Plots
backend(:plotly)

include("ONS_Fixed_Links.jl")

function dataCollect(range::Int) #collects data from 0.1 - range/10 for PNC/D and 0.01 - range/100 for PR
    freqs = zeros(range, range, 8)
    for(i) in 1:range #PNC/PND 
        for(j) in 1:range #PNR 
            vals = (i/(range), j/(20*range)) #easier access
            coopFreq = runSimsReturn(B=2.0, C=0.5, gen=500, pnc=vals[1], pnd=vals[1], pr=vals[2], reps=50)
            #function run0SimsReturn(BEN::Float64, CL::Float64, gen::Int, pnc::Float64, pnd::Float64, pr::Float64)
            println(coopFreq[8], " ", round(vals[1]; digits = 3), "-PN, ", round(vals[2]; digits = 3), "-PR")
            freqs[i, j, :] = coopFreq
        end
    end
    return freqs
end

function hmap(data, range::Int, index::Int) #index 8 is coop frequency
    x_axis = String[]
    y_axis = String[]
    for(i) in 1:range #PNC/PND times 10
        push!(x_axis, (string(round(i/(20*range); digits = 3))))
        push!(y_axis, (string(round(i/(range); digits = 3))))
    end
    p = heatmap(x_axis, y_axis, data[:, :, index]; title = "Cooperation Frequencies",)
    gui(p)
end

#notebook for running below
data = dataCollect(10)
hmap(data, 10, 8)
