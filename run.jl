include("ONS_Fixed_Links.jl")

argTab = ArgParseSettings(description = "arguments and stuff, don't worry about it")
@add_arg_table argTab begin
    "--cLink"
        arg_type = Float64
        default = 0.0
    "--gens"
        arg_type = Int
        default = 500
    "--pnc"
        arg_type = Float64
        default = 0.5
    "--pnd"
        arg_type = Float64
        default = 0.5
    "--pr"
        arg_type = Float64
        default = 0.001
end
parsedArgs = parse_args(ARGS, argTab)
currCostLink = parsedArgs["cLink"]
currGens = parsedArgs["gens"]
currPNC = parsedArgs["pnc"]
currPD = parsedArgs["pnd"]
currPR = parsedArgs["pr"]
for(b) in 2:3
    currBenefit = Float64(b)
    runSims(currBenefit, currCostLink, currGens, currPNC, currPD, currPR)
end