using JLD2
using ArgParse
using StatsBase
using FileIO
using Graphs
using GraphPlot

popSize = 6
edgeMatrix=rand([0,1],(popSize,popSize))
        edgeMatrix = edgeMatrix .* transpose(edgeMatrix) #ensures that connections are reciprocated
        for(i) in 1:popSize
            edgeMatrix[i, i] = 0 #cannot connect with themself
        end

edgeGraph = Graph(edgeMatrix)
    
allPath = 0.0
disconn = 0 #keeps track of disconnected components
    for(i) in 1:popSize
        dsp = dijkstra_shortest_paths(edgeGraph, i)
        paths = dsp.dists
        undefs = zeros(Int, 0) #keeps track of infinite and zero distances between nodes
        for(j) in 1:(size(paths)[1])
            if(paths[j] > (popSize+1) || paths[j] == 0)
                push!(undefs, j)
            end
        end
        print(i, " Undefs: ", undefs)
        if(size(undefs)==size(paths)) #if all shortest paths are zero or infinite, disconn increments
            disconn += 1
            println(" Disconnected Node")
        else
            for(k) in 1:size(undefs)[1]
                deleteat!(paths, undefs[k]) #removes all infinite and zero distances from paths
                undefs .-= 1 #avoids index out of bounds
            end
            println(" Paths: ", paths)
            allPath += mean(paths)
        end
    end
    allPath /= (popSize-disconn)
    
println(allPath)
gplot(edgeGraph, nodelabel=1:popSize)

meanConnComponents = 0
meanConnComponents += size(connected_components(edgeGraph))[1]
println(meanConnComponents)

function get_absolute_time()
    s = 0.0
    for i=1:10000000
        s += i
    end
    sleep(1)
end

@time get_absolute_time()

popSize = 100
popLocations = zeros(Float64, popSize)
            for(i) in 1:popSize
                popLocations[i] = i #sets popLocations to [1, 2, ... popSize]
            end
            edgeMatrix = zeros(Float64, popSize, popSize)
            for(i) in 1:popSize
                linkWeights = (((popLocations[i].-popLocations)).^(2)).^(-0.25) #inverse of the sqrt of distance between locations
                links = zeros(Float64, Int(popSize/2))
                sample!(popLocations, Weights(linkWeights), links, replace=false)
                for(j) in links
                    edgeMatrix[i, Int(j)] = 1    
                end
            end  

    dists = abs.(popLocations.-popLocations[1]) #.^(2).^(0.5) 
    println(dists)
    for(i) in 1:popSize
        if(dists[i] < 6)
            dists[i] = 1 #all distances 5 or under are set to 1
        else
            dists[i] -= 5
        end
    end
    dists = dists.^(-0.5)
    
println(dists)
