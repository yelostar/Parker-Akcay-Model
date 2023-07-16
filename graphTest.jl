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
