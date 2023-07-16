#!/usr/bin/env julia

using JLD2
using ArgParse
using StatsBase
using FileIO
using Graphs
#using Plots

mutable struct NetworkParameters

    #measurement data
    meanCoopFreq::Float64
    meanProbNeighborCoop::Float64
    meanProbNeighborDef::Float64
    meanProbRandomCoop::Float64
    meanProbRandomDef::Float64
    meanDegree::Float64
    meanAssortment::Float64
    meanCoopDefDistance::Float64
    meanDistInclusion::Float64
    meanFitness::Float64
    meanShortestPaths::Float64
    meanConnComponents::Float64

    #Node characteristics
    popPNC::Array{Float64, 1}
    popPND::Array{Float64, 1}
    popPRC::Array{Float64, 1}
    popPRD::Array{Float64, 1}
    popStrategies::Array{Int64, 1}
    popPayoff::Array{Float64, 1}
    popFitness::Array{Float64, 1}

    #structural variables
    numGens::Int64
    popSize::Int64
    edgeMatrix::Array{Int64, 2}

    #social variables
    cost::Float64
    benefit::Float64
    synergism::Float64
    linkCost::Float64
    muS::Float64
    muP::Float64    
    delta::Float64
    pnd::Bool
    prd::Bool

    #evolving links variables
    sigmapn::Float64
    sigmapr::Float64

    function NetworkParameters(b::Float64, c::Float64, d::Float64, cL::Float64, gen::Int, pn::Float64, pnd::Bool, pr::Float64, prd::Bool, muP::Float64, delta::Float64, sigmapn::Float64, sigmapr::Float64)

        popSize = 100
        popPNC = zeros(Float64, popSize)
        popPNC[:] .= pn
        popPND = zeros(Float64, popSize)
        popPND[:] .= pn
        popPRC = zeros(Float64, popSize)
        popPRC[:] .= pr
        popPRD = zeros(Float64, popSize)
        popPRD[:] .= pr
        popStrategies = zeros(Int64, popSize)
        popStrategies[2:2:popSize] .= 1
        popFitness = zeros(Float64, popSize)
        popFitness[:] .= 1.0

        #edgeMatrix = zeros(Int64, 100, 100)
        edgeMatrix=rand([0,1],(popSize,popSize))
        edgeMatrix = edgeMatrix .* transpose(edgeMatrix) #ensures that connections are reciprocated
        for(i) in 1:popSize
            edgeMatrix[i, i] = 0 #cannot connect with themself
        end

        cost = c
        synergism = d
        benefit = b
        linkCost = cL
        muS = muP #changing strategies

        new(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, popPNC, popPND, popPRC, popPRD, popStrategies, zeros(Float64, popSize), popFitness, gen, popSize, edgeMatrix, cost, benefit, synergism, linkCost, muS, muP, delta, pnd, prd, sigmapn, sigmapr)
    end
end

#measurement functions
function coopRatio(network::NetworkParameters)
    coopCount = 0
    for(i) in 1:network.popSize
        if(network.popStrategies[i] == 1)
            coopCount += 1
        end
    end
    coopCount /= network.popSize
    network.meanCoopFreq += coopCount
end

function probNeighbor(network::NetworkParameters)
    pNCTotal = 0.0
    pNDTotal = 0.0
    for(i) in 1:network.popSize
        pNCTotal += network.popPNC[i]
        pNDTotal += network.popPND[i]
    end
    pNCTotal /= network.popSize
    pNDTotal /= network.popSize
    network.meanProbNeighborCoop += pNCTotal
    network.meanProbNeighborDef += pNDTotal
end

function probRandom(network::NetworkParameters)
    pRCTotal = 0.0
    pRDTotal = 0.0
    for(i) in 1:network.popSize
        pRCTotal += network.popPRC[i]
        pRDTotal += network.popPRD[i]
    end
    pRCTotal /= network.popSize
    pRDTotal /= network.popSize
    network.meanProbRandomCoop += pRCTotal
    network.meanProbRandomDef += pRDTotal
end

function degrees(network::NetworkParameters)
    degTotal = 0
    assmtTotal = 0
    fitnessTotal = 0
    coopCount = 0.0 
    for(i) in 1:network.popSize
        if(network.popStrategies[i]==1)
            coopCount+=1.0
        end
        fitnessTotal+=network.popFitness[i]
    end
    coopCount = coopCount/network.popSize
    for(i) in 1:network.popSize
        degCounter = 0.0
        assmtCounter = 0.0
        for(ii) in 1:network.popSize
            if(network.edgeMatrix[i, ii] != 0)
                degCounter += 1.0
                for(iii) in 1:network.popSize
                    if(network.edgeMatrix[ii, iii] != 0 && network.edgeMatrix[i, iii] != 0)
                        assmtCounter += 1.0
                    end
                end
            end
        end
        degTotal += degCounter
        if(network.popStrategies[i]==1)
            assmtTotal+= (assmtCounter/degCounter)-(coopCount)
        else
            assmtTotal+= (assmtCounter/degCounter)-(1-coopCount)
        end
    end
    degTotal /= network.popSize
    fitnessTotal /= network.popSize
    assmtTotal /= network.popSize
    network.meanDegree += degTotal
    network.meanFitness += fitnessTotal
    network.meanAssortment += assmtTotal
end

function distance(network::NetworkParameters)
    distanceTotal = 0.0
    included = network.popSize
    for(i) in 1:network.popSize
        found = false
        usualSuspects = zeros(Int64, network.popSize)
        oldSuspects = zeros(Int64, network.popSize)
        distCount = 0
        usualSuspects[i] = 1
        while(!found)
            distCount += 1
            oldSuspects .= usualSuspects
            for(s) in 1:network.popSize
                if(oldSuspects[s] == 1)
                    for(ii) in 1:network.popSize
                        if(network.edgeMatrix[s, ii] != 0)
                            usualSuspects[ii] = 1
                        end
                    end
                end
            end
            if(oldSuspects == usualSuspects)
                found = true
                distCount = NaN
                included -= 1
            else
                for(ii) in 1:network.popSize
                    if(usualSuspects[ii] == 1 && network.popStrategies[ii]!=network.popStrategies[i])
                        found = true
                    end
                end
            end
        end
        if(distCount == distCount)
            distanceTotal += distCount
        end
    end
    if(included!=0)
        distanceTotal /= included
    end
    network.meanCoopDefDistance += distanceTotal
    network.meanDistInclusion += included/network.popSize
end

#evolution functions
function death(network::NetworkParameters)
    deadID = sample(1:network.popSize)
    network.edgeMatrix[deadID, :] .= 0
    network.edgeMatrix[:, deadID] .= 0
    network.popPayoff[deadID] = 0
    network.popFitness[deadID] = 1
    deadID
end

function findMom(network::NetworkParameters, kID::Int64)
    network.popFitness[kID] = 0
    fitWeights = weights(network.popFitness)
    momIndex = sample(1:network.popSize, fitWeights)
    momIndex
end

function birth(network::NetworkParameters, child::Int64, parent::Int64)
    
    network.popStrategies[child] = network.popStrategies[parent]
    if(rand()<network.muS)
        network.popStrategies[child] -= 1
        network.popStrategies[child] *= -1
    end
    
    network.popPNC[child] = network.popPNC[parent]
    if(rand()<network.muP)
        network.popPNC[child] += randn()*network.sigmapn
        network.popPNC[child] = clamp(network.popPNC[child], 0, 1)
    end
    
    if(network.pnd) #allows separate evolution of PND if pnd==true
        network.popPND[child] = network.popPND[parent]
        if(rand()<network.muP)
            network.popPND[child] += randn()*network.sigmapn
            network.popPND[child] = clamp(network.popPND[child], 0, 1)
        end
    else
        network.popPND[child] = network.popPNC[child]
    end

    network.popPRC[child] = network.popPRC[parent]
    if(rand()<network.muP)
        network.popPRC[child] += randn()*network.sigmapr
        network.popPRC[child] = clamp(network.popPRC[child], 0, 1)
    end

    if(network.prd) #allows separate evolution of PRD if prd==true
        network.popPRD[child] = network.popPRD[parent]
        if(rand()<network.muP)
            network.popPRD[child] += randn()*network.sigmapr
            network.popPRD[child] = clamp(network.popPRD[child], 0, 1)
        end
    else
        network.popPRD[child] = network.popPRC[child]
    end
    #if(rand() < network.popPNC[child])
    network.edgeMatrix[parent, child] = 1
    network.edgeMatrix[child, parent] = 1
    #end 
    #Pb != 1 now

    for(i) in 1:network.popSize
        if(i != child && network.edgeMatrix[i, child] == 0)
            if(network.edgeMatrix[i, parent] != 0)
                if(network.popStrategies[i] == 1) #PNC MODE
                    if(rand() < network.popPNC[child])
                        network.edgeMatrix[i, child] = 1
                        network.edgeMatrix[child, i] = 1
                    end
                else
                    if(rand() < network.popPND[child]) #PND MODE
                        network.edgeMatrix[i, child] = 1
                        network.edgeMatrix[child, i] = 1
                    end
                end
            else
                if(network.popStrategies[i] == 1) #PRC MODE
                    if(rand() < network.popPRC[child])
                        network.edgeMatrix[i, child] = 1
                        network.edgeMatrix[child, i] = 1
                    end
                else
                    if(rand() < network.popPRD[child]) #PRD MODE
                        network.edgeMatrix[i, child] = 1
                        network.edgeMatrix[child, i] = 1
                    end
                end
                #if(rand() < network.popPR[child]) #original code
                network.edgeMatrix[i, child] = 1
                network.edgeMatrix[child, i] = 1
                #end
            end
        end
    end
end

function cooperate(network::NetworkParameters)
    degs = getDegree(network)
    for(i) in 1:network.popSize
        if(network.popStrategies[i] == 1)
            B = network.benefit/degs[i]
            D = network.synergism/degs[i]

            for(ii) in 1:network.popSize
                if(network.edgeMatrix[i, ii] != 0)
                    network.popPayoff[ii] += B
                    if(network.popStrategies[ii] == 1)
                        network.popPayoff[ii] += D/degs[ii]
                    end
                end
            end

            network.popPayoff[i] -= network.cost
        end
        network.popPayoff[i] -= network.linkCost * degs[i]
    end
end

function resolveFitnesses(network::NetworkParameters)
    for(i) in 1:network.popSize
        network.popFitness[i] = (1.0 + network.delta) ^ network.popPayoff[i]
    end
    network.popPayoff[:] .= 0.0
end

function getDegree(network::NetworkParameters) #made less efficient by 2 in edgeMatrix; can use sum for 1/0 vals
    degGetter = zeros(100)
    for(i) in 1:100
        for(ii) in 1:100
            if(network.edgeMatrix[i,ii] != 0)
                degGetter[i] += 1
            end
        end
        #degGetter[i] = sum(edgeMatrix[i,:])
    end

    degGetter
end

function graphCalc(network::NetworkParameters) #computes all calculations involving Graphs.jl, including shortest paths & connected components
    edgeGraph = Graph(network.edgeMatrix)

    allPath = 0.0
    disconn = 0 #keeps track of disconnected components
    for(i) in 1:network.popSize
        dsp = dijkstra_shortest_paths(edgeGraph, i)
        paths = dsp.dists
        undefs = zeros(Int, 0) #keeps track of infinite and zero distances between nodes
        for(j) in 1:(size(paths)[1])
            if(paths[j] > (network.popSize+1) || paths[j] == 0)
                push!(undefs, j)
            end
        end
        if(size(undefs)==size(paths)) #if all shortest paths are zero or infinite, disconn increments
            disconn += 1
        else
            for(k) in 1:size(undefs)[1]
                deleteat!(paths, undefs[k]) #removes all infinite and zero distances from paths
                undefs .-= 1 #avoids index out of bounds
            end
            allPath += mean(paths)
        end
    end
    network.meanShortestPaths += (allPath / (network.popSize-disconn))

    network.meanConnComponents += size(connected_components(edgeGraph))[1] #connected_components(edgeGraph) returns a vector of vectors of each connected component
end

function runSimsReturn(;B::Float64=2.0, C::Float64=0.5, D::Float64=0.0, CL::Float64=0.0, gen::Int=500, pn::Float64=0.5, pnd::Bool=false, pr::Float64=0.01, prd::Bool=false, muP::Float64=0.001, delta::Float64=0.1, sigmapn::Float64=0.05, sigmapr::Float64=0.01, reps::Int64=50)
    dataArray = zeros(10) 
    repSims = reps
    for(x) in 1:repSims

        #initializes globalstuff structure with generic constructor
        network = NetworkParameters(B, C, D, CL, gen, pn, pnd, pr, prd, muP, delta, sigmapn, sigmapr) #if pnd/prd = true, then defector & cooperator probabilities will evolve separately for pn/pr

        #checks efficiency of simulation while running it
        for(g) in 1:(network.numGens * network.popSize)

            childID = death(network)
            parentID = findMom(network, childID)
            birth(network, childID, parentID)
            if(g > (20))
                cooperate(network)
            end
            #cooperate(network)
            resolveFitnesses(network)

            if(g > (network.numGens * network.popSize / 5) && (g % network.popSize) == 0)
                coopRatio(network)
                probNeighbor(network)
                probRandom(network)
                degrees(network)
                graphCalc(network)
                #distance(network)
            end

        end

        #divides meanCooperationRatio by last 400 generations to get a true mean, then outputs
        network.meanProbNeighborCoop /= (network.numGens*0.8)
        network.meanProbNeighborDef /= (network.numGens*0.8)
        network.meanProbRandomCoop /= (network.numGens*0.8)
        network.meanProbRandomDef /= (network.numGens*0.8)
        network.meanDegree /= (network.numGens*0.8)
        network.meanAssortment /= (network.numGens*0.8)
        network.meanCoopFreq /= (network.numGens*0.8)
        network.meanCoopDefDistance /= (network.popSize*network.numGens*0.8)
        network.meanDistInclusion /= (network.popSize*network.numGens*0.8)
        network.meanFitness /= (network.numGens*0.8)
        network.meanShortestPath /= (network.numGens*0.8)
        network.meanConnComponents /= (network.numGens*0.8)

        dataArray[1] += network.meanProbNeighborCoop
        dataArray[2] += network.meanProbNeighborDef
        dataArray[3] += network.meanProbRandomCoop
        dataArray[4] += network.meanProbRandomDef
        dataArray[5] += network.meanDegree
        dataArray[6] += network.meanAssortment
        dataArray[7] += network.meanCoopDefDistance
        dataArray[8] += network.meanDistInclusion
        dataArray[9] += network.meanCoopFreq
        dataArray[10] += network.meanFitness
        dataArray[11] += network.meanShortestPath
        dataArray[12] += network.meanConnComponents
    end
    dataArray[:] ./= Float64(repSims)
    return dataArray
    #save("sim_PNCD$(pnc)_$(pnd)_PR$(pr)_CL$(CL)_B$(BEN)_G$(gen).jld2", "parameters", [CL, BEN], "meanPNI", dataArray[1], "meanPNR", dataArray[2], "meanPR", dataArray[3], "meanDegree", dataArray[4], "meanAssortment", dataArray[5], "meanDistanceFromDefToCoop", dataArray[6], "meanDistanceInclusion", dataArray[7], "meanCooperationRatio", dataArray[8])
end

#=profiling
using Profile
#net = NetworkParameters(1.0,0.1)
Profile.clear()
runSims(0.1, 1.0)
@profile runSims(0.1, 1.0)
#Profile.print()
#runSims(0.1, 1.0)=#

#test = Graph(edgeMatrix)
#gplot(test, nodelabel=1:10)
