#!/usr/bin/env julia

using JLD2
using ArgParse
using StatsBase
using FileIO
#using Plots

mutable struct NetworkParameters

    #measurement data
    meanCoopFreq::Float64
    meanProbNeighborCoop::Float64
    meanProbNeighborDef::Float64
    meanProbRandom::Float64
    meanDegree::Float64
    meanAssortment::Float64
    meanCoopDefDistance::Float64
    meanDistInclusion::Float64

    #Node characteristics
    popPNC::Array{Float64, 1}
    popPND::Array{Float64, 1}
    popPR::Array{Float64, 1}
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
    mu::Float64
    delta::Float64

    function NetworkParameters(b::Float64, cL::Float64, gen::Int)

        numGens = gen
        popSize = 100
        popPNC = zeros(Float64, popSize)
        popPNC[:] .= 0.5
        popPND = zeros(Float64, popSize)
        popPND[:] .= 0.5
        popPR = zeros(Float64, popSize)
        popPR[:] .= .0001
        popStrategies = zeros(Int64, popSize)
        popStrategies[2:2:popSize] .= 1
        popFitness = zeros(Float64, popSize)
        popFitness[:] .= 1.0

        edgeMatrix = zeros(Int64, 100, 100)
        cost = 0.5
        synergism = 0.0
        benefit = b
        linkCost = cL
        mu = .01
        delta = 0.1 #EDIT 0.5

        new(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, popPNC, popPND, popPR, popStrategies, zeros(Float64, popSize), popFitness, numGens, popSize, edgeMatrix, cost, benefit, synergism, linkCost, mu, delta)
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
    pRTotal = 0.0
    for(i) in 1:network.popSize
        pRTotal += network.popPR[i]
    end
    pRTotal /= network.popSize
    network.meanProbRandom += pRTotal
end

function degrees(network::NetworkParameters)
    degTotal = 0
    assmtTotal = 0
    coopCount = 0.0
    for(i) in 1:network.popSize
        if(network.popStrategies[i]==1)
            coopCount+=1.0
        end
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
    assmtTotal /= network.popSize
    network.meanDegree += degTotal
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
    if(rand()<network.mu)
        network.popStrategies[child] -= 1
        network.popStrategies[child] *= -1
    end
    network.popPNC[child] = network.popPNC[parent]
    if(rand()<network.mu)
        network.popPNC[child] += randn()/100
        network.popPNC[child] = clamp(network.popPNC[child], 0, 1)
    end

    network.popPND[child] = network.popPND[parent]
    if(rand()<network.mu)
        network.popPND[child] += randn()/100
        network.popPND[child] = clamp(network.popPND[child], 0, 1)
    end

    network.popPR[child] = network.popPR[parent]
    if(rand()<network.mu)
        network.popPR[child] += randn()/100
        network.popPR[child] = clamp(network.popPR[child], 0, 1)
    end
    if(rand() < network.popPNC[child])
        network.edgeMatrix[parent, child] = 1
        network.edgeMatrix[child, parent] = 1
    end #Pb != 1 now

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
                if(rand() < network.popPR[child])
                    network.edgeMatrix[i, child] = 2
                    network.edgeMatrix[child, i] = 2
                end
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

function runSims(CL::Float64, BEN::Float64, gen::Int)
    dataArray = zeros(8)
    repSims = 10
    for(x) in 1:repSims

        #initializes globalstuff structure with generic constructor
        network = NetworkParameters(BEN, CL, gen)

        #checks efficiency of simulation while running it
        for(g) in 1:(network.numGens * network.popSize)

            childID = death(network)
            parentID = findMom(network, childID)
            birth(network, childID, parentID)
            if(g > (network.numGens * network.popSize / 5))
                cooperate(network)
            end
            resolveFitnesses(network)

            if(g > (network.numGens * network.popSize / 5) && (g % network.popSize) == 0)
                coopRatio(network)
                #probNeighbor(network)
                #probRandom(network)
                degrees(network)
                #distance(network)
            end

        end

        #divides meanCooperationRatio by last 400 generations to get a true mean, then outputs
        network.meanProbNeighborCoop /= (network.popSize*network.numGens*0.8)
        network.meanProbNeighborDef /= (network.popSize*network.numGens*0.8)
        network.meanProbRandom /= (network.popSize*network.numGens*0.8)
        network.meanDegree /= (network.popSize*network.numGens*0.8)
        network.meanAssortment /= (network.popSize*network.numGens*0.8)
        network.meanCoopFreq /= (network.popSize*network.numGens*0.8)
        network.meanCoopDefDistance /= (network.popSize*network.numGens*0.8)
        network.meanDistInclusion /= (network.popSize*network.numGens*0.8)

        dataArray[1] += network.meanProbNeighborCoop
        dataArray[2] += network.meanProbNeighborDef
        dataArray[3] += network.meanProbRandom
        dataArray[4] += network.meanDegree
        dataArray[5] += network.meanAssortment
        dataArray[6] += network.meanCoopDefDistance
        dataArray[7] += network.meanDistInclusion
        dataArray[8] += network.meanCoopFreq
    end
    dataArray[:] ./= Float64(repSims)
    #EDIT NAME
    save("newAssortmentCD_CL$(CL)_B$(BEN)_G$(gen).jld2", "parameters", [CL, BEN], "meanPNI", dataArray[1], "meanPNR", dataArray[2], "meanPR", dataArray[3], "meanDegree", dataArray[4], "meanAssortment", dataArray[5], "meanDistanceFromDefToCoop", dataArray[6], "meanDistanceInclusion", dataArray[7], "meanCooperationRatio", dataArray[8])
end

argTab = ArgParseSettings(description = "arguments and stuff, don't worry about it")
@add_arg_table argTab begin
    "--cLink"
        arg_type = Float64
        default = 0.0
    "--gens"
        arg_type = Int
        default = 500
end
parsedArgs = parse_args(ARGS, argTab)
currCostLink = parsedArgs["cLink"]
currGens = parsedArgs["gens"]
for(b) in 2:3
    currBenefit = Float64(b)
    runSims(currCostLink, currBenefit, currGens)
end


#=profiling
using Profile
#net = NetworkParameters(1.0,0.1)
Profile.clear()
runSims(0.1, 1.0)
@profile runSims(0.1, 1.0)
#Profile.print()
#runSims(0.1, 1.0)=#