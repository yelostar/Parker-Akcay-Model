using JLD2
using ArgParse
using StatsBase
using FileIO
using Graphs

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
    meanConnCompSize::Float64
    meanLargestConnComp::Float64
    meanDistConnection::Float64

    #Node characteristics
    popPNC::Array{Float64, 1}
    popPND::Array{Float64, 1}
    popPRC::Array{Float64, 1}
    popPRD::Array{Float64, 1}
    popStrategies::Array{Int64, 1}
    popPayoff::Array{Float64, 1}
    popFitness::Array{Float64, 1}
    popLocations::Array{Float64, 1}

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

    #simulation determining
    pnd::Bool
    prd::Bool
    distInherit::Bool

    #evolving links variables
    sigmapn::Float64
    sigmapr::Float64

    #distance inherit variables
    distFactor::Float64

    function NetworkParameters(b::Float64, c::Float64, d::Float64, cL::Float64, gen::Int, distInherit::Bool, distFactor::Float64, pn::Float64, pnd::Bool, pr::Float64, prd::Bool, muP::Float64, delta::Float64, sigmapn::Float64, sigmapr::Float64)

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
        
        popLocations = zeros(Float64, popSize)
        for(i) in 1:popSize
            popLocations[i] = i #sets popLocations to [1, 2, ... popSize]
        end
        if(distInherit==true)
            edgeMatrix = zeros(popSize, popSize)
            for(i) in 1:popSize
                linkWeights = distFactor.^(min.(abs.(popLocations.-popLocations[i]), abs.(popLocations.-popLocations[i].-popSize),abs.(popLocations.-popLocations[i].+popSize))) #creates weighted list of locations based on distFactor
                links = zeros(Float64, Int(popSize/2)) 
                sample!(popLocations, Weights(linkWeights), links, replace=false) #picks 50 links for each individual with probabilities weighted by locations
                for(j) in links
                    edgeMatrix[i, Int(j)] = 1    
                end
            end       
        else
            edgeMatrix=rand([0,1],(popSize,popSize))
        end
        edgeMatrix = edgeMatrix .* transpose(edgeMatrix) #ensures that connections are reciprocated
        for(i) in 1:popSize
            edgeMatrix[i, i] = 0 #cannot connect with themself
        end

        cost = c
        synergism = d
        benefit = b
        linkCost = cL
        muS = muP #changing strategies
        new(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, popPNC, popPND, popPRC, popPRD, popStrategies, zeros(Float64, popSize), popFitness, popLocations, gen, popSize, edgeMatrix, cost, benefit, synergism, linkCost, muS, muP, delta, pnd, prd, distInherit, sigmapn, sigmapr, distFactor)
    end
end

function getDegree(network::NetworkParameters) #2 in edgeMatrix removed
    degGetter = zeros(100)
    for(i) in 1:100
        degGetter[i] = sum(network.edgeMatrix[i, :])
    end

    degGetter
end

function cooperate(network::NetworkParameters)
    degs = getDegree(network)
    for(i) in 1:network.popSize
        if(network.popStrategies[i] == 1)
            B = network.benefit/degs[i]
            D = network.synergism/degs[i]
            network.popPayoff.+= B.*network.edgeMatrix[i, :]
            network.popPayoff.+= (D.*network.edgeMatrix[i, :].*network.popStrategies./degs)
            network.popPayoff[i] -= network.cost
        end
        network.popPayoff[i] -= network.linkCost * degs[i]
    end
end

function oldcoop(network::NetworkParameters)
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

function rfold(network::NetworkParameters)
    for(i) in 1:network.popSize
        network.popFitness[i] = (1.0 + network.delta) ^ network.popPayoff[i]
    end
    network.popPayoff[:] .= 0.0
end

function resolveLocs(network::NetworkParameters, child, parent)
    deathDist = network.popLocations[child]-network.popLocations[parent]

    for(i) in 1:network.popSize #adjusts location of each individual in between parent and child's old location
        if((network.popLocations[i] > network.popLocations[parent] && network.popLocations[i] < network.popLocations[parent] + deathDist) || (network.popLocations[i] < network.popLocations[parent] && network.popLocations[i] > network.popLocations[parent] + deathDist))
            network.popLocations[i] += deathDist/abs(deathDist) #adds or subtracts one based on sign of deathDist
        end
    end

    network.popLocations[child] = network.popLocations[parent] + deathDist/abs(deathDist) #child's location set to parent +/- 1
end

function resolveFitnesses(network::NetworkParameters)
    network.popFitness.=(1.0 + network.delta).^network.popPayoff
    network.popPayoff[:] .= 0.0
end

function neighborMom(network::NetworkParameters, kID::Int64, neighborRange::Int64=1)
    indexes = zeros(Int64, 2*neighborRange)
    fitnesses = zeros(2*neighborRange)
    for(i) in 1:neighborRange
        indexes[i] = (kID-i+99)%100+1 #index minus 1 through neighborRange, loops around 100
        indexes[i+neighborRange] = (kID+i+99)%100+1 #index plus 1 through neighborRange, loops around 100
        fitnesses[i] = network.popFitness[indexes[i]] #same for fitnesses
        if(indexes[i]==indexes[i+neighborRange]) #necessary for when neighborRange = popSize/2
            break
        end
        fitnesses[i+neighborRange] = network.popFitness[indexes[i+neighborRange]]
    end
    fitWeights = StatsBase.weights(fitnesses)
    momIndex = sample(indexes, fitWeights) 
    resolveLocs(network, kID, momIndex) #updates locations, ensures each individual has a unique location
    momIndex
end

function anyMom(network::NetworkParameters, kID::Int64, neighborRange::Int64=50)
    fitWeights = StatsBase.weights(network.popFitness)
    momIndex = sample(1:network.popSize, fitWeights)
    resolveLocs(network, kID, momIndex) #updates locations, ensures each individual has a unique location
    momIndex
end

function runSimsReturn(;B::Float64=2.0, C::Float64=0.5, D::Float64=0.0, CL::Float64=0.0, gen::Int=500, distInherit::Bool=false, distFactor::Float64=0.975, pn::Float64=0.5, pnd::Bool=false, pr::Float64=0.01, prd::Bool=false, muP::Float64=0.001, delta::Float64=0.1, sigmapn::Float64=0.05, sigmapr::Float64=0.01, reps::Int64=50)

        network = NetworkParameters(B, C, D, CL, gen, distInherit, distFactor, pn, pnd, pr, prd, muP, delta, sigmapn, sigmapr) #if pnd/prd = true, then defector & cooperator probabilities will evolve separately for pn/pr
        #network2 = NetworkParameters(B, C, D, CL, gen, distInherit, distFactor, pn, pnd, pr, prd, muP, delta, sigmapn, sigmapr) #if pnd/prd = true, then defector & cooperator probabilities will evolve separately for pn/pr
        #network2.edgeMatrix = network.edgeMatrix
        #cooperate(network2)
        #resolveFitnesses(network2)
        kid = 100
        #network.popFitness[100] = 6
        #network.popFitness[1] = 10
        #network.popFitness[5] = 5
        #network.popFitness[4] = 4
    
        network.popFitness[kid] = 0
        selections = zeros(100)
        selectionstwo = zeros(100)
        for i in 1:1000000
            parent = neighborMom(network, kid, 50)
            selections[parent]+=1
        end
        for i in 1:1000000
            parent = anyMom(network, kid, 3)
            selectionstwo[parent]+=1
        end
        print(selections-selectionstwo, " mean: ", mean((selections[2:100].-selectionstwo[2:100])))
        #print(selectionstwo)
        

        #oldcoop(network)
        #rfold(network)
        #println(network.popFitness)
    end

runSimsReturn()

