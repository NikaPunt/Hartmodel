# Code voor cellulaire automaat
using LightGraphs
using MetaGraphs
using DelimitedFiles
include("plot_graph.jl")

mutable struct CellulaireAutomaat
    mg::MetaGraph
    #properties nodes:
        #1: non-active and can be depolarised -> DI
        #2: active and can fire in next period -> APD
        #3: active and can't fire -> APD
        #4: non-active and can't be depolarised -> DI
    fireDict::Dict{Int64, Set{Int64}}    

end

function constructGraph(filename_vertices::String,filename_edges::String,delimiter::Char)
    input_vertices = readdlm(filename_vertices,delimiter)
    input_edges = readdlm(filename_edges,delimiter)
    # define metagraph
    mg = MetaGraph(SimpleGraph())
    # add vertices
    for i in range(2,stop=size(input_vertices,1))
        # add one vertex
        add_vertex!(mg)
        # add other labels
        for j in range(2,stop=size(input_vertices,2))
            set_prop!(mg,Int(input_vertices[i,1]),Symbol(input_vertices[1,j]), input_vertices[i,j])
        end
    end
    # add edges
    for i in range(2,stop=size(input_edges,1))
        add_edge!(mg,Int(input_edges[i,1]),Int(input_edges[i,2]))
        # add labels
        # println(input_edges[i,1]," ",input_edges[i,2])
        for j in range(3,stop=size(input_edges,2))
            set_prop!(mg,Int(input_edges[i,1]),Int(input_edges[i,2]),Symbol(input_edges[1,j]),input_edges[i,j])
        end
    end
    println("Graph check!")
    return mg
end

#returns a nodefillc
function coloringGraph(celAutom::CellulaireAutomaat)
    nodecolor = [colorant"black", colorant"maroon",colorant"red4",colorant"darkred",colorant"firebrick4",
    colorant"firebrick",colorant"orangered3",colorant"orangered",colorant"darkorange1",colorant"orange",
    colorant"darkgoldenrod2",colorant"darkgoldenrod1",colorant"goldenrod1",colorant"gold",colorant"yellow2",
    colorant"yellow"]
    #nv = number of vertices
    membership = ones(Int64, nv(celAutom.mg))
    ActiveList = collect(celAutom.fireDict[2,3])
    for i in range(1, size(ActiveList, 1))
        membership[i] = 2
    end
    return nodefillc = nodecolor[membership]
end

function plotGraph2(celAutom::CellulaireAutomaat,i::Int64)
    nodefillc = coloringGraph(celAutom)
    loc_x,loc_y,loc_z = get_coordinates(celAutom.mg)
    g1=gplot(celAutom.mg,loc_x,loc_y,nodefillc=nodefillc)
    draw(PNG("plotjes_test2/frame$i.png", 16cm, 16cm), g1)
end

function updateFireDict(celAutom::CellulaireAutomaat)
    for i in range(1, stop=length(celAutom.fireDict[4]))
        push!(celAutom.fireDict[1], pop!(celAutom.fireDict[4]))
    end
    for i in range(1, stop=length(celAutom.fireDict[3]))
        push!(celAutom.fireDict[4], pop!(celAutom.fireDict[3]))
    end
    for i in range(1, stop=length(celAutom.fireDict[2]))
        fireElement=pop!(celAutom.fireDict[2])
        push!(celAutom.fireDict[3],fireElement)
        for i in range(1,stop=length(neighbors(fireElement)))
            if neighbors(fireElement)[i] in celAutom.fireDict[1]
                push!(celAutom.fireDict[2],neighbors(fireElement))
            end
        end
    end
end

function main()
    graph = constructGraph("data_vertices_2D.dat", "data_edges_2D.dat", ',')
    fireDict = Dict{Int64, Set{Int64}}
    fireDict[1] = setdiff(Set{Int64}(collect(vertices(graph))),Set{Int64}([1 55 203]))
    fireDict[2] = Set{Int64}([1 55 203])
    fireDict[3] = Set{Int64}([])
    fireDict[4] = Set{Int64}([])
    celAutom = CellulaireAutomaat((mg = graph, fireDict = fireDict)...)
    plotGraph2(celAutom,0)
    for i in range(1,stop=15)
        updatefireDict(celAutom)
        #print(i)
        plotGraph2(celAutom,i)
        #print(i)
    end
end

main()
