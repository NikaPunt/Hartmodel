# Code voor cellulaire automaat
using LightGraphs
using MetaGraphs
using DelimitedFiles
include("plot_graph.jl")

mutable struct CellulaireAutomaat
    mg::MetaGraph
    fireList::Set{Int64}
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
    copiedFireList = collect(celAutom.fireList)
    for i in range(1, size(copiedFireList, 1))
        membership[i] = 2
    end
    return nodefillc = nodecolor[membership]
end

function plotGraph2(celAutom::CellulaireAutomaat,i::Int64)
    nodefillc = coloringGraph(celAutom)
    loc_x,loc_y,loc_z = get_coordinates(celAutom.mg)
    g1=gplot(celAutom.mg,loc_x,loc_y,nodefillc=nodefillc)
    draw(PNG("frame$i.png", 16cm, 16cm), g1)
end

function updateFireList(celAutom::CellulaireAutomaat)
    vuurlijst= Set{Int64}([])
    for i in range(1,stop=length(celAutom.fireList))
        buren = Set{Int64}(neighbors(celAutom.mg, pop!(celAutom.fireList)))
        union(vuurlijst, buren)
    end
    celAutom.fireList = vuurlijst
end

function main()
    graph = constructGraph("data_vertices_2D.dat", "data_edges_2D.dat", ',')
    celAutom = CellulaireAutomaat((mg = graph, fireList = Set{Int64}([1, 50 ,60, 850]))...)
    plotGraph2(celAutom,0)
    for i in range(1,stop=15)
        updateFireList(celAutom)
        #print(i)
        plotGraph2(celAutom,i)
        #print(i)
    end
end

main()
