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
#    ts::Int64 #time sample
#    td::Int64 #time display
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
    for node in filter_vertices(celAutom.mg,:state,2)
        membership[node]=2
    end
    for node in filter_vertices(celAutom.mg,:state,3)
        membership[node]=3
    end
    return nodefillc = nodecolor[membership]
end

function plotGraph2(celAutom::CellulaireAutomaat,i::Int64)
    nodefillc = coloringGraph(celAutom)
    loc_x,loc_y,loc_z = get_coordinates(celAutom.mg)
    g1=gplot(celAutom.mg,loc_x,loc_y,nodefillc=nodefillc)
    draw(PNG("plotjes_test3/frame$i.png", 16cm, 16cm), g1)
end

function updateFireDict(celAutom::CellulaireAutomaat)
    for node in filter_vertices(celAutom.mg,:state,4)
        set_prop!(celAutom.mg,node,:state,1)
    end
    for node in filter_vertices(celAutom.mg,:state,3)
        set_prop!(celAutom.mg,node,:state,4)
    end
    looparray=collect(filter_vertices(celAutom.mg,:state,2))
    for node in looparray
        set_prop!(celAutom.mg,node,:state,3)
        buren=neighbors(celAutom.mg,node)
        for buur in buren
            if get_prop(celAutom.mg,buur,:state)==1
                set_prop!(celAutom.mg,buur,:state,2)
            end
        end
    end
end

function main()
    graph = constructGraph("data_vertices_test2.dat", "data_edges_test2.dat", ',')
    celAutom = CellulaireAutomaat(graph)
    for node in collect(vertices(graph))
        set_prop!(graph,node,:state,1)
    end
    startwaarden = [1, 50 ,650, 780]
    for node in startwaarden
        set_prop!(graph,node,:state,2)
    end
    plotGraph2(celAutom,0)
    for i in range(1,stop=50)
        updateFireDict(celAutom)
        #print(i)
        plotGraph2(celAutom,i)
        #print(i)
    end
end

main()
