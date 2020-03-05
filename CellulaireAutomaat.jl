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
    # tDisp::Int64 #display time
    # tSamp::Int64 #sample time
    time::Int64
    δt::Int64
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
function coloringEdge(celAutom::CellulaireAutomaat)
    edgecolor = [colorant"white", colorant"maroon",colorant"red4",colorant"darkred",colorant"firebrick4",
    colorant"firebrick"]
    #ne = number of edges
    membership=ones(Int64,ne(celAutom.mg))
    j=1
    for edge in collect(edges(celAutom.mg))
    ltransition=get_prop(celAutom.mg,edge,:ltransition)
        if ltransition==-1
            membership[j]=6
            j+=1
        else
            htransition=get_prop(celAutom.mg,edge,:htransition)
            dx = 1 # get_prop(celAutom.mg,edge,:dx
            membership[j]=ceil(Int64,(ltransition+htransition)/dx*5)+1
            j+=1
        end
    end
    return colors=edgecolor[membership]
end

function plotGraph2(celAutom::CellulaireAutomaat,i::Int64,folder::String)
    nodefillc = coloringGraph(celAutom)
    edgefillc =coloringEdge(celAutom)
    loc_x,loc_y,loc_z = get_coordinates(celAutom.mg)
    g1=gplot(celAutom.mg,loc_x,loc_y,nodefillc=nodefillc,edgestrokec=edgefillc)
    draw(PNG("$folder/frame$i.png", 16cm, 16cm), g1)
    for edge in collect(filter_edges(celAutom.mg,:htransition,-1))
        set_prop!(celAutom.mg,edge,:ltransition,0)
        set_prop!(celAutom.mg,edge,:htransition,0)
    end
end


function state1to2(celAutom::CellulaireAutomaat)
    C = Set{Any}(collect(edges(celAutom.mg))) #alle edges
    D = Set{Any}(collect(filter_edges(celAutom.mg,:ltransition,0))) #edges met ltransition nul
    E = Set{Any}(collect(filter_edges(celAutom.mg,:htransition,0))) #edges met htransition nul

    #Itereer over alle edges waar er bij minstens één kant een transitie is.
    for edge in collect(setdiff(C, intersect(D, E)))
        #ltransition = kant van de edge met de lager genummerde vertex
        #htransition = kant van de edge met de hoger genummerde vertex
        if get_prop(celAutom.mg,edge,:ltransition) == 0
            #Als ltransition nul is nemen we de htransition.
            transitionProp= get_prop(celAutom.mg,edge,:htransition)
            #Conduction velocity CV
            CV = get_prop(celAutom.mg,edge.dst,:CV)
            #De dx in de file
            dx = get_prop(celAutom.mg,edge,Symbol(":dx"))
            if transitionProp + CV < dx
                set_prop!(celAutom.mg,edge,:htransition,transitionProp + CV)
            else
                overschot = transitionProp + CV - dx
                for node in collect(neighbors(celAutom.mg,edge.src))
                    if get_prop(celAutom.mg,node,:state)==2
                        overschot2 = 0
                    else
                        overschot2=overschot
                    end
                    if node < edge.src
                        set_prop!(celAutom.mg,node,edge.src,:htransition,
                            max(overschot2,get_prop(celAutom.mg,node,edge.src,
                            :htransition)))
                    else
                        set_prop!(celAutom.mg,node,edge.src,:ltransition,
                            max(overschot2,get_prop(celAutom.mg,node,edge.src,
                            :ltransition)))
                    end
                    set_prop!(celAutom.mg,edge.src,:state,2)
                    set_APD!(celAutom.mg,edge.src)
                    set_prop!(celAutom.mg,edge.src,:tcounter,0)
                end
                set_prop!(celAutom.mg,edge,:htransition,0)
            end
        elseif get_prop(celAutom.mg,edge,:htransition) == 0
            transitionProp= get_prop(celAutom.mg,edge,:ltransition)
            CV = get_prop(celAutom.mg,edge.src,:CV)
            dx = get_prop(celAutom.mg,edge,Symbol(":dx"))
            if transitionProp + CV < dx
                set_prop!(celAutom.mg,edge,:ltransition,transitionProp+ CV)
            else
                overschot = transitionProp+ CV - dx
                for node in collect(neighbors(celAutom.mg,edge.dst))
                    if get_prop(celAutom.mg,node,:state)==2
                        overschot2 = 0
                    else
                        overschot2=overschot
                    end
                    if node < edge.dst
                        set_prop!(celAutom.mg,node,edge.dst,:htransition,
                            max(overschot2,get_prop(celAutom.mg,node,edge.dst,
                            :htransition)))
                    else
                        set_prop!(celAutom.mg,node,edge.dst,:ltransition,
                            max(overschot2,get_prop(celAutom.mg,node,edge.dst,
                            :ltransition)))
                    end
                    set_prop!(celAutom.mg,edge.dst,:state,2)
                    set_APD!(celAutom.mg, edge.dst)
                    set_prop!(celAutom.mg,edge.dst,:tcounter,0)
                end
                set_prop!(celAutom.mg,edge,:ltransition,0)
            end
        else
            set_prop!(celAutom.mg,edge,:ltransition,-1)
            set_prop!(celAutom.mg,edge,:htransition,-1)
        end
    end
end

function state2to1(celAutom::CellulaireAutomaat)
    for node in collect(filter_vertices(celAutom.mg,:state,2))
        if get_prop(celAutom.mg,node,:tcounter)>=get_prop(celAutom.mg,node,:APD)
            set_prop!(celAutom.mg,node,:state,1)
            set_prop!(celAutom.mg,node,:tcounter,0)
        end
    end
end

function updateState(celAutom::CellulaireAutomaat)
    state1to2(celAutom)
    state2to1(celAutom)

    # for node in filter_vertices(celAutom.mg,:state,4)
    #     set_prop!(celAutom.mg,node,:state,1)
    # end
    # for node in filter_vertices(celAutom.mg,:state,3)
    #     set_prop!(celAutom.mg,node,:state,4)
    # end
    # looparray=collect(filter_vertices(celAutom.mg,:state,2))
    # for node in looparray
    #     set_prop!(celAutom.mg,node,:state,3)
    #     buren=neighbors(celAutom.mg,node)
    #     for buur in buren
    #         if get_prop(celAutom.mg,buur,:state)==1
    #             set_prop!(celAutom.mg,buur,:state,2)
    #         end
    #     end
    # end
end

function createCellulaireAutomaat(graph::MetaGraph, startwaarden::Array{Int64,1})
    celAutom = CellulaireAutomaat((mg=graph,time=0)...)#,tDisp = 100, tSamp = 100)...)
    for node in collect(vertices(graph))
        set_prop!(graph,node,:state,1)
        set_prop!(graph,node,:CV,0.9365)
        set_prop!(graph,node,:tcounter,0)
    end
    for edge in collect(edges(graph))
        set_prop!(graph,edge,:ltransition,0)
        set_prop!(graph,edge,:htransition,0)
    end
    for node in startwaarden
        set_prop!(graph,node,:state,2)
        CV = get_prop(graph,node,:CV)
        for buur in collect(neighbors(graph, node))
            if buur < node
                set_prop!(graph,buur,node,:htransition,CV)
            else
                set_prop!(graph,node,buur,:ltransition,CV)
            end
        end
    end
    return celAutom
end

#Creates the frames inside a folder, creates the folder if non-existent
function createFrames(folderName::String,amount::Int64,celAutom::CellulaireAutomaat)
    #this calculates timesteps until we have the amount of necessary frames
    #mkdir(folderName)
    plotGraph2(celAutom,0,folderName)
    for i in range(1,stop=amount)
        updateState(celAutom)
        plotGraph2(celAutom,i,folderName)
        for node in collect(vertices(celAutom.mg))
            set_prop!(celAutom.mg,node,:tcounter,get_prop(celAutom.mg,node,:tcounter)+1)
        end
    end
end

#We make the assumption ARI=APD
function set_APD!(celAutom::CellulaireAutomaat, node)
    #Implementation formula in ms
    DI= get_prop(celAutom.mg, node, :tcounter)
    ARI_ss = 242
    a = 404
    b = 36
    ARI = ARI_ss - a*exp(-DI/b)
    #ms -> t.u.
    ARI=ARI/celAutom.δt
    set_prop!(celAutom.mg,node,:APD, ARI)
end

function main()
    start = time()
    graph = constructGraph("data_vertices_2D.dat", "data_edges_2D.dat", ',')
    startwaarden = [1, 50, 124, 245, 378, 472, 596, 632, 780, 875]
    celAutom = createCellulaireAutomaat(graph, startwaarden)
    folder="plotjes_test4"
    createFrames(folder,50,celAutom)
    elapsed = time() - start
    println(elapsed)
end

main()
