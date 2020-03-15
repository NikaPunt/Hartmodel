# Code voor cellulaire automaat
using LightGraphs
using MetaGraphs
using DelimitedFiles
include("plot_graph.jl")

mutable struct CellulaireAutomaat
    mg::MetaGraph
    #properties nodes:
        #1: non-active and can be depolarised -> DI (black)
        #2: active and fires in begin period -> APD (red)
        #3: non-active and can't be depolarised -> DI (green)
    time::Int64
    δt::Int64 #the time unit in ms
    δx::Float64 #the space unit in cm
end

##
#   constructGraph will construct a metagraph with several given properties
#   that are described in the header of the input files.
#
#   @param  (String) filename_vertices
#           This is the name of the file to be read containing the vertices
#   @param  (String) filename_edges
#           This is the name of the file to be read containing the edges
#   @param  (Char)  delimiter
#           This is the delimiter symbol that separates the streamed string
#           into parameter values.
#   @post   Returns a MetaGraph object with the properties described in
#           filename_edges and filename_vertices.
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
##
#   coloringGraph will return an array of colors that can color the nodes of
#   the MetaGraph embedded in the given CellulaireAutomaat. This array can be
#   used as the nodefillc array used in the gplot function.
#
#   @param  (CellulaireAutomaat) celAutom
#           An object of the class CellulaireAutomaat that carries an object of
#           type MetaGraph.
#   @post   Returns nodefillc, an array of type Array{RGB{Normed{UInt8, 8}},1}
#           with the needed coloring for the gplot function.
#   @post   nodefillc will assign the color "black" to nodes in state one,
#           "maroon" to nodes in state two and "red4" to nodes in state three.
function coloringGraph(celAutom::CellulaireAutomaat)
    nodecolor = [colorant"black", colorant"red",colorant"green"]
    #,colorant"darkred",colorant"firebrick4",colorant"firebrick",colorant"orangered3",colorant"orangered",colorant"darkorange1",colorant"orange",colorant"darkgoldenrod2",colorant"darkgoldenrod1",colorant"goldenrod1",colorant"gold",colorant"yellow2", colorant"yellow"]
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
##
#   coloringGraph will return an array of colors that can color the edges of
#   the MetaGraph embedded in the given CellulaireAutomaat. This array can be
#   used as the edgestrokec array used in the gplot function.
#
#   @param  (CellulaireAutomaat) celAutom
#           An object of the class CellulaireAutomaat that carries an object of
#           type MetaGraph.
#   @post   Returns colors, an array of type Array{RGB{Normed{UInt8, 8}},1}
#           with the needed coloring for the gplot function.
#   @post   The array colors will TODO
function coloringEdge(celAutom::CellulaireAutomaat)
    edgecolor = [colorant"white", colorant"maroon",colorant"red4",colorant"darkred",colorant"firebrick4",
    colorant"firebrick",colorant"green"]
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
            dx = 1 # get_prop(celAutom.mg,edge,:dx)
            membership[j]=ceil(Int64,(ltransition+htransition)/dx*5)+1
            j+=1
        end
    end
    for i in range(1,stop=ne(celAutom.mg))
        if membership[i] in [1 2 3 4 5 6]
        else
            println("Probleempje")
            membership[i]=7
        end
    end
    return colors=edgecolor[membership]
end
##
#   plotGraph2 uses the MetaGraph object embedded in the given
#   CellulaireAutomaat object (celAutom) to create an image of the (colored)
#   MetaGraph object.
#
#   @param  (CellulaireAutomaat) celAutom
#           An object of the class CellulaireAutomaat that carries an object of
#           type MetaGraph.
#   @param  (Int64) i
#           This is the number assigned to the end of the name of the
#           picture file.
#   @param  (String) folder
#           this is the name of the folder used to store the pictures. If it
#           exists it doesn't attempt to make a new one. If it doesn't exist,
#           the directory will be created.
#   @post   Generates an image of the colored MetaGraph object using the gplot
#           function and stored as $folder/frame$i.png.
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
##
#   This function can be called upon when ltransition or htransition equals to 0
#   and computes the transition of the wavefront in one time unit
#
#   @param  (CellulaireAutomaat) celAutom
#           An object of the class CellulaireAutomaat that carries an object of
#           type MetaGraph.
#   @param  (LightGraphs.SimpleGraphs.SimpleEdge{Int64}) edge
#           The edge over which the ltransition or htransition changes
#   @param  (Bool) direction
#           true when there is an ltransition
#           false when there is an htransition
#   @post   The conduction fraction in the edge will be updated.
function makeTransition!(celAutom::CellulaireAutomaat, edge::LightGraphs.SimpleGraphs.SimpleEdge{Int64},direction::Bool)
    if direction
        transitionSide = [:ltransition, :htransition]
        edgeSide = [edge.src, edge.dst]
    else
        transitionSide = [:htransition, :ltransition]
        edgeSide = [edge.dst, edge.src]
    end
    #Als ltransition nul is nemen we de htransition.
    transitionProp= get_prop(celAutom.mg,edge,transitionSide[2])
    #Conduction velocity CV
    CV = get_prop(celAutom.mg,edgeSide[2],:CV)
    #De dx in de file
    dx = get_prop(celAutom.mg,edge,Symbol(":dx"))
    #We look at the fraction of the dx in the next period and if it is lower then 100%
    if transitionProp + CV/dx < 1
        set_prop!(celAutom.mg,edge,transitionSide[2],transitionProp + CV/dx)
    else
        overschot = transitionProp*dx + CV - dx #overschot in s.u.
        for node in collect(setdiff(neighbors(celAutom.mg,edgeSide[1]),edgeSide[2]))
            if get_prop(celAutom.mg,node,:state)!=1
                overschot2 = 0
            else
                #Fraction in time in other edge = overschot/CV
                #We get the distance by multiplying with the CV of the other edge
                overschot2=overschot/CV*get_prop(celAutom.mg,edgeSide[1],:CV)
            end
            if node < edgeSide[1]
                set_prop!(celAutom.mg,node,edgeSide[1],:htransition,
                    max(overschot2/dx,get_prop(celAutom.mg,node,edgeSide[1],
                    :htransition)))
            else
                set_prop!(celAutom.mg,node,edgeSide[1],:ltransition,
                    max(overschot2/dx,get_prop(celAutom.mg,node,edgeSide[1],
                    :ltransition)))
            end
        end
        if get_prop(celAutom.mg, edgeSide[1], :state)==1
            set_prop!(celAutom.mg,edgeSide[1],:state,2)
            set_APD!(celAutom,edgeSide[1])
            set_CV!(celAutom,edgeSide[1])
            set_prop!(celAutom.mg,edgeSide[1],:tcounter,0)
        end
        set_prop!(celAutom.mg,edge,transitionSide[2],0)
    end
end

##
#   state1to2! looks at the MetaGraph object embedded in the given
#   CellulaireAutomaat object (celAutom) and evaluates every single node.
#   If a node meets all the requirements to switch from state one to two,
#   state1to2! will set the state property from those nodes to the value '2.'
#
#   @param  (CellulaireAutomaat) celAutom
#           An object of the class CellulaireAutomaat that carries an object of
#           type MetaGraph.
#   @post   If a node in state 1 TODO
function state1to2!(celAutom::CellulaireAutomaat)
    C = Set{Any}(collect(edges(celAutom.mg))) #alle edges
    D = Set{Any}(collect(filter_edges(celAutom.mg,:ltransition,0))) #edges met ltransition nul
    E = Set{Any}(collect(filter_edges(celAutom.mg,:htransition,0))) #edges met htransition nul

    #Itereer over alle edges waar er bij minstens één kant een transitie is.

    for edge in collect(setdiff(C, intersect(D, E)))
        #ltransition = kant van de edge met de lager genummerde vertex
        #htransition = kant van de edge met de hoger genummerde vertex
        if get_prop(celAutom.mg,edge,:ltransition) == 0
            makeTransition!(celAutom, edge, true)
        elseif get_prop(celAutom.mg,edge,:htransition) == 0
            makeTransition!(celAutom, edge, false)
        else
            set_prop!(celAutom.mg,edge,:ltransition,-1)
            set_prop!(celAutom.mg,edge,:htransition,-1)
        end
    end
end

##
#   state2to1! looks at the MetaGraph object embedded in the given
#   CellulaireAutomaat object (celAutom) and evaluates every single node.
#   If a node meets all the requirements to switch from state two to one,
#   state2to1! will set the state property from those nodes to the value '1.'
#
#   @param  (CellulaireAutomaat) celAutom
#           An object of the class CellulaireAutomaat that carries an object of
#           type MetaGraph.
#   @post   If a node in state 2 TODO
function state2to3!(celAutom::CellulaireAutomaat)
    for node in collect(filter_vertices(celAutom.mg,:state,2))
        if get_prop(celAutom.mg,node,:tcounter)>=get_prop(celAutom.mg,node,:APD)
            set_prop!(celAutom.mg,node,:state,3)
            set_prop!(celAutom.mg,node,:tcounter,0)
        end
    end
end
##
#   state3to1! looks at the MetaGraph object embedded in the given
#   CellulaireAutomaat object (celAutom) and evaluates every single node.
#   If a node meets all the requirements to switch from state three to one,
#   state3to1! will set the state property from those nodes to the value '1.'
#
#   @param  (CellulaireAutomaat) celAutom
#           An object of the class CellulaireAutomaat that carries an object of
#           type MetaGraph.
#   @post   If a node in state 3 has a tcounter larger than 18.5 ms,
#           then we set the state to 1
function state3to1!(celAutom::CellulaireAutomaat)
    for node in collect(filter_vertices(celAutom.mg,:state,3))
        #18.5 ms / δt time units
        if get_prop(celAutom.mg,node,:tcounter)>=(18.5/celAutom.δt)
            set_prop!(celAutom.mg,node,:state,1)
        end
    end
end
##
#   updateState! will update the states of all the nodes of the embedded
#   MetaGraph object in the CellulaireAutomaat object (celAutom).
#
#   @param  (CellulaireAutomaat) celAutom
#           An object of the class CellulaireAutomaat that carries an object of
#           type MetaGraph.
#   @effect This function calls state1to2! first, then state2to3! and finally
#           state3to1!.
function updateState!(celAutom::CellulaireAutomaat)
    state1to2!(celAutom)
    state2to3!(celAutom)
    state3to1!(celAutom)
end
##
#   createCellulaireAutomaat will create an instance of CellulaireAutomaat,
#   with the given graph, and node(s) with initial excitation.
#
#   @param  (MetaGraph) graph
#           The graph to be embedded in the new cellular automaton.
#   @param  (Array{Int64, 1}) startwaarden
#           The values of the nodes that have an initial excitation.
#   @post   The embedded MetaGraph has the property state imposed on the nodes.
#           The state of all nodes except the values in startwaarden are set to
#           '1' by default. Whilst the rest are set to '2.'
#   @post   The embedded MetaGraph has the property tcounter imposed on the nodes.
#           (tcounter - time counter) this is set to one hundred in the start
#           to put the model in its steady state.
#   @post   The embedded MetaGraph has the property CV imposed on the nodes.
#           (CV - conduction velocity) TODO
#   @post   The embedded MetaGraph has the property APD imposed on the nodes.
#           (APD - action potential duration) TODO
#   @post   The embedded MetaGraph has the property ltransition imposed on the edges.
#           (ltransition - lower transition) This is the fraction of the
#           conduction in the edge going from the lower numbered node to
#           the higher numbered node.
#   @post   The embedded MetaGraph has the property htransition imposed on the edges.
#           (ltransition - higher transition) This is the fraction of the
#           conduction in the edge going from the higher numbered node to
#           the lower numbered node.
function createCellulaireAutomaat(graph::MetaGraph, startwaarden::Array{Int64,1})
    celAutom = CellulaireAutomaat((mg=graph,time=0,δt=5,δx=1)...)#,tDisp = 100, tSamp = 100)...)
    for node in collect(vertices(graph))
        set_prop!(graph,node,:state,1)
        set_prop!(graph,node,:CV,0.09365)
        set_prop!(graph,node,:tcounter,100)
        set_prop!(graph,node,:APD,242/celAutom.δt)
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
##
#   createFrames will create frames as much as the given amount. This function
#   calls the updateState function to update everything to the next timestep.
#   Then it calls plotGraph2 to create the image.
#
#   @param  (String) folderName
#           The name of the folder to be made. If it
#           exists it doesn't attempt to make a new one. If it doesn't exist,
#           the directory will be created.
#   @param  (Int64) amountFrames
#           The amount of frames to make.
#   @param  (Int64) amountCalcs
#           The amount of calculations to be made
#   @param  (CellulaireAutomaat) celAutom
#           An object of the class CellulaireAutomaat that carries an object of
#           type MetaGraph.
#   @post   #amount frames will be made. Each time by calling an updateState to
#           go to the next state and then calling plotGraph2 to plot the graph.
function createFrames(folderName::String, amountFrames::Int64, amountCalcs::Int64,
                        celAutom::CellulaireAutomaat)
    #this calculates timesteps until we have the amount of necessary frames
    try
        mkdir(folderName)
    catch
    end
    plotGraph2(celAutom,0,folderName)
    printIndex=floor(amountCalcs/amountFrames)
    for i in range(1,stop=amountCalcs)
        updateState!(celAutom)
        if mod(i,printIndex)==0
            plotGraph2(celAutom,Int64(i/printIndex),folderName)
        end
        for node in collect(vertices(celAutom.mg))
            set_prop!(celAutom.mg,node,:tcounter,get_prop(celAutom.mg,node,:tcounter)+1)
        end
        celAutom.time+=celAutom.δt
    end
end
##
#   This function sets the APD of the given node. The used formula was taken
#   from a paper's function fit. We make the assumption
#   Absolute Restitution Interval = Action Potential Duration.
#
#   @param  (CellulaireAutomaat) celAutom
#           An object of the class CellulaireAutomaat that carries an object of
#           type MetaGraph.
#   @param  (Int64) node
#           The value of the node in the graph, whose APD needs to be changed.
#   @post   The given node will have its APD property changed to the new value
#           of: ARI = ARI_ss - a*exp(-DI*celAutom.δt/b)
function set_APD!(celAutom::CellulaireAutomaat, node::Int64)
    DI= get_prop(celAutom.mg, node, :tcounter)
    ARI_ss = 242
    a = 404
    b = 36
    ARI = ARI_ss - a*exp(-DI*celAutom.δt/b)#Implementation formula in ms
    ARI=ARI/celAutom.δt#ms -> t.u
    set_prop!(celAutom.mg,node,:APD, ARI)
end
##
#   This function sets the CV of the given node.
#
#   @param  (CellulaireAutomaat) celAutom
#           An object of the class CellulaireAutomaat that carries an object of
#           type MetaGraph.
#   @param  (Int64) node
#           The value of the node in the graph, whose APD needs to be changed.
#   @post   The given node will have its CV property changed to the new value
#           of: CV = 70.03-52.12*e^(-DI/87.6)
function set_CV!(celAutom::CellulaireAutomaat, node::Int64)
    DI=get_prop(celAutom.mg, node, :tcounter)/celAutom.δt*1000#DI in sec
    CV=70.03-52.12*exp(-DI/87.6)#cm/sec
    CV=CV*celAutom.δt/1000#transition to cm/t.u.
    CV=CV*celAutom.δx#transition to s.u./t.u.
    set_prop!(celAutom.mg, node,:CV, CV)
end
##
function main()
    start = time()
    graph = constructGraph("data_vertices_2D.dat", "data_edges_2D.dat", ',')
    startwaarden = [1, 50, 124, 245, 378, 472, 596, 632, 780, 875]
    celAutom = createCellulaireAutomaat(graph, startwaarden)
    folder="plotjes_test4"
    createFrames(folder,250,250,celAutom)
    elapsed = time() - start
    println(elapsed)
end

main()
