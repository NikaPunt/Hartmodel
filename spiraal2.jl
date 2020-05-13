# Code voor cellulaire automaat
using LightGraphs
using MetaGraphs
using DelimitedFiles
using DataStructures
include("plot_graph.jl")

mutable struct CellulaireAutomaat
    mg::MetaGraph
    #properties nodes:
        #1: non-active and can be depolarised -> DI (black)
        #2: active and fires in begin period -> APD (red)
        #3: non-active and can't be depolarised -> DI (green)
    time::Int64#Time since running time in t.u.
    δt::Int64 #the time unit in ms -> 1 t.u. = δt ms

    δx::Float64 #the space unit in cm -> 1 s.u. = δx cm

    ARI_ss_epi::Float64
    ARI_ss_endo::Float64

    a_epi::Float64
    a_endo::Float64

    b_epi::Float64
    b_endo::Float64

    #transitions of the active edges
    edgesA::PriorityQueue{Tuple{Int64,Int64}, Float64}

    #exit nodes for purkinje
    LVexit::Int64
    RVexit::Int64
    RVexit_time::Int64

    #steady state CV
    CV_ss
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
    #purkinje
    #add_edge!(mg,1297,13124)
    #set_prop!(mg,1297,13124,:anisotropy,1)
    #set_prop!(mg,1297,13124,:dx,0.01)
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
    nodecolor = [colorant"black", colorant"firebrick",colorant"blue"]
    #,colorant"darkred",colorant"firebrick4",colorant"firebrick",colorant"orangered3",colorant"orangered",colorant"darkorange1",colorant"orange",colorant"darkgoldenrod2",colorant"darkgoldenrod1",colorant"goldenrod1",colorant"gold",colorant"yellow2", colorant"yellow"]
    #nv = number of vertices
    membership = ones(Int64, nv(celAutom.mg))
    for i in range(1,stop=nv(celAutom.mg))
        membership[i]=get_prop(celAutom.mg,i,:state)
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
#   @post   The edges will be colored in a spectrum from white (not excited at all)
#           to dark red (completely excited),
#           depending on the ltransition and/or htransition in the edge.
#           Edges recovering from excitation will be colored green.
function coloringEdge(celAutom::CellulaireAutomaat)
    edgecolor = [colorant"black", colorant"red"]
    #ne = number of edges
    membership=ones(Int64,ne(celAutom.mg))
    j=1
    for edge in collect(edges(celAutom.mg))
        try
            celAutom.edgesA[tuple(edge.src,edge.dst)]
            membership[j]=2
        catch
        end
        try
            celAutom.edgesA[tuple(edge.dst,edge.src)]
            membership[j]=2
        catch
        end
        j+=1
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
#   @param  (Int64) dim
#           the dimension of the graph. Can be either 2 or 3.
#   @post   Generates an image of the colored MetaGraph object using the gplot
#           function and stored as $folder/frame$i.png.
#   @post   If dim = 2 this will generate a plot of a 2D graph.
#   @post   If dim = 3 this will generate a plot of a 3D graph.
function plotGraph2(celAutom::CellulaireAutomaat,i::Int64,folder::String,y_start::Float64,y_stop::Float64)
    nodes=get_area(celAutom.mg,-1000000.0,1000000.0,y_start,y_stop,-1000000.0,10000000.0)
    nodefillc = coloringGraph(celAutom)[nodes]
    I=[]
    j=1
    for edge in collect(edges(celAutom.mg))
        if (Int64(edge.src) in nodes) && (Int64(edge.dst) in nodes)
            append!(I,j)
        end
        j=j+1;
    end
    edgefillc =coloringEdge(celAutom)[I]
    loc_x,loc_y,loc_z = get_coordinates(celAutom.mg)
    print(nv(celAutom.mg))
    print(size(loc_z))
    g1=gplot(celAutom.mg,loc_x[nodes],loc_z[nodes],nodefillc=nodefillc,edgestrokec=edgefillc)
    draw(PNG("$folder/frame$i.png", 16cm, 16cm), g1)
end

## TODO specificatie saai
function plotGraph2(celAutom::CellulaireAutomaat,i::Int64,folder::String,y_start::Float64,y_stop::Float64)
    nodes=get_area(celAutom.mg,-1000000.0,1000000.0,y_start,y_stop,-1000000.0,10000000.0)
    nodefillc = coloringGraph(celAutom)[nodes]
    edges=[]
    j=1
    for edge in collect(edges(celAutom.mg))
        if (edge.src in nodes) && (edge.dst in nodes)
            append!(edges,j)
        end
        j=j+1;
    end
    edgefillc =coloringEdge(celAutom)[edges]
    loc_x,loc_y,loc_z = get_coordinates(celAutom.mg)[nodes,nodes,nodes]
    g1=gplot(celAutom.mg,loc_x,loc_z,nodefillc=nodefillc,edgestrokec=edgefillc)
    draw(PNG("$folder/frame$i.png", 16cm, 16cm), g1)
end
##
#   canTransitionInto will return false if and only if conduction can not travel
#   to the given node. This is a simple check of the node's state.
#TODO add time fraction so that when node goes from state 3 to 1 in that time,
# it can be activated -> @param
function canPassCurrentTo(celAutom::CellulaireAutomaat, node::Int64,timeFraction::Float64)
    if get_prop(celAutom.mg,node,:state)==1
        return true
    end
timeState3=get_prop(celAutom.mg,node,:APD)/4/celAutom.δt
    if get_prop(celAutom.mg,node,:state)==3
        return (get_prop(celAutom.mg,node,:tcounter)+timeFraction>=timeState3)
    end
    timeState2and3= timeState3+get_prop(celAutom.mg,node,:APD)
    return (get_prop(celAutom.mg,node,:tcounter)+timeFraction>=timeState2and3)
end
##
#   makeTransition will add the CV*anistropy/dx to every edge with a running
#   current. (No attention given to surpluses yet.)
#
#   @param  (CellulaireAutomaat) celAutom
#           An object of the class CellulaireAutomaat that carries an object of
#           type MetaGraph and a priority queue with all the edges that must undergo
#           the transition
#   @post   The conduction fraction in the edge will be updated with
#           old_conduction_fraction + CV*ani/dx
function makeTransition!(celAutom::CellulaireAutomaat)
    for key in keys(celAutom.edgesA)
        if (get_prop(celAutom.mg, key[1],:celtype)==2) && (get_prop(celAutom.mg, key[2], :celtype)==2)
            CV = celAutom.CV_ss
        else
            CV = get_prop(celAutom.mg, key[1], :CV)
        end
        ani = get_prop(celAutom.mg, key[1], key[2], :anisotropy)
        dx = get_prop(celAutom.mg, key[1],key[2],:dx)
        celAutom.edgesA[key]+= CV*ani/dx
    end
end

##
#   activate sets the node on state 2 and the tcounter on 0, recalculates CV
#   APD
#
#   @param  (CellulaireAutomaat) celAutom
#           An object of the class CellulaireAutomaat that carries an object of
#           type MetaGraph and a priority queue with all the edges that must undergo
#           the transition
#   @param  (Int64) node
#           The (number of the) node which has to be activated
#   @post   Will set the property state of node on 2 and the property tcounter of
#           node on 0
function activate!(celAutom::CellulaireAutomaat, node::Int64, timeFraction::Float64)
    set_prop!(celAutom.mg, node, :state, 2)
    set_prop!(celAutom.mg, node, :tcounter, get_prop(celAutom.mg, node, :tcounter)+timeFraction)
    set_APD!(celAutom,node)
    set_CV!(celAutom,node)
    set_prop!(celAutom.mg, node, :tcounter, -timeFraction)
end
##
#   canActivateEdge returns true if and only if there is no running current from
#   there to here. This means that a current can go from here to there.
#
#   @param  (CellulaireAutomaat) celAutom
#   @param  (Int64) here
#           The node from which the edge activation should start from
#   @param  (Int64) there
#           The other side of the edge
function canActivateEdge(celAutom::CellulaireAutomaat, here::Int64, there::Int64)
    try
        celAutom.edgesA[tuple(there, here)]
        return false
    catch
        return true
    end

end
##
#   handleSurplus is where we use the PriorityQueue embedded in the mutable struct
#   to handle the surplus in any edges that have a conduction fraction above 1.
#
#   @param  (CellulaireAutomaat) celAutom
#           An object of the class CellulaireAutomaat that carries an object of
#           type MetaGraph and a priority queue with all the edges that must undergo
#           the transition
#   @post   for all edges with a transition more than 1, we take that surplus and
#           add it to the neighbouring edges
#   TODO
#   2.  De fractie van de tijdstap berekenen om toe te voegen aan tcounter voor
#       de APD met de juiste edge

function handleSurplus(celAutom::CellulaireAutomaat)
    #als grootste transition groter is dan 1
    while (!isempty(celAutom.edgesA))&&peek(celAutom.edgesA)[2] > 1
        surplus=peek(celAutom.edgesA)[2]-1#fraction surplus with anistropy
        key=dequeue!(celAutom.edgesA)
        #properties of current edge
        dx=get_prop(celAutom.mg, key[1],key[2],:dx)
        ani=get_prop(celAutom.mg,key[1],key[2],:anisotropy)
        #If its purkinje channel
        if (get_prop(celAutom.mg, key[1],:celtype)==2) && (get_prop(celAutom.mg, key[2], :celtype)==2)
            CV = celAutom.CV_ss
        #if not
        else
            CV = get_prop(celAutom.mg, key[1], :CV)
        end
        #fraction of the time step to get to the node
        timeFraction=1-(surplus*dx)/(CV*ani)
        if canPassCurrentTo(celAutom, key[2],timeFraction)
            #activate node
            activate!(celAutom, key[2],timeFraction)
            #go to next edges and set transition
            for node in collect(setdiff(neighbors(celAutom.mg,key[2]),key[1]))
                if canActivateEdge(celAutom,key[2],node)
                    currentTransition=0
                    try
                        currentTransition+=celAutom.edgesA[tuple(key[2],node)]
                    catch
                    end
                    #get different properties of new edge
                    dx2=get_prop(celAutom.mg, key[2],node,:dx)
                    ani2=get_prop(celAutom.mg,key[2],node,:anisotropy)
                    #If its purkinje channel
                    if (get_prop(celAutom.mg, key[2],:celtype)==2) && (get_prop(celAutom.mg, node, :celtype)==2)
                        CV2 = celAutom.CV_ss
                    #if not
                    else
                        CV2 = get_prop(celAutom.mg, key[2], :CV)
                    end
                    #calculate the new transition
                    newTransition = surplus*dx/dx2*ani2/ani*CV2/CV
                    #set the transition on the max of the current and new transition
                    try
                        celAutom.edgesA[tuple(key[2],node)]=max(newTransition,currentTransition)
                    catch
                        enqueue!(celAutom.edgesA,(key[2], node),max(newTransition,currentTransition))
                    end
                else
                    try
                        delete!(celAutom.edgesA,(node,key[2]))
                    catch
                    end
                end
            end
        end
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
#   @post   If a node in state 1 has neighbouring edge(s) that have a ltransition
#           or htransition high enough to excite the node, the state of this node
#           will change from '1' to '2'.
function state1to2!(celAutom::CellulaireAutomaat)
    #bij alles cv bij (met anistropie) in uw priority
    makeTransition!(celAutom)
    #handleSurplus
    handleSurplus(celAutom)
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
#   @post   If a node has been in state 2 for a time tcounter longer than its APD,
#           its state will change from '2' to '3'. tcounter will be reset to
#           get_prop(celAutom.mg,node,:tcounter)-get_prop(celAutom.mg,node,:APD)

function state2to3!(celAutom::CellulaireAutomaat)
    for node in collect(filter_vertices(celAutom.mg,:state,2))
        if get_prop(celAutom.mg,node,:tcounter)+1>=get_prop(celAutom.mg,node,:APD)
            set_prop!(celAutom.mg,node,:state,3)
            set_prop!(celAutom.mg,node,:tcounter,get_prop(celAutom.mg,node,:tcounter)-get_prop(celAutom.mg,node,:APD))
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
#   @post   If a node in state 3 has a tcounter larger than timeState3,
#           then we set the state to 1
function state3to1!(celAutom::CellulaireAutomaat)
    for node in collect(filter_vertices(celAutom.mg,:state,3))
        timeState3=get_prop(celAutom.mg,node,:APD)/4/celAutom.δt
        if get_prop(celAutom.mg,node,:tcounter)+1>=timeState3
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
#   @param (Array{Int64, 1}) stopwaarden
#           The values of the nodes that can't be depolarised in te inital condition
#
#   The next six parameters are parameters of the formula APD = ARI_ss - a*exp(-DI*celAutom.δt/b)
#   @param  (Float64) ARI_ss_epi
#           The steady state ARI of the epicardium
#   @param  (Float64) ARI_ss_endo
#           The steady state ARI of the endocardium
#   @param  (Float64) a_epi
#           The a coefficient in the formula for the epicardium
#   @param  (Float64) a_endo
#           The a coefficient in the formula for the endocardium
#   @param  (Float64) b_epi
#           The b coefficient in the formula for the epicardium
#   @param  (Float64) b_endo
#           The b coefficient in the formula for the endocardium
#   @param  (Int64) RVexit_time
#           The time (in time units) when the node RVexit will get artifically
#           excited.
#   @pre    stopwaarden and startwaarden don't have common values
#   @post   The embedded MetaGraph has the property state imposed on the nodes.
#           The state of all nodes except the values in startwaarden are set to
#           '1' by default. Whilst the values in startwaarden are set to '2.'
#   @post   The embedded MetaGraph has the property tcounter imposed on the nodes.
#           (tcounter - time counter) this is set to one hundred in the start
#           to put the model in its steady state.
#   @post   The embedded MetaGraph has the property CV imposed on the nodes.
#           (CV - conduction velocity) The default setting of this property is TODO
#   @post   The embedded MetaGraph has the property APD imposed on the nodes.
#           (APD - action potential duration) The default setting of this property is TODO
#   @post   The embedded MetaGraph has the property ltransition imposed on the edges.
#           (ltransition - lower transition) This is the fraction of the
#           conduction in the edge going from the lower numbered node to
#           the higher numbered node.
#   @post   The embedded MetaGraph has the property htransition imposed on the edges.
#           (ltransition - higher transition) This is the fraction of the
#           conduction in the edge going from the higher numbered node to
#           the lower numbered node.
function createCellulaireAutomaat(graph::MetaGraph, dt::Int64, δx::Float64, startwaarden::Array{Any,1},
                                    stopwaarden::Array{Any,1}, ARI_ss_endo::Float64,
                                    ARI_ss_epi::Float64, a_epi::Float64, a_endo::Float64,
                                    b_epi::Float64, b_endo::Float64,
                                    LVexit::Int64, RVexit::Int64, RVexit_time::Int64,
                                    CV_ss::Float64)
    Priority = PriorityQueue{Tuple{Int64,Int64}, Float64}(Base.Order.Reverse)
    for node in startwaarden
        for buur in collect(neighbors(graph, node))
            #anders problemen met edges
            if !(buur in startwaarden) && !(buur == LVexit) &&!(buur in stopwaarden) &&!(buur==RVexit)
                enqueue!(Priority,(node, buur),0)
            end
        end
    end
#    for buur in collect(neighbors(graph, LVexit))
        #anders problemen met edges
#        if !(buur in startwaarden) && !(buur == LVexit) &&!(buur in stopwaarden) &&!(buur==RVexit)
#            enqueue!(Priority,(LVexit, buur),0)
#        end
#    end
    celAutom = CellulaireAutomaat((mg=graph,time=0,δt=dt,δx=δx, ARI_ss_endo=ARI_ss_endo,
                                ARI_ss_epi=ARI_ss_epi, a_epi=a_epi, a_endo=a_endo,
                                b_epi=b_epi, b_endo=b_endo, edgesA=Priority,
                                LVexit=LVexit, RVexit = RVexit,RVexit_time=RVexit_time,
                                CV_ss=CV_ss)...)
    for node in collect(vertices(graph))
        set_prop!(graph,node,:state,1)
        #CL negative but will invoke the methods so that their minimum values will be induced
        set_prop!(graph, node,:APD,0)
        set_prop!(graph,node,:tcounter,-1000)
        set_APD!(celAutom,node)
        set_CV!(celAutom,node)
    end

    for node in startwaarden
        set_prop!(graph,node,:state,2)
        set_prop!(graph,node,:tcounter,0)
    end
    set_prop!(graph,LVexit,:state,2)
    set_prop!(graph,LVexit,:tcounter,0)
    for node in stopwaarden
        set_prop!(graph,node,:state,2)
        set_prop!(graph,node,:tcounter,0)
    end
    return celAutom
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
#   @param  (Int64) dim
#           the dimension of the graph. Can be either 2 or 3.
#   @post   Generates an image of the colored MetaGraph object using the gplot
#           function and stored as $folder/frame$i.png.
#   @post   If dim = 2 this will generate a plot of a 2D graph.
#   @post   If dim = 3 this will generate a plot of a 3D graph.
function plotGraph2(celAutom::CellulaireAutomaat,i::Int64,folder::String, dim::Int64)
    if dim == 3
        loc_x,loc_y,loc_z = get_coordinates(celAutom.mg)
        nodefillc = coloringGraph(celAutom)
        #edgefillc =coloringEdge(celAutom)
        g2 = graphplot(celAutom.mg,dim=3, x = loc_x, y=loc_y, z=loc_z,
                        markercolor = nodefillc, linecolor = :black,
                        markersize = 4, linewidth=2, curves=false)
        p = plot(g2,show=true,grid=false,showaxis=false)
        png(p, "$folder/3Dframe$i.png")
    elseif dim == 2
        nodefillc = coloringGraph(celAutom)
        edgefillc =coloringEdge(celAutom)
        loc_x,loc_y,loc_z = get_coordinates(celAutom.mg)
        g1=gplot(celAutom.mg,loc_x,loc_y,nodefillc=nodefillc,edgestrokec=edgefillc)
        draw(PNG("$folder/frame$i.png", 16cm, 16cm), g1)
    end
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
#   @param  (Int64) RVexit_time
#           The time (in time units) when the node RVexit will get artifically
#           excited.
#   @post   #amount frames will be made. Each time by calling an updateState to
#           go to the next state and then calling plotGraph2 to plot the graph.
#   @post   Each frame will update celAutom.time by incrementing it by 1 time unit.
#           If celautom.time reaches RVexit_time, then the node RVexit will go in
#           state 2.
function createFrames(folderName::String, amountFrames::Int64, amountCalcs::Int64,
                        celAutom::CellulaireAutomaat, dim::Int64)
    #this calculates timesteps until we have the amount of necessary frames
    try
        mkdir(folderName)
    catch
    end
    plotGraph2(celAutom,0,folderName,dim)
    printIndex=floor(amountCalcs/amountFrames)
    for i in range(1,stop=amountCalcs)
        updateState!(celAutom)
        if mod(i,printIndex)==0
            plotGraph2(celAutom,Int64(i/printIndex),folderName, dim)
        end
        for node in collect(vertices(celAutom.mg))
            set_prop!(celAutom.mg,node,:tcounter,get_prop(celAutom.mg,node,:tcounter)+1)
        end
        #TODO title plot with celAutom.time*celAutom.δt
        celAutom.time+=1#1 time step further
        #If the time is right, we will excite RVexit.
        #if celAutom.time==celAutom.RVexit_time
    #        set_prop!(celAutom.mg,celAutom.RVexit,:state,2)
#            for buur in collect(neighbors(celAutom.mg, celAutom.RVexit))
#                #anders problemen met edges
#                if !(get_prop(celAutom.mg,buur,:state)==2)
#                    enqueue!(celAutom.edgesA,(celAutom.RVexit, buur),0)
#                end
#            end
#        end
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
    #CL = DI + APD (previous cycle) to ms
    CL = (get_prop(celAutom.mg, node, :tcounter) + get_prop(celAutom.mg, node, :APD))*celAutom.δt
    #minimum value of CL for APD, this only happens when invoked
    #under conditions which are not normal for a heart
    if CL<300
        CL=300
    end
    T = get_prop(celAutom.mg, node, :Temp)

    #epi APD
    ARI_epi = celAutom.ARI_ss_epi - celAutom.a_epi*exp(-CL/celAutom.b_epi)#Implementation formula in ms
    ARI_epi=ARI_epi/celAutom.δt#ms -> t.u

    #endo APD
    ARI_endo = celAutom.ARI_ss_endo - celAutom.a_endo*exp(-CL/celAutom.b_endo)#Implementation formula in ms
    ARI_endo=ARI_endo/celAutom.δt#ms -> t.u

    #Real APD
    ARI = ARI_endo*(1-T) + ARI_epi*T
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
    DI=get_prop(celAutom.mg, node, :tcounter)*celAutom.δt#DI in milisec
    #minimum value of DI for CV, this only happens when invoked
    #under conditions which are not normal for a heart
    if DI<40
        DI=40
    end
    CV=70.03-52.12*exp(-DI/87.6)#cm/sec
    CV=CV*celAutom.δt/1000#transition to cm/t.u.
    CV=CV/celAutom.δx#transition to s.u./t.u.
    set_prop!(celAutom.mg, node,:CV, CV)

end
##
function get_area(g::MetaGraph,x_start::Float64,x_stop::Float64,y_start::Float64,y_stop::Float64,z_start::Float64,z_stop::Float64)
    I = []
    for i in range(1,stop=nv(g))
         if x_start <= get_prop(g,i,:loc_x) <= x_stop
             if y_start <= get_prop(g,i,:loc_y) <= y_stop
                 if z_start <= get_prop(g,i,:loc_z) <= z_stop
                     append!(I, i)
                 end
             end
         end
    end
    return I
end
##
function main()
    start=time()
    graph = constructGraph("data_vertices.dat", "data_edges.dat", ',')

    startwaarden=get_area(graph, -100.0,115.0, 110.0,125.0,-1000.0,1000.0)
    stopwaarden = get_area(graph, -100.0,115.0,125.0,140.0,-1000.0,1000.0)


    #epi APD
    ARI_ss_epi = Float64(392.61)
    a_epi = Float64(339.39)
    b_epi = Float64(520)

    #endo APD
    ARI_ss_endo = Float64(505.74)
    a_endo = Float64(485.4)
    b_endo = Float64(501)

    #timestep in ms
    dt=20
    #spacestep in cm
    δx = 1/10

    #CV_ss
    CV_ss_mps = 2 #m/s
    CV_ss_sups = CV_ss_mps*100/δx #s.u./s
    CV_ss_suptu= CV_ss_sups/dt #s.u./t.u.

    LVexit = Int64(8427)
    RVexit = Int64(5837)
    time_ms = Int64(20) #na 20 miliseconden moet RVexit oplichten
    RVexit_time = Int64(ceil(time_ms/dt)) #aantal tijdstappen voordat RVexit moet oplichten.

    celAutom = createCellulaireAutomaat(graph, dt, δx, startwaarden,stopwaarden,
                        ARI_ss_endo, ARI_ss_epi, a_epi, a_endo, b_epi, b_endo,
                        LVexit, RVexit, RVexit_time, CV_ss_suptu)

    folder="plotjesxyDiepte"

    dim = 2
    createFrames(folder,100,100,celAutom, dim)
    elapsed = time() - start
    println(elapsed)
end

main()
