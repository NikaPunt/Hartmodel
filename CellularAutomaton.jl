# Code for cellular automaton
using DelimitedFiles
using DataStructures
include("plot_graph.jl")
include("CellularAutomatonMutableStruct.jl")

##
#   createFrames will create frames as much as the given amount. This function
#   calls the calculateTimeStep! function to update everything to the next timestep.
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
#   @param  (CellularAutomaton) celAutom
#           An object of the class CellularAutomaton that carries an object of
#           type MetaGraph.
#   @post   amountCals frames will be made. For this we call plotGraph2.
#   @post   The function calculateTimeStep! will be called amountCalcs times.
function createFrames(folderName::String, amountFrames::Int64, amountCalcs::Int64,
                        celAutom::CellularAutomaton)
    #this calculates timesteps until we have the amount of necessary frames
    try
        mkdir(folderName)
    catch
    end
    plotGraph2(celAutom,0,folderName)
    printIndex=floor(amountCalcs/amountFrames)
    for i in range(1,stop=amountCalcs)
        calculateTimeStep!(celAutom)
        if mod(i,printIndex)==0
            plotGraph2(celAutom,Int64(i/printIndex),folderName)
        end
    end
end
##
#   calculateTimeStep! will update the states of all the nodes of the embedded
#   MetaGraph object in the CellularAutomaton object (celAutom).
#
#   @param  (CellularAutomaton) celAutom
#           An object of the class CellularAutomaton that carries an object of
#           type MetaGraph.
#   @effect This function calls state1to2! first, then state2to3! and
#           state3to1! and finally increaseTimeStep!().
function calculateTimeStep!(celAutom::CellularAutomaton)
    state1to2!(celAutom)
    state2to3!(celAutom)
    state3to1!(celAutom)
    increaseTimeStep!(celAutom)
end
##
#   state1to2! looks at the MetaGraph object embedded in the given
#   CellularAutomaton object (celAutom) and updates the active edges
#   for the next time step.
#
#   @param  (CellularAutomaton) celAutom
#           An object of the class CellularAutomaton that carries an object of
#           type MetaGraph.
#   @effect All the active edges will increase in priority with the fraction
#           they will have traveled in the next time step
#           |makeTransition!(celAutom)
#   @effect All priorities af the active edges that are greater then 1 will be
#           adjusted properly
#           |handleSurplus(celAutom)
function state1to2!(celAutom::CellularAutomaton)
    makeTransition!(celAutom)
    handleSurplus(celAutom)
end
##
#   state2to3! looks at the MetaGraph object embedded in the given
#   CellularAutomaton object (celAutom) and evaluates every single node.
#   If a node meets all the requirements to switch from state two to three,
#   state2to3! will set the state property from those nodes to the value '3'.
#
#   @param  (CellularAutomaton) celAutom
#           An object of the class CellularAutomaton that carries an object of
#           type MetaGraph.
#   @post   If the node's tcounter has surpassed its APD, its :state property
#           will change from 2 to 3. tcounter will be reset to
#           get_prop(celAutom.mg,node,:tcounter)-get_prop(celAutom.mg,node,:APD)
function state2to3!(celAutom::CellularAutomaton)
    for node in collect(filter_vertices(celAutom.mg,:state,2))
        if get_prop(celAutom.mg,node,:tcounter)+1>=get_prop(celAutom.mg,node,:APD)
            set_prop!(celAutom.mg,node,:state,3)
            set_prop!(celAutom.mg,node,:tcounter,get_prop(celAutom.mg,node,:tcounter)-get_prop(celAutom.mg,node,:APD))
        end
    end
end
##
#   state3to1! looks at the MetaGraph object embedded in the given
#   CellularAutomaton object (celAutom) and evaluates every single node.
#   If a node meets all the requirements to switch from state three to one,
#   state3to1! will set the state property from those nodes to the value '1.'
#
#   @param  (CellularAutomaton) celAutom
#           An object of the class CellularAutomaton that carries an object of
#           type MetaGraph.
#   @post   If a node in state 3 has a tcounter larger than get_timeState3(),
#           then we set the state to 1
function state3to1!(celAutom::CellularAutomaton)
    for node in collect(filter_vertices(celAutom.mg,:state,3))
        timeState3=get_timestate3(celAutom,node)
        if get_prop(celAutom.mg,node,:tcounter)+1>=timeState3
            set_prop!(celAutom.mg,node,:state,1)
        end
    end
end
##
#   increaseTimeStep! will increase the tcounter of each node and the overall program
#   runtime. It also does a check to see if it must excite RVexit.
#
#   @param  (CellularAutomaton) celAutom
#           An object of the class CellularAutomaton that carries an object of
#           type MetaGraph.
#   @post   Each node will have its tcounter incremented by 1.
#   @post   celAutom.time will be incremented by 1.
#   @post   if celAutom.time reaches celAutom.RVexit_time and celAutom.purkinje_flag
#           is set to true, then we call exciteRVexit!()
function increaseTimeStep!(celAutom::CellularAutomaton)
    for node in collect(vertices(celAutom.mg))
        set_prop!(celAutom.mg,node,:tcounter,get_prop(celAutom.mg,node,:tcounter)+1)
    end
    celAutom.time+=1#1 time step further
    #If the time is right, we will excite RVexit.
    if celAutom.time==celAutom.RVexit_time && celAutom.purkinje_flag == true
        try
            exciteRVexit!(celAutom)
        catch
            println("!!! ATTENTION !!!")
            println("Current reaches RVexit before excitation can happen!")
            println("Cannot depolarise RVexit during its APD!")
        end
    end
end
##
#   exciteRVexit! will excite the node RVexit.
#
#   @param  (CellularAutomaton) celAutom
#           An object of the class CellularAutomaton that carries an object of
#           type MetaGraph.
#   @post   The :state property of RVexit is set to 2.
#   @post   For each excitable neighbor (state 1) of RVexit we will add the key
#           tuple(celAutom.RVexit, neighbor) to celAutom.edgesA with value 0.
function exciteRVexit!(celAutom::CellularAutomaton)
    set_prop!(celAutom.mg,celAutom.RVexit,:state,2)
    for buur in collect(neighbors(celAutom.mg, celAutom.RVexit))
        #anders problemen met edges
        if (get_prop(celAutom.mg,buur,:state)==1)
            enqueue!(celAutom.edgesA,(celAutom.RVexit, buur),0)
        end
    end
end
##
#   makeTransition will add the CV*anistropy/dx to every edge with a running
#   current. (No attention given to surpluses yet.)
#
#   @param  (CellularAutomaton) celAutom
#           An object of the class CellularAutomaton that carries an object of
#           type MetaGraph and a priority queue with all the edges that must undergo
#           the transition
#   @post   The conduction fraction in the edge will be updated with
#           old_conduction_fraction + CV*ani/dx,
#           the CV will be the property of the startnode or
#           be celAutom.CV_purkinje when we are traveling a purkinje edge
function makeTransition!(celAutom::CellularAutomaton)
    for key in keys(celAutom.edgesA)
        if (get_prop(celAutom.mg, key[1],:celtype)==2) && (get_prop(celAutom.mg, key[2], :celtype)==2)
            CV = celAutom.purkinje_CV_multiplier*get_prop(celAutom.mg, key[1], :CV)
        else
            CV = get_prop(celAutom.mg, key[1], :CV)
        end
        ani = get_prop(celAutom.mg, key[1], key[2], :anisotropy)
        dx = get_prop(celAutom.mg, key[1],key[2],:dx)
        celAutom.edgesA[key]+= CV*ani/dx
    end
end
##
#   handleSurplus is where we use the PriorityQueue embedded in the mutable struct
#   to handle the surplus in any edges that have a conduction fraction above 1.
#
#   @param  (CellularAutomaton) celAutom
#           An object of the class CellularAutomaton that carries an object of
#           type MetaGraph and a priority queue with all the edges that must undergo
#           the transition
#   @post   for all edges with a transition more than 1, we take that surplus and
#           add it to the neighbouring edges
function handleSurplus(celAutom::CellularAutomaton)
    #If the fraction is greater then 1
    while (!isempty(celAutom.edgesA))&&peek(celAutom.edgesA)[2] > 1
        surplus=peek(celAutom.edgesA)[2]-1#fraction surplus with anistropy
        key=dequeue!(celAutom.edgesA)
        #properties of current edge
        dx=get_prop(celAutom.mg, key[1],key[2],:dx)
        ani=get_prop(celAutom.mg,key[1],key[2],:anisotropy)
        #If its purkinje channel
        if (get_prop(celAutom.mg, key[1],:celtype)==2) && (get_prop(celAutom.mg, key[2], :celtype)==2)
            CV = celAutom.purkinje_CV_multiplier*get_prop(celAutom.mg, key[1], :CV)
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
                        CV2 = celAutom.purkinje_CV_multiplier*get_prop(celAutom.mg, key[2], :CV)
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
#   canTransitionInto will return false if and only if conduction can not travel
#   to the given node.
#
#   @param  (CellularAutomaton) celAutom
#           An object of the class CellularAutomaton
#   @param (Int64) node
#           the node to be checked
#   @param (Float64) timeFraction
#           The fraction of the timestep it takes to get to the node
#   @return true if and only if the time of arrival at the node is later than
#           the time it can be activated
function canPassCurrentTo(celAutom::CellularAutomaton, node::Int64,timeFraction::Float64)
    if get_prop(celAutom.mg,node,:state)==1
        return true
    end
    timeState3=get_timestate3(celAutom,node)
    if get_prop(celAutom.mg,node,:state)==3
        return (get_prop(celAutom.mg,node,:tcounter)+timeFraction>=timeState3)
    end
    timeState2and3= timeState3+get_prop(celAutom.mg,node,:APD)
    return (get_prop(celAutom.mg,node,:tcounter)+timeFraction>=timeState2and3)
end
##
#   activate sets the node on state 2 and the tcounter on -timeFraction, recalculates CV
#   APD
#
#   @param  (CellularAutomaton) celAutom
#           An object of the class CellularAutomaton
#   @param  (Int64) node
#           The (number of the) node which has to be activated
#   @param (Float64) timeFraction
#           The fraction of the timestep it takes to get to the node
#   @pre    the given node must be able to activate
#           |canPassCurrentTo()
#   @post   The CV and APD will be calculated for the time of arrival:
#           get_prop(celAutom.mg, node, :tcounter)+timeFraction
#   @post   Will set the property state of node on 2 and the property tcounter of
#           node on -timeFraction
function activate!(celAutom::CellularAutomaton, node::Int64, timeFraction::Float64)
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
#   @param  (CellularAutomaton) celAutom
#   @param  (Int64) here
#           The node from which the edge activation should start from
#   @param  (Int64) there
#           The other side of the edge
function canActivateEdge(celAutom::CellularAutomaton, here::Int64, there::Int64)
    try
        celAutom.edgesA[tuple(there, here)]
        return false
    catch
        return true
    end

end
##
#   returns the time a node needs to be in state 3
#
#   @param  (CellularAutomaton) celAutom
#           An object of the class CellularAutomaton
#   @param (Int64) node
#           The (number of the) node
#   @pre    the node must be in state 3
function get_timestate3(celAutom::CellularAutomaton,node::Int64)
    return get_prop(celAutom.mg,node,:APD)/4/celAutom.δt
end
##
#   This function sets the APD of the given node. The used formula was taken
#   from a paper's function fit.
#
#   @param  (CellularAutomaton) celAutom
#           An object of the class CellularAutomaton that carries an object of
#           type MetaGraph.
#   @param  (Int64) node
#           The value of the node in the graph, whose APD needs to be changed.
#   @post   The given node will have its APD property changed to the new value
#           of:  APD_endo*(1-T) + APD_epi*T
#           with
#              • T 0.5 when we have a fixed APD and the :Temp when variable (celAutom.variable_APD_flag)
#              • APD_epi = celAutom.APD_ss_epi - celAutom.a_epi*exp(-CL/celAutom.b_epi)
#              • APD_endo = celAutom.APD_ss_endo - celAutom.a_endo*exp(-CL/celAutom.b_endo)
function set_APD!(celAutom::CellularAutomaton, node::Int64)
    #CL = DI + APD (previous cycle) to ms
    CL = (get_prop(celAutom.mg, node, :tcounter) + get_prop(celAutom.mg, node, :APD))*celAutom.δt
    #minimum value of CL for APD, this only happens when invoked
    #under conditions which are not normal for a heart
    if CL<300
        CL=300
    end
    if celAutom.variable_APD_flag
        T = get_prop(celAutom.mg, node, :Temp)
    else
        T=0.5
    end

    #epi APD
    APD_epi = celAutom.APD_ss_epi - celAutom.a_epi*exp(-CL/celAutom.b_epi)#Implementation formula in ms
    APD_epi=APD_epi/celAutom.δt#ms -> t.u

    #endo APD
    APD_endo = celAutom.APD_ss_endo - celAutom.a_endo*exp(-CL/celAutom.b_endo)#Implementation formula in ms
    APD_endo=APD_endo/celAutom.δt#ms -> t.u

    #Real APD
    APD = APD_endo*(1-T) + APD_epi*T
    set_prop!(celAutom.mg,node,:APD, APD)
end
##
#   This function sets the CV of the given node.
#
#   @param  (CellularAutomaton) celAutom
#           An object of the class CellularAutomaton that carries an object of
#           type MetaGraph.
#   @param  (Int64) node
#           The value of the node in the graph, whose APD needs to be changed.
#   @post   The given node will have its CV property changed to the new value
#           of: CV = celAutom.CV_SS-52.12*e^(-DI/87.6)
function set_CV!(celAutom::CellularAutomaton, node::Int64)
    DI=get_prop(celAutom.mg, node, :tcounter)*celAutom.δt#DI in milisec
    #minimum value of DI for CV, this only happens when invoked
    #under conditions which are not normal for a heart
    if DI<40
        DI=40
    end
    CV=celAutom.CV_SS-52.12*exp(-DI/87.6)#cm/sec
    CV=CV*celAutom.δt/1000#transition to cm/t.u.
    CV=CV/celAutom.δx#transition to s.u./t.u.
    set_prop!(celAutom.mg, node,:CV, CV)
end
