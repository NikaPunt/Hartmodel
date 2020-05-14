mutable struct CellularAutomaton
    mg::MetaGraph
    #properties nodes:
        #1: non-active and can be depolarised -> DI (black)
        #2: active and fires in begin period -> APD (red)
        #3: non-active and can't be depolarised -> DI (green)
    time::Int64#Time since running time in t.u.
    δt::Int64 #the time unit in ms -> 1 t.u. = δt ms

    δx::Float64 #the space unit in cm -> 1 s.u. = δx cm

    #Values of the APD-restiutioncurve
    APD_ss_epi::Float64
    APD_ss_endo::Float64

    a_epi::Float64
    a_endo::Float64

    b_epi::Float64
    b_endo::Float64

    #Values of the CV_restitutioncurve
    CV_SS::Float64#70.03

    #transitions of the active edges
    edgesA::PriorityQueue{Tuple{Int64,Int64}, Float64}

    #flag: do we use purkinje nodes?
    purkinje_flag::Bool
    #exit nodes for purkinje
    LVexit::Int64
    RVexit::Int64
    #time to start the excitation of RV_exit
    RVexit_time::Int64
    #multiplier value for the CV of a purkinje channel
    purkinje_CV_multiplier::Float64

    #flag: do we use variable APD-restitution
    variable_APD_flag::Bool

    aanzicht::Int64
end

##
#   createCellularAutomaton will create an instance of CellularAutomaton,
#   with the given graph, and node(s) with initial excitation.
#
#   @param  (MetaGraph) graph
#           The graph to be embedded in the new cellular automaton.
#   @param  (Array{Int64, 1}) startwaarden
#           The values of the nodes that have an initial excitation.
#   @param (Array{Int64, 1}) stopwaarden
#           The values of the nodes that can't be depolarised in te inital condition
#
#   The next six parameters are parameters of the formula APD = APD_ss - a*exp(-DI*celAutom.δt/b)
#   @param  (Float64) APD_ss_epi
#           The steady state APD of the epicardium
#   @param  (Float64) APD_ss_endo
#           The steady state APD of the endocardium
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
#   @param  (Float64) CV_SS
#           The CV that we want in a node in steady state.
#   @param  (Float64) CV_purkinje
#           The constant CV between two nodes with :celtype property 2.
#   @param  (Bool) purkinje_flag
#           A flag that represents whether we need to excite LVexit and RVexit.
#   @param  (Bool)  variable_APD_flag
#           A flag that represents whether the APD should be variable or not.
#   @pre    stopwaarden and startwaarden don't have common values
#   @post   The embedded MetaGraph has the property state imposed on the nodes.
#           The state of all nodes except the values in startwaarden are set to '1' by default.
#           Whilst the state of nodes in startwaarden are set to '2' (and if
#           purkinje_flag is true, also the state of LVexit).
#   @post   If purkinje_flag is set to true, the tcounter property of all nodes
#           except the ones in startwaarden and LVexit is set to -1000.
#           If purkinje_flag is false, the tcounter of all nodes, except for the
#           ones in startwaarden, is set to 200/δt.
#   @post   The tcounter of all nodes in startwaarden and stopwaarden is set to 0.
#   @post   The state property of all nodes in startwaarden and stopwaarden is set to 2.
#   @post   If purkinje_flag is set to true, the state of LVexit is set to 2.
#   @post   For each node in startvalues we collect the neighbors of that node.
#           We check whether or not those neighbors belong to startwaarden or
#           stopwaarden (if purkinje_flag is true we also check if they're not
#           LVexit or RVexit). If not, we put tuple(node, neighbor) as a key in
#           our priority queue with the value 0. This tuple therefore also
#           respects the direction the current runs.
#   @post   If purkinje_flag is set to true, we collect the neighbors of LVexit.
#           We check whether or not those neighbors belong to startwaarden,
#           stopwaarden or is RVexit. If not, we put tuple(node, neighbor) as a
#           key in our priority queue with the value 0.
function createCellularAutomaton(graph::MetaGraph, δt::Int64, δx::Float64, startwaarden::Array{Any,1},
                                    stopwaarden::Array{Any,1}, APD_ss_endo::Float64,
                                    APD_ss_epi::Float64, a_epi::Float64, a_endo::Float64,
                                    b_epi::Float64, b_endo::Float64,
                                    LVexit::Int64, RVexit::Int64, RVexit_time::Int64, purkinje_CV_multiplier::Float64,
                                    CV_SS::Float64, purkinje_flag::Bool,variable_APD_flag::Bool, aanzicht::Int64)
    Priority = PriorityQueue{Tuple{Int64,Int64}, Float64}(Base.Order.Reverse)
    for node in startwaarden
        for buur in collect(neighbors(graph, node))
            #anders problemen met edges
            if purkinje_flag
                if !(buur in startwaarden) &&!(buur in stopwaarden) &&!(buur==RVexit) && !(buur == LVexit)
                    enqueue!(Priority,(node, buur),0)
                end
            else
                if !(buur in startwaarden) &&!(buur in stopwaarden)
                    enqueue!(Priority,(node, buur),0)
                end
            end

        end
    end
    #if purkinje_flag == true we can excite LVexit.
    if purkinje_flag
        for buur in collect(neighbors(graph, LVexit))
            #anders problemen met edges
            if !(buur in startwaarden) &&!(buur in stopwaarden) &&!(buur==RVexit)
                enqueue!(Priority,(LVexit, buur),0)
            end
        end
    end
    celAutom = CellularAutomaton((mg=graph,time=0,δt=δt,δx=δx, APD_ss_epi=APD_ss_epi,
                                APD_ss_endo=APD_ss_endo, a_epi=a_epi, a_endo=a_endo,
                                b_epi=b_epi, b_endo=b_endo, CV_SS=CV_SS, edgesA=Priority,
                                purkinje_flag=purkinje_flag, LVexit=LVexit, RVexit = RVexit,
                                RVexit_time=RVexit_time, purkinje_CV_multiplier=purkinje_CV_multiplier,
                                variable_APD_flag=variable_APD_flag, aanzicht=aanzicht)...)
    for node in collect(vertices(graph))
        set_prop!(graph,node,:state,1)
        set_prop!(graph, node,:APD,0)
        if purkinje_flag
            set_prop!(graph,node,:tcounter,200/δt)
        else
            #CL negative but will invoke the methods so that their minimum values will be induced
            set_prop!(graph,node,:tcounter,-1000)
        end
        set_APD!(celAutom,node)
        set_CV!(celAutom,node)
    end

    for node in startwaarden
        set_prop!(graph,node,:state,2)
        set_prop!(graph,node,:tcounter,0)
    end
    #if purkinje_flag == true we can excite LVexit.
    if purkinje_flag
        set_prop!(graph,LVexit,:state,2)
        set_prop!(graph,LVexit,:tcounter,0)
    end
    for node in stopwaarden
        set_prop!(graph,node,:state,2)
        set_prop!(graph,node,:tcounter,0)
    end
    return celAutom
end
