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
end
