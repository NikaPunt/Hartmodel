include("CellularAutomaton.jl")
function Main()
    start=time()
    graph = constructGraph("data_vertices.dat", "data_edges.dat", ',')

    startwaarden=get_area(graph, -100.0,115.0, 110.0,125.0,-1000.0,1000.0)
    stopwaarden = get_area(graph, -100.0,115.0,125.0,140.0,-1000.0,1000.0)


    #epi APD
    APD_ss_epi = Float64(392.61)
    a_epi = Float64(339.39)
    b_epi = Float64(520)

    #endo APD
    APD_ss_endo = Float64(505.74)
    a_endo = Float64(485.4)
    b_endo = Float64(501)

    #timestep in ms
    dt=20
    #spacestep in cm
    δx = 1/10

    #CV_purkinje
    CV_purkinje_mps = 2 #m/s
    CV_purkinje_sups = CV_purkinje_mps*100/δx #s.u./s
    CV_purkinje_suptu= CV_purkinje_sups/dt #s.u./t.u.

    LVexit = Int64(8427)
    RVexit = Int64(5837)
    time_ms = Int64(20) #na 20 miliseconden moet RVexit oplichten
    RVexit_time = Int64(ceil(time_ms/dt)) #aantal tijdstappen voordat RVexit moet oplichten.
    purkinjeflag = true

    variable_APD_flag = true

    celAutom = createCellulaireAutomaat(graph, dt, δx, startwaarden,stopwaarden,
                        APD_ss_endo, APD_ss_epi, a_epi, a_endo, b_epi, b_endo,
                        LVexit, RVexit, RVexit_time, CV_purkinje_suptu,
                        purkinjeflag, variable_APD_flag)

    folder="plotjesxyDiepte"

    createFrames(folder,100,100,celAutom)
    elapsed = time() - start
    println(elapsed)
end

Main()
