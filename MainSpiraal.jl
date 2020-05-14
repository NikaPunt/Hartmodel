include("CellularAutomaton.jl")
start=time()
graph = constructGraph("data_vertices.dat", "data_edges.dat", ',')

startwaarden=get_area(graph, -100.0,115.0, 110.0,125.0,-1000.0,1000.0)
stopwaarden =get_area(graph, -100.0,115.0,125.0,140.0,-1000.0,1000.0)


#epi APD
APD_ss_epi = Float64(392.61)
a_epi = Float64(339.39)
b_epi = Float64(520)

#endo APD
APD_ss_endo = Float64(505.74)
a_endo = Float64(485.4)
b_endo = Float64(501)

#Steady state value of CV
CV_SS = 70.03

#timestep in ms
δt=20
#spacestep in cm
δx = 1/10

#the multiplier for the CV
purkinje_CV_multiplier = Float64(2)

LVexit = Int64(8427)
RVexit = Int64(5837)
time_ms = Int64(20) #na 20 miliseconden moet RVexit oplichten
RVexit_time = Int64(ceil(time_ms/δt)) #aantal tijdstappen voordat RVexit moet oplichten.
purkinjeflag = false

variable_APD_flag = true

celAutom = createCellularAutomaton(graph, δt, δx, startwaarden,stopwaarden,
                    APD_ss_endo, APD_ss_epi, a_epi, a_endo, b_epi, b_endo,
                    LVexit, RVexit, RVexit_time, purkinje_CV_multiplier,
                    CV_SS, purkinjeflag, variable_APD_flag)

folder="plotjesxyDiepte"

createFrames(folder,200,200,celAutom)
elapsed = time() - start
println(elapsed)
