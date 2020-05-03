# Code voor cellulaire automaat
using LightGraphs
using MetaGraphs
using DelimitedFiles
using DataStructures
include("CellulaireAutomaat.jl")
graph = constructGraph("nieuwdata_vertices.dat", "nieuwdata_edges.dat", ',')
startwaarden=[1297]
stopwaarden = []

#epi APD
ARI_ss_epi = Float64(392.61)
a_epi = Float64(339.39)
b_epi = Float64(520)

#endo APD
ARI_ss_endo = Float64(505.74)
a_endo = Float64(485.4)
b_endo = Float64(501)

compileT=[]
for timestep in range(1,stop=100)
    dt=Int64(timestep)
    celAutom = createCellulaireAutomaat(graph,dt, startwaarden,stopwaarden,
                    ARI_ss_endo, ARI_ss_epi, a_epi, a_endo, b_epi, b_endo)
    compileTries=[]
    for i in range(1,stop=10)
        start = time()
        for i in range(1,stop=ceil(1000/timestep))
            updateState!(celAutom)
            for node in collect(vertices(celAutom.mg))
                set_prop!(celAutom.mg,node,:tcounter,get_prop(celAutom.mg,node,:tcounter)+1)
            end
            celAutom.time+=1
        end
        elapsed = time() - start
        append!(compileTries, elapsed)
    end
    append!(compileT, sum(compileTries)/10)
end
println(compileT)
plot(1:100, compileT, title = "Snelheid programma",xlabel ="Grootte tijdstap" ,ylabel="Compileertijd", lw = 3)
