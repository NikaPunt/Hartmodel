################################################################################
########
######## Main file that calls necessary functions specified in functions_AP.jl
######## and plot_graph.jl and thus performs simulation. Make sure these two
######## scripts are compiled in Julia before executing the Main.jl
########
######## questions: lore.leenknegt@kuleuven.be
######## code version 1.0 Lore Leenknegt (12/02/2020) git:
################################################################################
# load in scripts with all functions
include("functions_AP.jl")
cd("C:/Users/nikap/OneDrive/Desktop/code_students")
# load in data file
graph = construct_graph("data_vertices_test.dat","data_edges_test.dat",',')

# set details model and simulation
model = Model((k = 8.0,a = 0.15,ϵ0 = 0.002,μ1 = 0.2,μ2 = 0.3 ,dim = 3,D = 20.0)...)
simulation = Simulation((δx = 1, δt = 0.1, t_sample = 1000, t_display = 100, Ntime = 5000, reaction = true, diffusion = true)...)
stimulus = Stimulus((t_stimulation=1,val=0.9,x_start=0.0*simulation.δx,x_stop=30.0*simulation.δx,y_start=1.0*simulation.δx,y_stop=1.0*simulation.δx,z_start=1.0*simulation.δx,z_stop=1.0*simulation.δx)...)
block = Block((add_block=false,val=1.0,t_block=1,x_start=0.0*simulation.δx,x_stop=0.0*simulation.δx,y_start=0.0*simulation.δx,y_stop=0.0*simulation.δx,z_start=0.0*simulation.δx,z_stop=0.0*simulation.δx)...)
check_cfl(model,simulation)


# print and write info
println("#####################################################################")
println("k = ", model.k, "\n", "a = ",model.a, "\n", "ϵ0 = ",model.ϵ0,"\n","μ1 = ", model.μ1, "\n","μ2 = ",
 model.μ2, "\n", "Ntime = ", simulation.Ntime, "\n", "t_sample = ", simulation.t_sample, "\n","t_display = ", simulation.t_display,"\n","dim = ",model.dim,"\n","D = ",model.D)
info_1 = ["k","a","ϵ0","μ1","μ2","Ntime","t_sample","t_display","Nnodes","dim","D"]
info_2 = [model.k, model.a, model.ϵ0, model.μ1, model.μ2, simulation.Ntime, simulation.t_sample, simulation.t_display, nv(graph), model.dim, model.D]
writedlm("output/info_simulation.txt",[info_1 info_2], ' ')
println("#####################################################################")
println("Simulation initialized, start solving differential equation now...")

# do simulation
simulate!(model,simulation,graph,stimulus,block)

# plotting to test output
# plot_uv("output/sampling_data1.txt")
plot_utvt(50,simulation)

println("Finished!")
println("#####################################################################")
