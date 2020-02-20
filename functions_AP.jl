################################################################################
####### All functions necessary for doing the AP simulation in Julia
#######
####### Warning:  Make sure all packages are installed and compiled in REPL
#######
################################################################################

# load packages
using LightGraphs
using MetaGraphs
using DelimitedFiles
include("plot_graph.jl")

# create type structure to set model- and simulation- specific details

mutable struct Model
    k::Float64
    a::Float64
    ϵ0::Float64
    μ1::Float64
    μ2::Float64
    dim::Int
    D::Float64 # diffusion coefficient in isotropic case
end

mutable struct Simulation
    δx::Float64 # place discretizations in isotropic case
    δt::Float64 # time step in ms
    t_sample::Int # ever t_sample integration steps written to file
    t_display::Int # every t_display integration steps a plot
    Ntime::Int # number of integration steps
    reaction::Bool #1: reaction part included, 0: pure diffusion
    diffusion::Bool # if true diffusion part included
end

# create type structure to store information to add a stimulus or block a region
mutable struct Stimulus
    t_stimulation::Int # time at which stimulus is given
    val::Float64 # value at which the u is set
    x_start::Float64 # region of stimulation
    x_stop::Float64
    y_start::Float64
    y_stop::Float64
    z_start::Float64
    z_stop::Float64
end

mutable struct Block
    add_block::Bool #true if block is added, false if no block
    t_block::Int #time at which region is blocked
    val::Float64 #value at which the v is set
    x_start::Float64 #region of block
    x_stop::Float64
    y_start::Float64
    y_stop::Float64
    z_start::Float64
    z_stop::Float64
end

################################################################################
########
######## function that checks if the CFL condition is satisfied
########
######## input :m, object of type model containing model-specific parameters
######## s, object of type simulation containing the simulation-specific param.
########
######## output : when the condition is not satisfied δt will be changed
########
########
################################################################################
# check if CFL condition is satisfied
function check_cfl(m::Model,s::Simulation)
    if s.δt > (s.δx)^2/(2*m.dim*m.D)
        dt_cfl = 0.95*(s.δx)^2/(2*m.dim*m.D)
        s.δt = 1/ceil(1/dt_cfl)
        #aanpassen dt en zeggen, dt = 0.95*factor
        println("Timestep δt decreased to ", s.δt)
    end
end

################################################################################
########
######## function that performs the diffusion part of the algorithm to solve
######## the PDE. It does this by calculating the laplacian (5-point stencil).
######## input :g, a metagraph contructed before.
########
######## output : The labels in the graph will be updated.
########
################################################################################
function diffusion_step!(g::MetaGraph)
    vex = collect(vertices(g))
    # loop over vertices
    for i in range(1,stop=size(vex,1))
        # generate u-values and list of neigbours for one vertex
        u_i = get_prop(g,vex[i],Symbol(":u"))
        n = neighbors(g,vex[i])
        du = 0
        # loop over neigbours
        for x in range(1,stop=size(n,1))
            # generate u value neighbour
            u_x = get_prop(g,n[x],Symbol(":u"))
            # generate the coupling factor by getting properties of edge between vertex and its neighbour
            δx = get_prop(g,vex[i],n[x],Symbol(":dx"))
            ani = get_prop(g,vex[i],n[x],Symbol(":anisotropy"))
            D = get_prop(g,vex[i],n[x],Symbol(":diffcoef"))
            f = ani*D/(δx^2)
            # add the calculated term of the laplacian
            du += f*(u_x-u_i)
        end
        set_prop!(g,vex[i],Symbol(":ut"),du)
    end

    #update
    for i in range(1,stop=size(vex,1))
        du =  get_prop(g,vex[i],Symbol(":ut"))
        u = get_prop(g,vex[i],Symbol(":u"))
        u += simulation.δt*du
        set_prop!(g,vex[i],Symbol(":u"),u)
    end

end

################################################################################
########
######## function that performs the reaction part of the algorithm to solve
######## the PDE.
######## input :g, a metagraph contructed before
######## m, a model that contains the necessary parameters
########
######## output : The labels in the graph will be updated.
########
################################################################################
function reaction_step!(g::MetaGraph,m::Model)
    vex = collect(vertices(g))
    for i in range(1,stop=nv(g))
        # u part
        u = get_prop(g,vex[i],Symbol(":u"))
        v = get_prop(g,vex[i],Symbol(":v"))
        du  = - m.k*u*(u - (m.a))*(u-1) - u*v
        # v part
        ϵ = m.ϵ0 + (m.μ1*v/(u + m.μ2))
        dv = ϵ*(-v - m.k*u*(u-m.a-1))

        u += simulation.δt*du
        v += simulation.δt*dv
        set_prop!(g,vex[i],Symbol(":u"),u)
        set_prop!(g,vex[i],Symbol(":v"),v)
    end
end

################################################################################
########
######## function that performs the simulation according to the AP model. The
######## parameters of the model are specified in the Type Model and the
######## details of the simulation are in the Type Simulation.
######## input : m::Model, object of type Model containing info of model
######## s::Simulation, object of type Simulation containing info of simulation
######## g, a metagraph contructed before.
########
######## output : The labels in the graph will be updated every timestep. An
######## output file with the sampled values of u will be generated for every
######## timestep. Frames will be saved of the metagraph plotted in space.
########
################################################################################
function simulate!(m::Model,s::Simulation,g::MetaGraph,stimulus::Stimulus,block::Block)
    # initialisation
    vex = collect(vertices(g))
    v_ecg = zeros(s.Ntime)
    u_old = 0
    v_old = 0
    if !s.reaction
        for i in range(1, stop=nv(g))
            set_prop!(g,vex[i],Symbol(":vt"),0)
        end
    end
    if !s.diffusion
        for i in range(1, stop=nv(g))
            set_prop!(g,vex[i],Symbol(":ut"),0)
        end
    end
    # main loop over time step
    for t in range(1, stop=s.Ntime)
        if t==stimulus.t_stimulation
            add_stimulus!(g,stimulus)
        end
        if block.add_block && t==block.t_block
            add_block!(g,block)
        end
        # calculate laplacian/diffusion part for this timestep
        if s.diffusion
            diffusion_step!(g)
        end
        # calculate reaction part of this timestep
        if s.reaction
            reaction_step!(g,m)
        end

        # sampling
        if t%s.t_sample == 0
            println(t)
            # write to file
            writedlm("output/sampling_data$t.txt",["nodes" "u" "v"],',')
            for i in range(1,stop=nv(g))
                open("output/sampling_data$t.txt","a") do io
                    writedlm(io,[vex[i] get_prop(g,vex[i],Symbol(":u")) get_prop(g,vex[i],Symbol(":v"))],',')
                end
            end
        end
        # calculate ecg
        v_ecg[t] = ecg(g)
        # plotting
        if t%s.t_display == 0
            # potential on graph
            if m.dim == 3
                plot_graph_extended(g,m,t)
            else
                plot_the_graph(g,t)
                p1 = plot(load("frame$t.png"),show=true,grid=false,showaxis=false)
                display(p1)
            end
            # ECG
            x_ax = (t-s.t_display+1):t
            #p2 = plot(collect(x_ax),v_ecg[x_ax],label = "ECG $t", grid=false,show=true)
            # display(p2)
            #png(p2,"ECG$t.png")
        end

    end
end
################################################################################
########
######## function that reads in the data files and constructs a graph
######## input : filename_vertices::String, filenames of the data files
######## containing information of the vertices;
######## filename_edges::String, filename with the info of the edges.
######## delimiter::Char, specify how the different points in a row are
######## separated, for instance ',' or '\t' or ' '
########
######## Warning: make sure the files are in the same path as from where you
######## call this function.
########
######## output : MetaGraph containing all the vertices and edges with labels.
########
################################################################################
function construct_graph(filename_vertices::String,filename_edges::String,delimiter::Char)
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
################################################################################
########
######## function that adds a stimulus to the system. It will set the u variable
######## to one and v (recovery) to zero in an area specified by the input.
########
######## input : g, metagraph object
######## the range specified by x_start, x_stop, y_start, y_stop, z_start, z_stop
######## Warning : the range of the coordinates start at one!
########
######## output : The labels of the metagraph will be updated (u and v)
########
################################################################################
function add_stimulus!(g::MetaGraph,stimulus::Stimulus)
    vex_range = get_area(g,stimulus.x_start,stimulus.x_stop,stimulus.y_start,stimulus.y_stop,stimulus.z_start,stimulus.z_stop)
    println("stimulus added ", vex_range[1], " to ",vex_range[end])
    for i in range(1,stop=size(vex_range,1))
        set_prop!(g,vex_range[i],Symbol(":u"),stimulus.val)
        set_prop!(g,vex_range[i],Symbol(":v"),0)
    end
end

################################################################################
########
######## function that adds a block to the system. It will set the u variable
######## to zeros and v (recovery) to one in an area specified by the input.
########
######## input : g, metagraph object
######## the range specified by x_start, x_stop, y_start, y_stop, z_start, z_stop
########
######## output : The labels of the metagraph will be updated (u and v)
########
################################################################################
function add_block!(g::MetaGraph,block::Block)
    vex_range = get_area(g,block.x_start,block.x_stop,block.y_start,block.y_stop,block.z_start,block.z_stop)
    println("block added ", vex_range[1], " to ",vex_range[end])
    for i in range(1,stop=size(vex_range,1))
        set_prop!(g,vex_range[i],Symbol(":u"),0)
        set_prop!(g,vex_range[i],Symbol(":v"),block.val)
    end
end

################################################################################
########
######## function that returns the indices of the vertices that span the
######## range/area specified by the inputs.
########
######## input : g, metagraph object
######## the range specified by x_start, x_stop, y_start, y_stop, z_start, z_stop
########
######## output : array I of indices of vertices
########
################################################################################
function get_area(g::MetaGraph,x_start::Float64,x_stop::Float64,y_start::Float64,y_stop::Float64,z_start::Float64,z_stop::Float64)
    I = []
    for i in range(1,stop=nv(g))
         if x_start <= get_prop(g,i,Symbol(":loc_x")) <= x_stop
             if y_start <= get_prop(g,i,Symbol(":loc_y")) <= y_stop
                 if z_start <= get_prop(g,i,Symbol(":loc_z")) <= z_stop
                     append!(I, i)
                 end
             end
         end
    end
    return I
end

