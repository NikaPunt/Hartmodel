################################################################################
########
######## Script with all functions necessary for plotting. This is loaded in
######## functions_AP.jl
########
########
################################################################################

using LightGraphs
using MetaGraphs
using GraphPlot
using Colors
using Compose
using Cairo
using Fontconfig
using Images
using Plots
using GraphRecipes
# theme(:juno)


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
    δt::Float64
    t_sample::Int
    t_display::Int
    Ntime::Int
    reaction::Bool #1: reaction part included, 0: pure diffusion
    diffusion::Bool
end

################################################################################
########
######## function that assigns color ID to the u-values in the metagraph
########
######## input :mg, a metagraph
######## levels, an array that contains all the levels (so ranges of the colors)
########
######## output : An array with for each vertex a color ID.
########
################################################################################
function colorize(mg::MetaGraph,levels::Array{Float64,1})
    col = zeros(Int,nv(mg))
    for vi in collect(vertices(mg))
        u = get_prop(mg,vi,:u)
        col[Int(vi)] = searchsortedfirst(levels,u)
    end
    return col
end
################################################################################
########
######## function that returns the coordinates of the vertices.
########
######## input :mg, a metagraph
########
######## output: 3 lists that contain the x-,y- and z-values resp of each vertex
########
################################################################################
function get_coordinates(mg::MetaGraph)
    loc_x = zeros(nv(mg))
    loc_y = zeros(nv(mg))
    loc_z = zeros(nv(mg))
    for i in range(1,stop=nv(mg))
        loc_x[i] = get_prop(mg,i,:loc_x)
        loc_y[i] = get_prop(mg,i,:loc_y)
        loc_z[i] = get_prop(mg,i,:loc_z)
    end
    return loc_x,loc_y,loc_z
end
################################################################################
########
######## function that generates a frame of the graph at one timestep. It saves
######## the gplot as a PNG. This is done because this can then afterwards be
######## plotted and used in a gif. This is not supported for gplot.
########
######## input :mg, a metagraph
########
######## output : frame.png
########
################################################################################
function plot_the_graph(mg::MetaGraph,t::Int32)
    levels = [-2.0,-1.0,0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,5.0,10.0]
    membership = colorize(mg,levels)
    c1 = colorant"red"
    c2 = colorant"yellow"
    nodecolor = [colorant"black", colorant"maroon",colorant"red4",colorant"darkred",colorant"firebrick4",
    colorant"firebrick",colorant"orangered3",colorant"orangered",colorant"darkorange1",colorant"orange",
    colorant"darkgoldenrod2",colorant"darkgoldenrod1",colorant"goldenrod1",colorant"gold",colorant"yellow2",
     colorant"yellow"]
    nodefillc = nodecolor[membership]
    loc_x,loc_y,loc_z = get_coordinates(mg)
    g1 = gplot(mg,loc_x,loc_y,nodefillc=nodefillc)
    draw(PNG("frame$t.png", 16cm, 16cm), g1)
end


################################################################################
########
######## function that plots the graph in 3D at one timestep.
########
######## input :mg, a metagraph
########
######## output : PNG and plot
########
################################################################################
function plot_graph_extended(g::MetaGraph,m::Model,t)
    X,Y,Z = get_coordinates(g)
    levels = collect(0:.1:1.0)
    membership = colorize(g,levels)
    c = range(colorant"blue", stop=colorant"red", length=size(levels,1)+1)
    markerc = c[membership]
    g2 = graphplot(g,dim=m.dim, x = X, y=Y, z=Z, markercolor=markerc,linecolor=:black,curves=false)
    p = plot(g2,show=true,grid=false,showaxis=false)
    display(p)
    png(p, "3Dframe$t.png")
    # draw(PNG("3Dframe$t.png", 16cm, 16cm), g2)
end

################################################################################
########
######## function that calculates the voltage output measured by ECG, i.e.
######## the sum of the u values of all the vertices
########
######## input :mg, a metagraph
########
######## output : a Float sum_u
########
################################################################################
function ecg(g::MetaGraph)
    sum_u = 0
    for i in range(1,stop=nv(g))
        sum_u += get_prop(g,i,:u)
    end
    return sum_u
end

################################################################################
########
######## function that plots the u and v values over the range of vertices
######## (so in space, for one fixed time)
########
######## input :filename_sampling, the filename of the data file that contains
######## the u,v for all the vertices (normally sampling_data$i.txt)
########
######## output : two plots, one for u and for v
########
################################################################################
function plot_uv(filename_sampling::String)
    input_vertices = readdlm(filename_sampling,',',skipstart = 1)
    p1 = plot(input_vertices[:,1],input_vertices[:,2], color="blue",label="u")
    xlabel!("vertices")
    ylabel!("u")
    p2 = plot(input_vertices[:,1],input_vertices[:,3], color="red",label="v")
    xlabel!("vertices")
    ylabel!("v")
    display(p1)
    display(p2)
end

################################################################################
########
######## function that plots the u and v values of one vertex in time
########
######## input : vertex, the vertex you want the u v values from
######## s, an object of type Simulation
########
######## output : two plots, one for u and for v
########
################################################################################
function plot_utvt(vertex::Int,s::Simulation)
    N = s.Ntime
    dt = s.t_sample
    ut = zeros(Int(N/dt))
    vt = zeros(Int(N/dt))
    cnt = 1
    while cnt <= Int(N/dt)
        t = cnt*dt
        #println(t)
        t = Int(t)
        filename = "output/sampling_data$t.txt"
        data = readdlm(filename,',',skipstart=1)
        ut[cnt] = data[vertex,2]
        vt[cnt] = data[vertex,3]
        cnt += 1
    end
    t = range(1,stop=s.Ntime,length=Int(N/dt))
    p1 = plot(t,ut,color="blue",label="u")
    xlabel!("time t")
    ylabel!("u")
    p2 = plot(t,vt,color="red",label="v")
    xlabel!("time t")
    ylabel!("v")
    display(p1)
    display(p2)
    png(p1, "AP.png")
end
