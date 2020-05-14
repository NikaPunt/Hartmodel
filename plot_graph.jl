################################################################################
########
######## Script with all functions necessary for plotting.
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
include("CellularAutomatonMutableStruct.jl")

##
# function that returns the coordinates of the vertices.
#
#   @param (MetaGraph) mg
#   @return 3 lists that contain the x-,y- and z-values resp of each vertex
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
##
#   plotGraph2 uses the MetaGraph object embedded in the given
#   CellularAutomaton object (celAutom) to create an image of the (colored)
#   MetaGraph object.
#
#   @param  (CellularAutomaton) celAutom
#           An object of the class CellularAutomaton that carries an object of
#           type MetaGraph.
#   @param  (Int64) i
#           This is the number assigned to the end of the name of the
#           picture file.
#   @param  (String) folder
#           this is the name of the folder used to store the pictures. If it
#           exists it doesn't attempt to make a new one. If it doesn't exist,
#           the directory will be created.
#   @post   Generates an image of the colored MetaGraph object using the gplot
#           function and stored as $folder/frame$i.png.
#   @post   this will generate a plot of a 2D graph.
function plotGraph2(celAutom::CellularAutomaton,i::Int64,folder::String)
    nodefillc = coloringNodes(celAutom)
    edgefillc =coloringEdge(celAutom)
    loc_x,loc_y,loc_z = get_coordinates(celAutom.mg)
    g1=gplot(celAutom.mg,loc_x,loc_y,nodefillc=nodefillc,edgestrokec=edgefillc)
    draw(PNG("$folder/frame$i.png", 16cm, 16cm), g1)
end
##
#   coloringNodes will return an array of colors that can color the nodes of
#   the MetaGraph embedded in the given CellularAutomaton. This array can be
#   used as the nodefillc array used in the gplot function.
#
#   @param  (CellularAutomaton) celAutom
#           An object of the class CellularAutomaton that carries an object of
#           type MetaGraph.
#   @post   nodefillc will assign the color "black" to nodes in state one,
#           "maroon" to nodes in state two and "red4" to nodes in state three.
#   @return Returns "nodefillc", an array of type Array{RGB{Normed{UInt8, 8}},1}
#           with the needed coloring for the gplot function.
function coloringNodes(celAutom::CellularAutomaton)
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
#   coloringEdge will return an array of colors that can color the edges of
#   the MetaGraph embedded in the given CellularAutomaton. This array can be
#   used as the edgestrokec array used in the gplot function.
#
#   @param  (CellularAutomaton) celAutom
#           An object of the class CellularAutomaton that carries an object of
#           type MetaGraph.
#   @post   The edges will be colored in a spectrum from white (not excited at all)
#           to dark red (completely excited),
#           depending on the ltransition and/or htransition in the edge.
#           Edges recovering from excitation will be colored green.
#   @return Returns colors, an array of type Array{RGB{Normed{UInt8, 8}},1}
#           with the needed coloring for the gplot function.
function coloringEdge(celAutom::CellularAutomaton)
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
#   @post   The metagraph object will get all the properties described in the
#           file names.
#   @return Returns a MetaGraph object with the properties described in
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
    println("Graph check!")
    return mg
end
##
# function that returns the indices of the vertices that span the
# range/area specified by the inputs.
#
#   @param (Metagraph) object
#   @param the range specified by x_start, x_stop, y_start, y_stop, z_start, z_stop
#   @return array I of indices of vertices
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
