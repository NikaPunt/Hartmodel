# Code voor cellulaire automaat
using LightGraphs
using MetaGraphs
using DelimitedFiles
include("plot_graph.jl")

function plot_graph_extended(mg::MetaGraph,i::Int64,folder::String)
    loc_x,loc_y,loc_z = get_coordinates(mg)
    #nodefillc = coloringGraph(celAutom)
    #edgefillc =coloringEdge(celAutom)
    g2 = graphplot(mg,dim=3, x = loc_x, y=loc_y, z=loc_z,markercolor = colorant"red",linecolor = colorant"green",markersize = 0.01, curves=false)
    p = plot(g2,show=true,grid=false,showaxis=false)
    display(p)
    png(p, "$folder/3Dframe$i.png")
    #draw(PNG("$folder/3Dframe$i.png", 16cm, 16cm), g2)

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
#   @post   Returns a MetaGraph object with the properties described in
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
function main()
    start = time()
    graph = constructGraph("data_vertices.dat", "data_edges.dat", ',')
    folder="plotjes_test3D"
    plot_graph_extended(graph,1,folder)
    elapsed = time() - start
    println(elapsed)
end

main()
