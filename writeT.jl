
println("####################################")
open("data_vertices_2D_met_T.dat", "w") do fileT
    open("data_vertices_2D.dat") do fileNT
        for ln in eachline(fileNT)
            if ln[1] == 'v'
                write(fileT, "$(ln),T\n")
            else
                write(fileT, "$(ln),1.0\n")
            end
        end
    end
end
