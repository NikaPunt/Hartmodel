f= open("data_vertices_test2.dat","w")
write(f,"vertices,:celtype,:loc_x,:loc_y,:loc_z,:u,:v\n")
for i in range(1,stop=50)
    for j in range(1,stop=50)
        i=float(i)
        j=float(j)
        write(f,"$((i-1)*50+j),1.0,$i,$j,0.0,0.0,0.0\n")
    end
end
close(f)
