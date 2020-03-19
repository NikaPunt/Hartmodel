f= open("data_edges_test2.dat","w")
write(f,"edge1,edge2,:anisotropy,:dx,:diffcoef\n")
for i in range(1,stop=2450)
        if mod(i,50)==0
                i=float(i)
                write(f,"$i,$(i+50),1.0,1.0,1.0\n")
        else
                i=float(i)
                write(f,"$i,$(i+1),1.0,1.0,1.0\n")
                write(f,"$i,$(i+50),1.0,1.0,1.0\n")
                write(f,"$i,$(i+51),1.0,$float64(sqrt(2)),1.0\n")
        end
end
for i in range(2450,stop=2499)
        i=float(i)
        write(f,"$i,$(i+1),1.0,1.0,1.0\n")
end
close(f)
