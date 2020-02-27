function coloringEdge(celAutom::CellulaireAutomaat)

    #ne = number of edges
    colors = [colorant"white" for i in 1:ne(celAutom.mg)]

    j = 0
    for edge in collect(edges(celAutom.mg))
        j += 1
        prop = get_prop(celAutom.mg,edge,:htransition)
        CV = get_prop(celAutom.mg,edge.dst,:CV)
        dx = get_prop(celAutom.mg,edge,:dx)
        if (ltransition + htransition)/dx == 0
            colors[j]= colorant"white"
        elseif (ltransition + htransition)/dx < 0.2
            colors[j]= colorant"grey"
        elseif (ltransition + htransition)/dx < 0.4
            colors[j]= colorant"blue"
        elseif (ltransition + htransition)/dx < 0.6
            colors[j]= colorant"green"
        elseif (ltransition + htransition)/dx < 0.8
            colors[j]= colorant"orange"
        elseif (ltransition + htransition)/dx == 1
            colors[j]= colorant"red"
        end
    end
    return colors = colors
end
