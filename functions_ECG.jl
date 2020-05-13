# loading packages

using LinearAlgebra
using LightGraphs
using MetaGraphs
using DelimitedFiles
using Profile
using Distributions

# If you have questions about the entries of matrices such as tetrahedron, consult CODE_README.

mutable struct ECG
      type::String # "beam","heart3","heart12"
      elec::Dict{String,Array{Float64,1}}
      lead::Dict{String,Array{Float64,1}}
      loc_elec::Array{Float64,2}
      name_elec::Array{String,1}
      wavefront::Array{Float64,1}
      n_calcs::Int64
      p_0::Float64
      p::Array{Float64,1}
end


"""
    find_v_fire(tetraeder::Array{Any,1})\n
Returns the vertices that fire in a particular tetrahedron, as wel as the transition
fractions on the edges.

@param      (Array{Float64,1}) `vertices` \n
            Contains the vertices of the tetrahedron, in the order [A B C D]\n
@param      (PriorityQueue{Tuple{Int64,Int64}, Float64}) `edgesA`
            Contains as keys the active edges as tuples, such that the first edge
            in the tuple has fired. The values are the transition fraction  α of
            the particular edge.

@return     (Set{Int64}) `v_fire`\n
            Indices of the vertices that have fired\n
@return     (Array{Float64,1}) `etransition`\n
            Transition factors of the edges, in the order [AB BD DA CA CB CD]
"""
function find_v_fire(vertices::Array{Int64,1},edgesA::PriorityQueue{Tuple{Int64,Int64}, Float64})
      # all possible edges, in order AB BD DA CA CB CD
      combin_edges = [1 2
                        2 4
                        1 4
                        1 3
                        2 3
                        3 4]
      v_fire = Set()
      etransition = zeros(6) # transition value of edges, in order AB BD DA CA CB CD

      for i = 1:6
            if haskey(edgesA,(vertices[combin_edges[i,1]],vertices[combin_edges[i,2]]))
                  α = edgesA[(vertices[combin_edges[i,1]],vertices[combin_edges[i,2]])]
                  if α !=(0)
                        union!(v_fire,vertices[combin_edges[i,1]]) # first vertex has fired
                  end
                  etransition[i] = α
            elseif haskey(edgesA,(vertices[combin_edges[i,2]],vertices[combin_edges[i,1]]))
                  α = edgesA[(vertices[combin_edges[i,2]],vertices[combin_edges[i,1]])]
                  if α !=(0)
                        union!(v_fire,vertices[combin_edges[i,2]]) # second vertex has fired
                  end
                  etransition[i] = α
            else etransition[i] = 0 # edge is not active
            end
      end

      return (v_fire,etransition)
end


"""
    get_permutation(tetrahedron, edgen, v_fire)
Gives the indices of the order in which the vertices and edges of the tetrahedron
should be renamed such that the active vertices are:\n
* in the case of 1 active vertice: on the 'A' position;\n
* in the case of 2 active vertices: on the 'A' and 'B' position;\n
* in the case of 3 active vertices: on the 'B', 'C' and 'D' position.\n
This in order to easily calculate the vector normal to the wave front.\n

@param      (Array{Any,1}) `vertices` \n
            Contains the vertices of the tetrahedron, in the order [A B C D]\n
@param      (Set{Any}) `v_fire`\n
            Indices of the vertices that have fired\n

@return     (Array{Int64,1}) `indices_tetra`\n
            Contains the indices of the rows of 'tetrahedron' such
            that the tetrahedron with the original rows in this given order
            suffices the conditions specified before.\n
@return     (Array{Float64,1}) `indices_edgen`\n
            Contains the indices of the elements of 'etransition' such
            that the tetrahedron with the original rows in this given order
            suffices the conditions specified before.
"""
function get_permutation(vertices::Array{Any,1},v_fire::Set{Any})

    # alle mogelijke rotaties. Eerste lijst zijn de vertices in nieuwe volgorde ABCD, met oude nummers A=1;B=2;C=3;D=3
    # Tweede lijst bevat de nieuwe volgorde van de edges, in nieuwe volgorde [AB BD DA CA CB CD]
    # Oude labels edges: AB=1,BD=2,DA=3,CA=4,CB=5,CD=6
    rotation1 = ([1 2 3 4],[1 2 3 4 5 6]) # doet niets
    rotation2 = ([2 4 3 1],[2 3 1 5 6 4])
    rotation3 = ([3 1 2 4],[4 3 6 5 1 2])
    rotation4 = ([4 2 1 3],[2 5 6 3 1 4])
    rotation5 = ([4 1 3 2],[3 1 2 6 4 5])
    rotation6 = ([4 3 2 1],[6 4 3 2 5 1])

    if length(v_fire) == 1 # 1 vertex firet, onderscheid per vertex
      if vertices[1] in v_fire # A
            return rotation1
      elseif vertices[2] in v_fire # B
            return rotation2
      elseif vertices[3] in v_fire # C
            return rotation3
      elseif vertices[4] in v_fire # D
            return rotation4
      end

    elseif length(v_fire) == 2 # 2 vertices firen
      if vertices[1] in v_fire # A firet
             if vertices[2] in v_fire # AB
                   return rotation1
             elseif vertices[3] in v_fire # AC
                   return rotation3
             elseif vertices[4] in v_fire # AD
                   return rotation5
             end
      elseif vertices[2] in v_fire # B firet zonder dat A firet
            if vertices[3] in v_fire # BC
                return rotation4
          elseif vertices[4] in v_fire # BD
                return rotation2
            end
      else # CD, enige overgebleven optie
            return rotation6
        end

 elseif length(v_fire) == 3 # 3 vertices firen
       if !(vertices[1] in v_fire) # BCD (A zit er niet in)
             return rotation1
      elseif !(vertices[2] in v_fire) # ACD (B zit er niet in)
            return rotation2
      elseif !(vertices[3] in v_fire) # ABD (C zit er niet in)
            return rotation3
      elseif !(vertices[4] in v_fire) # ABC (D zit er niet in)
            return rotation4
        end
    end
end


"""
    compute_NS(tetraeder::Array{Any,1},etransition::Array{Any,1},v_fire::Set{Int64})
Computes for a given tetrahedron the surface normal vector of the current wave front.
In the case 2 vertices (A and B) have fired, the returned surface normal vector is the mean of
the surface normal vectors of the triangles obtained by intersecting AC AD BD (triangle 1)
and BD AC BC (triangle 2) respectively on the points where the wave front is situated.

@param      (Array{Any,1})  `tetraeder` \n
            One row of the 'heart', containing the vertices of the tetrahedron and
            the surface normal vectors of the faces. \n
@param      (Array{Float64,1}) `etransition` \n
            Transition factors of the edges, in the order [AB BD DA CA CB CD]\n
@param      (Set{Int64}) `v_fire`\n
            Set of vertices that have fired\n

@return     (Array{Float64,2}) \n
            Surface normal vector of the wave front in the given tetrahedron. \n
"""
function compute_NS(tetraeder::Array{Any,1},etransition::Array{Float64,1},v_fire::Set{Any})
      # convert row with info tetrahedron (='tetraeder') to matrix type 'tetrahedron' (used to be the function convert_tetrahedron)
      # 1st: retrieve the surface vectors S_A, S_B, S_C, S_D and put them as rows in a (4,3)-matrix
      tetra_surf = vcat(tetraeder[6:8]',tetraeder[9:11]',tetraeder[12:14]',tetraeder[15:17]')
      # concatenate column with the indices of the edges of the tetrahedron in the orde ABCD with the rows with the corresponding surface vectors
      tetrahedron = hcat(tetraeder[2:5],tetra_surf)

      if length(v_fire)!=0 && length(v_fire)!=4 # at least one vertex fires, not all
            # get the permutation of the vertices and the edges of the tetrahedron matrix
            (indices_tetra,indices_etrans) = get_permutation(tetrahedron[:,1],v_fire)
            tetrahedron2 = zeros(4,4)
            etransition2 = zeros(6)
            # make a new tetrahedron matrix, with the rows in the 'rotated' order
            for i = 1:4
                  tetrahedron2[i,:] = tetrahedron[indices_tetra[i],:]
            end
            # make a new etransition array, with the elements in the 'rotated' order
            for i = 1:6
                  etransition2[i] = etransition[indices_etrans[i]]
            end

            # Now tetrahedron2 and etransition2 matrices contain the info of the original tetrahedron,
            # but the vertices are renamed such that A, A and B or B, C and D are firing.
            # We distinguish 3 cases in computing the surface normal vector of the wavefront in the tetrahedron
            if length(v_fire)==1
                  return geval1(tetrahedron2,etransition2)
            elseif length(v_fire)==2
                  return geval2(tetrahedron2,etransition2)
            elseif length(v_fire)==3
                  return geval3(tetrahedron2,etransition2)
            end
      else
            return [0,0,0] # geen enkele vertex firet of golven doven elkaar uit op edges
      end
end


"""
    geval1(tetrahedron::Array{Any,2},etransition::Array{Float64,1})
Computes the surface normal vector of the wave front in the case that one vertex,
A, has fired.

@param      (Array{Float64,2}) `tetrahedron`\n
            Contains the info on the tetrahedron: each row starts with the vertex
            number, followed by the components of the surface normal vector of the
            face opposite to the vertex. The rows are ordered [A B C D] \n
@param      (Array{Float64,1}) `etransition` \n
            Transition factors of the edges, in the order [AB BD DA CA CB CD] \n

@return     (Array{Float64,1}) \n
            Surface normal vector of the wave front in the given tetrahedron.
"""
function geval1(tetrahedron::Array{Float64,2},etransition::Array{Float64,1})
      #we gaan op zoek naar α_AC, α_AD en α_AB in de geroteerde lijst
      #In de geroteerde lijst hebben we voor een vaste volgorde gezorgd zodat ze respectievelijk
      #op de 4de, 3de en 1ste plaats in de lijst staan.
      #nu gaan we op zoek naar de l_trans of h_trans voor elke edge

      α_AC = etransition[4]
      α_AD = etransition[3]
      α_AB = etransition[1]

      #de oppervlaktevectoren worden bijgehouden in tetrahedron.
      #hier hebben we nood aan S_B, S_C en S_D en die staan op de respectievelijk 2de, 3de en 4de plaats
      #de elementen van de vectoren bevinden zich op plaatsen 2 tot 4 in de overeenkomstige kolom
      S_B= tetrahedron[2,2:4]
      S_C= tetrahedron[3,2:4]
      S_D= tetrahedron[4,2:4]

      NS= -α_AC*α_AD*S_B-α_AB*α_AD*S_C-α_AC*α_AB*S_D #we berekenen NS volgens de gevonden formules

      return NS
end

"""
      geval2(tetrahedron::Array{Any,2},etransition::Array{Float64,1})
Computes the surface normal vector of the wave front in the case that two vertices,
A and B, have fired. The returned surface normal vector is the mean of
the surface normal vectors of the triangles obtained by intersecting AC AD BD (triangle 1)
and BD AC BC (triangle 2) respectively on the points where the wave front is situated.

@param      (Array{Float64,2}) `tetrahedron`\n
            Contains the info on the tetrahedron: each row starts with the vertex
            number, followed by the components of the surface normal vector of the
            face opposite to the vertex. The rows are ordered [A B C D] \n
@param      (Array{Float64,1}) `etransition` \n
            Transition factors of the edges, in the order [AB BD DA CA CB CD] \n

@return     (Array{Float64,1}) \n
            Surface normal vector of the wave front in the given tetrahedron.

"""

function geval2(tetrahedron::Array{Float64,2},etransition::Array{Float64,1})
      #we gaan op zoek naar α_AC, α_BD, α_BC en α_AD in de geroteerde lijst
      #In de geroteerde lijst hebben we voor een vaste volgorde gezorgd zodat ze respectievelijk
      #op de 4de, 2de, 5de en 3de plaats in de lijst staan.
      #nu gaan we op zoek naar de l_trans of h_trans voor elke edge

      α_AC = etransition[4]
      α_AD = etransition[3]
      α_BC = etransition[5]
      α_BD = etransition[2]

      #de oppervlaktevectoren worden bijgehouden in tetrahedron.
      #hier hebben we nood aan S_A, S_B, S_C en S_D en die staan op de respectievelijk 1ste, 2de, 3de en 4de plaats
      #de elementen van de vectoren bevinden zich op plaatsen 2 tot 4 in de overeenkomstige kolom

      S_A= tetrahedron[1,2:4]
      S_B= tetrahedron[2,2:4]
      S_C= tetrahedron[3,2:4]
      S_D= tetrahedron[4,2:4]

      #hier berekenen we nu NS volgens de gevonden formules
      NS= (-α_BD*α_BC)*S_A-(α_AC*α_AD)*S_B+(α_AC*α_BC-α_BC-α_AC)*S_D+(α_AD*α_BD-α_AD-α_BD)*S_C
      return NS
end

"""
      geval3(tetrahedron::Array{Any,2},etransition::Array{Float64,1})
Computes the surface normal vector of the wave front in the case that three vertices,
B,C and D, have fired.

@param      (Array{Float64,2}) `tetrahedron`\n
            Contains the info on the tetrahedron: each row starts with the vertex
            number, followed by the components of the surface normal vector of the
            face opposite to the vertex. The rows are ordered [A B C D] \n
@param      (Array{Float64,1}) `etransition` \n
            Transition factors of the edges, in the order [AB BD DA CA CB CD] \n

@return     (Array{Float64,1}) \n
            Surface normal vector of the wave front in the given tetrahedron.
"""
function geval3(tetrahedron::Array{Float64,2},etransition::Array{Float64,1})
      #in dit geval hebben we verondersteld dat alleen A NIET firet, hierdoor is de volgorde in de benaming van de edges omgedraaid.
      #we gaan op zoek naar α_CA, α_DA en α_BA  in de geroteerde lijst
      #In de geroteerde lijst hebben we voor een vaste volgorde gezorgd zodat ze respectievelijk
      #op de 4de, 3de en 1ste plaats in de lijst staan.
      #nu gaan we op zoek naar de l_trans of h_trans voor elke edge

      α_CA = etransition[4]
      α_DA = etransition[3]
      α_BA = etransition[1]

      #de oppervlaktevectoren worden bijgehouden in tetrahedron.
      #hier hebben we nood aan S_B, S_C en S_D en die staan op de respectievelijk 2de, 3de en 4de plaats
      #de elementen van de vectoren bevinden zich op plaatsen 2 tot 4 in de overeenkomstige kolom

      S_B= tetrahedron[2,2:4]
      S_C= tetrahedron[3,2:4]
      S_D= tetrahedron[4,2:4]


      #hier berekenen we NS aan de hand van de gevonden formules
      NS= (1-α_BA)*(1-α_DA)*S_C+(1-α_CA)*(1-α_BA)*S_D+(1-α_CA)*(1-α_DA)*S_B

      return NS
end


"""
    get_potential!(ECGstruct::ECG,celAutom::CellulaireAutomaat)
Iterates over all tetrahedrons in the heart and computes potential in the electrodes
contained in ECGstruct.

@param      (ECG) `ECGstruct`\n
            An object of the class ECG that contains the number of electrodes, as wel
            as a dictionary with the electrodes as keys and the potential at the
            electrode for each time step in an array as values. It carries also
            an array with the area of the wavefront at each timestep.
@param      (CellulaireAutomaat) `celAutom` \n
            An object of the class CellulaireAutomaat that carries an object of
            type MetaGraph representing the heart, as wel as the info on the tetrahedrons,
            the active edges along with the transition fractions on the edges and
            a dictionary with the column indices of the connecting projection vector (values)
            per electrode (keys) \n

@post       ECGstruct.wavefront is updated with the area of the wavefront at the
            present timestep.
@post       The values of all keys (electrodes) of ECGstruct.elec is updated with
            the potential at the given electrode at te present timestep.
@post       An expected potential in the future is added to represent the repolarisation
            potential. The potential at the present timestep in an electrode is
            spread out in a normal distribution way, with as mean celAutom.time+APD
            APD is the mean of the Action Potential Duration of the vertices
            that have fired in a given tetrahedron. The standard deviation is
            calculated such that approximately 98% of the area under the normal
            distribution curve is within n_spread from the mean.
            One can retrieve the electrode potential at the present time celAutom.time
            for a given electrode by calling ECGstruct.elec[electrode][celAutom.time].

"""
function get_potential!(ECGstruct::ECG,celAutom::CellulaireAutomaat)
      name_elec = collect(keys(ECGstruct.elec))
      n_elec = size(name_elec,1)
      #uitkomst = zeros(n_elec)
      wavefront = 0

      indices_elec = celAutom.elec # dictionary containing per electrode (key) the indices of the columns that contain the connecting projection vector for the particular electrode
      heart = celAutom.heart
      mg = celAutom.mg
      p_0 = ECGstruct.p_0
      p = ECGstruct.p
      n_spread = Int((length(p)-1)/2)

      for i in range(1,stop=size(heart,1)) # iterate over rows 'heart' = 'tetraeder'
            tetraeder = heart[i,:] # info tetrahedron in an array
            vertixen =round.(Int, tetraeder[1:5]) # turn numbers of vertices into integers
            tetraeder = vcat(vertixen,tetraeder[6:end]) # construct the tetraeder array again
            (v_fire,etransition) = find_v_fire(vertixen[2:end],celAutom.edgesA)

            if length(v_fire) != 0
                  NS=compute_NS(tetraeder,etransition,v_fire) # compute surface normal vector of the wavefront
                  wavefront += norm(NS)
                  APD = 0 # calculated as the mean of the APD's of the vertices that fire
                  for vertex in v_fire
                        APD += get_prop(celAutom.mg,vertex,:APD)
                  end
                  APD /= length(v_fire)
                  APD = floor(Int,APD) # convert APD into time units -> not necessary

                  for electrode in name_elec # for each given electrode
                        # dot product of surface normal vector and connecting projection vector for the particular electrode, stored in the tetraeder array
                        # hier zou Threads evt een probleem kunnen vormen
                        uitkomst = dot(NS,tetraeder[indices_elec[electrode]])
                        ECGstruct.elec[electrode][celAutom.time] += p_0*uitkomst

                        # T-wave
                        for i = 1:2*n_spread+1
                              if celAutom.time+APD+i-n_spread-1 <= ECGstruct.n_calcs # in order not to exceed boundary
                                    ECGstruct.elec[electrode][celAutom.time+APD+i-n_spread-1] -= p[i]*uitkomst
                                    # celAutom.time+APD+i-n_spread-1: go APD indices further than celAutom.time, then i-n_spread-1 is the deviation of the center celAutom.time+APD
                                    # see 'normalDistributionTWave.pdf' on github
                              end
                        end
                  end
            end
      end

      ECGstruct.wavefront[celAutom.time] = wavefront
      #display(ECGstruct.elec["RA"][150:160])
end


"""
      ecg3!(ECGstruct::ECG,celAutom::CellulaireAutomaat)
calculates the 3-lead potentials and stores them in ECGstruct.lead

@param      (ECG) `ECGstruct`\n
            An object of the class ECG that contains the number of electrodes, as wel
            as a dictionary with the electrodes as keys and the potential at the
            electrode for each time step in an array as values. It carries also
            an array with the area of the wavefront at each timestep.
@param      (CellulaireAutomaat) `celAutom` \n
            An object of the class CellulaireAutomaat that carries an object of
            type MetaGraph representing the heart, as wel as the info on the tetrahedrons,
            the active edges along with the transition fractions on the edges and
            a dictionary with the column indices of the connecting projection vector (values)
            per electrode (keys) \n

@post       ECGstruct.wavefront is updated with the area of the wavefront at the
            present timestep.
@post       The values of all keys (electrodes) of ECGstruct.elec is updated with
            the potential at the given electrode at te present timestep.
@post       An expected potential in the future is added to represent the repolarisation
            potential. See `get_potential!` for more info on the way this is calculated.
            One can retrieve the electrode potential at the present time celAutom.time
            for a given electrode by calling ECGstruct.elec[electrode][celAutom.time].
@post       ECGstruct.lead contains for every lead in the 3-lead system the corresponding
            potential at the present time step. One can retrieve the lead potential
            at the present time celAutom.time for a given lead by calling ECGstruct.lead[lead][celAutom.time].

"""
function ecg3!(ECGstruct::ECG,celAutom::CellulaireAutomaat)
      get_potential!(ECGstruct,celAutom)

      # compute lead potentials with given formulas
      ECGstruct.lead["I"]=ECGstruct.elec["LA"]-ECGstruct.elec["RA"]
      ECGstruct.lead["II"]=ECGstruct.elec["LL"]-ECGstruct.elec["RA"]
      ECGstruct.lead["III"]=ECGstruct.elec["LL"]-ECGstruct.elec["LA"]
end


"""
      ecg_beam!(celAutom::CellulaireAutomaat)
This function is especially written to calculate the EG for a beam.

@param      (ECG) `ECGstruct`\n
            An object of the class ECG that contains the number of electrodes, as wel
            as a dictionary with the electrodes as keys and the potential at the
            electrode for each time step in an array as values. It carries also
            an array with the area of the wavefront at each timestep.
@param      (CellulaireAutomaat) `celAutom` \n
            An object of the class CellulaireAutomaat that carries an object of
            type MetaGraph representing the beam, as wel as the info on the tetrahedrons,
            the active edges along with the transition fractions on the edges and
            a dictionary with the column indices of the connecting projection vector (values)
            per electrode (keys) \n

@post       ECGstruct.wavefront is updated with the area of the wavefront at the
            present timestep.
@post       The values of all keys (electrodes) of ECGstruct.elec is updated with
            the potential at the given electrode at te present timestep.
@post       An expected potential in the future is added to represent the repolarisation
            potential. See `get_potential!` for more info on the way this is calculated.
            One can retrieve the electrode potential at the present time celAutom.time
            for a given electrode by calling ECGstruct.elec[electrode][celAutom.time].
@post       ECGstruct.lead contains the lead "verschil" which is defined as
            the difference between the potentials in electrode 1 and 2 as well as
            the corresponding potential at the present time step. One can retrieve
            the lead potential at the present time celAutom.time for a given lead by calling ECGstruct.lead[lead][celAutom.time].
"""
function ecg_beam!(ECGstruct::ECG,celAutom::CellulaireAutomaat)
      get_potential!(ECGstruct,celAutom) # array with potentials in each electrode
      ECGstruct.lead["verschil"] = ECGstruct.elec["2"] - ECGstruct.elec["1"]
end


"""
      function ecg12!(ECGstruct::ECG,celAutom::CellulaireAutomaat)
calculates the 12-lead potentials and stores them in ECGstruct.lead

@param      (ECG) `ECGstruct`\n
            An object of the class ECG that contains the number of electrodes, as wel
            as a dictionary with the electrodes as keys and the potential at the
            electrode for each time step in an array as values. It carries also
            an array with the area of the wavefront at each timestep.
@param      (CellulaireAutomaat) `celAutom` \n
            An object of the class CellulaireAutomaat that carries an object of
            type MetaGraph representing the beam, as wel as the info on the tetrahedrons,
            the active edges along with the transition fractions on the edges and
            a dictionary with the column indices of the connecting projection vector (values)
            per electrode (keys) \n

@post       ECGstruct.wavefront is updated with the area of the wavefront at the
            present timestep.
@post       The values of all keys (electrodes) of ECGstruct.elec is updated with
            the potential at the given electrode at te present timestep.
@post       An expected potential in the future is added to represent the repolarisation
            potential. See `get_potential!` for more info on the way this is calculated.
            One can retrieve the electrode potential at the present time celAutom.time
            for a given electrode by calling ECGstruct.elec[electrode][celAutom.time].
@post       ECGstruct.lead contains for every lead in the 3-lead system the corresponding
            potential at the present time step. One can retrieve the lead potential
            at the present time celAutom.time for a given lead by calling ECGstruct.lead[lead][celAutom.time].

"""
function ecg12!(ECGstruct::ECG,celAutom::CellulaireAutomaat)
      get_potential!(ECGstruct,celAutom)
      # Wilkinson Central Terminal
      # this is stored with the leads in order not to have complications with get_potential!
      ECGstruct.lead["WCT"]=(ECGstruct.elec["LA"]+ECGstruct.elec["LL"]+ECGstruct.elec["RA"])/3
      # compute lead potentials with given formulas
      ECGstruct.lead["I"]=ECGstruct.elec["LA"]-ECGstruct.elec["RA"]
      ECGstruct.lead["II"]=ECGstruct.elec["LL"]-ECGstruct.elec["RA"]
      ECGstruct.lead["III"]=ECGstruct.elec["LL"]-ECGstruct.elec["LA"]
      ECGstruct.lead["aVL"]=ECGstruct.elec["LA"]-0.5*ECGstruct.elec["RA"]-0.5*ECGstruct.elec["LL"]
      ECGstruct.lead["aVR"]=ECGstruct.elec["RA"]-0.5*ECGstruct.elec["LA"]-0.5*ECGstruct.elec["LL"]
      ECGstruct.lead["aVF"]=ECGstruct.elec["LL"]-0.5*ECGstruct.elec["RA"]-0.5*ECGstruct.elec["LA"]
      leads = ["V1","V2","V3","V4","V5","V6"]
      for lead in leads
            ECGstruct.lead[lead]=ECGstruct.elec[lead]-ECGstruct.lead["WCT"]
      end
end


"""
"""
function ecg263!(ECGstruct::ECG,celAutom::CellulaireAutomaat)
      get_potential!(ECGstruct,celAutom)
end


"""
    initializeECG(type::String,amountCalcs::Int64)
Makes an array in which the data on the ECG is stored. The number of columns
is equal to the number of leads. Each row wil contain the value of the ECG in a
certain time step. The difference in time between the rows is celAutom.δt.

@param      (String) `type`
            "beam"      to calculate the EG in a beam
            "heart3"    to calculate a 3-lead ECG in a heart
            "heart12"   to calculate a 12-lead ECG in a heart
@param      (Int64) `amountCalcs`
            number of calculations the simulation consists of

@return     (Array{Int64,2}) `elec`
            Array with dimensions (amountCalcs,numberOfLeads), filled with zeros.
"""
function initializeECG(t::String,amountCalcs::Int64,celAutom::CellulaireAutomaat)
      elec_init = Dict{String,Array{Float64,1}}()
      lead_init = Dict{String,Array{Float64,1}}()
      loc_elec = zeros(Float64,3,1)

      if t == "beam"
            name_elecs = ["1","2"]
            name_leads=["verschil"]
      elseif t == "heart3"
            name_elecs = ["RA","LA","LL"]
            name_leads = ["I","II","III"]
      elseif t == "heart12"
            name_elecs = ["RA","LA","LL","V1","V2","V3","V4","V5","V6"]
            name_leads = ["I","II","III","aVL","aVR","aVF","V1","V2","V3","V4","V5","V6"]
      elseif t == "heart263"
            (loc_elec,name_elecs) = readdlm("electrode_pos.dat",',',header=true)
            name_elecs = reshape(collect(String,name_elecs),263)
            name_leads = cat(name_elecs[1:260],["I","II","III"],dims=1) # exclude LL, LA and RA
      end

      for i = 1:length(name_elecs)
            elec_init[name_elecs[i]] = float(zeros(amountCalcs))
      end

      for lead in name_leads
            lead_init[lead]=float(zeros(amountCalcs))
      end

      # this will be needed in get_potential!, but only needs to be calculated once
      # Make sure there are no other variables called p or p_0.
      p_0 = 1.0 # wild guess
      println("p_0 = ",p_0)
      # calculate p, which contains the normal distribution spread of the repolarisation wave
      σ = 50/celAutom.δt # see paper, unit = t.u.
      # number of timesteps on each side of t+APD that will be used to calculate the spread of the T-wave
      # suggested to take 3σ
      n_spread = floor(Int64,3*σ)
      d = Normal(0,σ)
      p = repeat([Float64(p_0)],2*n_spread+1)
      for i = 1:2*n_spread+1
            p[i] *= pdf(d,i-n_spread-1)
      end
      #println("p = ",p,"\nsom = ",sum(p),"\nlengte = ",length(p))


      return ECG((type=t,elec=elec_init,lead=lead_init,loc_elec=loc_elec,name_elec=name_elecs,wavefront=float(zeros(amountCalcs)),n_calcs=amountCalcs,p_0=p_0,p=p)...)
end


"""
"""
function updateECG!(ECGstruct::ECG,celAutom::CellulaireAutomaat)
      if ECGstruct.type == "beam"
            ecg_beam!(ECGstruct,celAutom)
      elseif ECGstruct.type == "heart3"
            ecg3!(ECGstruct,celAutom)
      elseif ECGstruct.type == "heart12"
            ecg12!(ECGstruct,celAutom)
      elseif ECGstruct.type == "heart263"
            ecg263!(ECGstruct,celAutom)
      end
end


"""
"""
function plotECG(ECGstruct::ECG,celAutom::CellulaireAutomaat,x_ax::Array{Int64,1},i::Int64)
      if ECGstruct.type == "beam"
            EGleads = hcat(ECGstruct.elec["1"],ECGstruct.elec["2"],ECGstruct.lead["verschil"])
            ecg_plot = plot(x_ax[1:i],EGleads[1:i,:],
                              title="EG beam",xlabel="Time (ms)",ylabel="Voltage (A.U.)",
                              xlims=(0,ECGstruct.n_calcs*celAutom.δt),ylims=(-3,3),
                              label=["1" "2" "Verschil"])
            display(ecg_plot)
            return ecg_plot
      elseif ECGstruct.type == "heart3"
            ECGleads = hcat(ECGstruct.lead["I"],ECGstruct.lead["II"],ECGstruct.lead["III"])
            ecg_plot = plot(x_ax[1:i],ECGleads[1:i,:],
                      title="ECG heart",xlabel="Time (ms)",ylabel="Voltage (A.U.)",
                      xlims=(0,ECGstruct.n_calcs*celAutom.δt),ylims=(-0.2,0.3),
                      label=["I" "II" "III"])
            display(ecg_plot)
            return ecg_plot
      elseif ECGstruct.type == "heart12"
            ecg_plot = repeat([plot(1)],4)
            ecg_plot[1] = plot(x_ax[1:i],hcat(ECGstruct.lead["I"],ECGstruct.lead["II"],ECGstruct.lead["III"])[1:i,:],
                      title="12-lead ECG heart (part 1)",xlabel="Time (ms)",ylabel="Voltage (A.U.)",
                      xlims=(0,ECGstruct.n_calcs*celAutom.δt),ylims=(-0.2,0.2),
                      label=["I" "II" "III"])
            ecg_plot[2] = plot(x_ax[1:i],hcat(ECGstruct.lead["aVL"],ECGstruct.lead["aVR"],ECGstruct.lead["aVF"])[1:i,:],
                      title="12-lead ECG heart (part 2)",xlabel="Time (ms)",ylabel="Voltage (A.U.)",
                      xlims=(0,ECGstruct.n_calcs*celAutom.δt),ylims=(-0.2,0.2),
                      label=["aVL" "aVR" "aVF"])
            ecg_plot[3] = plot(x_ax[1:i],hcat(ECGstruct.lead["V1"],ECGstruct.lead["V2"],ECGstruct.lead["V3"])[1:i,:],
                      title="12-lead ECG heart (part 3)",xlabel="Time (ms)",ylabel="Voltage (A.U.)",
                      xlims=(0,ECGstruct.n_calcs*celAutom.δt),ylims=(-0.8,0.6),
                      label=["V1" "V2" "V3"])
            ecg_plot[4] = plot(x_ax[1:i],hcat(ECGstruct.lead["V4"],ECGstruct.lead["V5"],ECGstruct.lead["V6"])[1:i,:],
                      title="12-lead ECG heart (part 4)",xlabel="Time (ms)",ylabel="Voltage (A.U.)",
                      xlims=(0,ECGstruct.n_calcs*celAutom.δt),ylims=(-0.4,0.4),
                      label=["V4" "V5" "V6"])
            display(ecg_plot[1])
            # we chose to visualise only the 3-lead ecg in runtime, for displaying all 12
            # lead potentials in one plot would not be clear.
            return ecg_plot
      elseif ECGstruct.type == "heart263"
            c_elec = zeros(length(ECGstruct.name_elecs))
            for i in range(1,length(ECGstruct.name_elecs))
                  c_elec[i] = ECGstruct.elec[ECGstruct.name_elecs[i]][celAutom.time]
            end
            x = ECGstruct.loc_elecs[1,:]
            y = ECGstruct.loc_elecs[2,:]
            z = ECGstruct.loc_elecs[3,:]
            plot_pot = scatter(x=x,y=y,z=z,mode="markers",
                        marker=attr(cmin=minimum(c_elec), cmax=maximum(c_elec), color=c_elec, colorscale="Bluered"))
            display(plot_pot)
            tu = celAutom.time
            png(plot_pot,"groepje3/heart263/body_surface$tu.png")
      end
end
