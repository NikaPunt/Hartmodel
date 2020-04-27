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
      ECG_calcs::Dict{String,Array{Float64,1}}
      n_calcs::Int64
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
    get_potential!(celAutom::CellulaireAutomaat,name_elec::Array{String,1})
Iterates over all tetrahedrons in the heart and computes potential in the electrodes
in name_elec.

@param      (CellulaireAutomaat) `celAutom` \n
            An object of the class CellulaireAutomaat that carries an object of
            type MetaGraph representing the heart, as wel as the info on the tetrahedrons,
            the active edges along with the transition fractions on the edges and
            a dictionary with the column indices of the connecting projection vector (values)
            per electrode (keys) \n
@param      (Array{String,1}) `name_elec` \n
            Contains all names of the electrodes in which the potential must be
            computed. The order of the electrodes in `name_elec` will be the order
            of the potentials in the outut.\n

@return     (Array{Float64,1})\n
            The potentials in the electrodes, in the order of `name_elec`.
"""
function get_potential!(ECGstruct::ECG,celAutom::CellulaireAutomaat)
      name_elec = collect(keys(ECGstruct.ECG_calcs))
      n_elec = size(name_elec,1)
      #uitkomst = zeros(n_elec)
      wavefront = 0

      p_0= 1 # wild guess
      # calculate p, which contains the normal distribution spread of the repolarisation wave
      n_spread = 2 # number of timesteps on each side of t+APD that will be used to calculate the spread of the T-wave
      #σ = (n_spread*celAutom.δt)/2 # see 'normalDistributionTWave.pdf' on github
      d = Normal()
      p = repeat([Float64(p_0)],2*n_spread+1)
      for i = 1:2*n_spread+1
            p[i] *= pdf(d,(i-n_spread-1)*celAutom.δt)
      end


      indices_elec = celAutom.elec # dictionary containing per electrode (key) the indices of the columns that contain the connecting projection vector for the particular electrode
      heart = celAutom.heart
      mg = celAutom.mg

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
                  APD = Int(floor(APD)) # convert APD into time units -> not necessary

                  for electrode in name_elec # for each given electrode
                        # dot product of surface normal vector and connecting projection vector for the particular electrode, stored in the tetraeder array
                        # hier zou Threads evt een probleem kunnen vormen
                        uitkomst = dot(NS,tetraeder[indices_elec[electrode]])
                        ECGstruct.ECG_calcs[electrode][celAutom.time] += p_0*uitkomst

                        for i = 1:2*n_spread+1
                              if celAutom.time+APD+i-n_spread-1 <= ECGstruct.n_calcs # in order not to exceed boundary
                                    ECGstruct.ECG_calcs[electrode][celAutom.time+APD+i-n_spread-1] -= p[i]*uitkomst
                                    # celAutom.time+APD+i-n_spread-1: go APD indices further than celAutom.time, then i-n_spread-1 is the deviation of the center celAutom.time+APD
                                    # see 'normalDistributionTWave.pdf' on github
                              end
                        end
                  end
            end
      end

      return wavefront
end


"""
      ecg3(celAutom::CellulaireAutomaat)
The function returns the 3 voltages for the ecg

@param      (CellulaireAutomaat) `celAutom` \n
            An object of the class CellulaireAutomaat that carries an object of
            type MetaGraph representing the heart, as wel as the info on the tetrahedrons,
            the active edges along with the transition fractions on the edges and
            a dictionary with the column indices of the connecting projection vector (values)
            per electrode (keys) \n

@return     (Array{Float64,1})\n
            The voltages in the order E1, E2, E3.

"""
function ecg3(celAutom::CellulaireAutomaat)
      electrodes = Array{String,1}
      electrodes = ["VR","VL","VF"]
      (pot_elec,wavefront) = get_potential(celAutom,electrodes) # array with potentials in each electrode
      # compute lead potentials with given formulas
      E1= pot_elec[2]-pot_elec[1]
      E2= pot_elec[3]-pot_elec[1]
      E3= pot_elec[3]-pot_elec[2]
      return ([E1,E2,E3],wavefront)

end


"""
      ecg_beam!(celAutom::CellulaireAutomaat)
This function is especially written to calculate the EG for a beam.

@param      (CellulaireAutomaat) `celAutom` \n
            An object of the class CellulaireAutomaat that carries an object of
            type MetaGraph representing the beam, as wel as the info on the tetrahedrons,
            the active edges along with the transition fractions on the edges and
            a dictionary with the column indices of the connecting projection vector (values)
            per electrode (keys) \n

@return     (Array{Float64,1})
            The potentials in the electrodes '1' and '2'
"""
function ecg_beam!(ECGstruct::ECG,celAutom::CellulaireAutomaat)
      wavefront = get_potential!(ECGstruct,celAutom) # array with potentials in each electrode

      return wavefront
end


"""
      function ecg12(heart,mg)
The function returns the 12 voltages that make up a 12-lead ECG, in the order
[E1 E2 E3 aVL aVR aVF V1 V2 V3 V4 V5 V6]

@param      (CellulaireAutomaat) `celAutom` \n
            An object of the class CellulaireAutomaat that carries an object of
            type MetaGraph representing the beam, as wel as the info on the tetrahedrons,
            the active edges along with the transition fractions on the edges and
            a dictionary with the column indices of the connecting projection vector (values)
            per electrode (keys) \n

@return     (Array{Float64,1})
            Contains the potential in the leads that make up a 12-lead ECG, in the
            order [E1 E2 E3 aVL aVR aVF V1 V2 V3 V4 V5 V6]

"""
function ecg12(celAutom::CellulaireAutomaat)
      electrodes = Array{String,1}
      electrodes = ["VR","VL","VF","V1","V2","V3","V4","V5","V6"]
      pot_elec = get_potential(celAutom,electrodes) # array with potentials in each electrode
      # compute lead potentials
      WCT = sum(pot_elec[1:3])
      ECG = zeros(12)
      ECG[1]= pot_elec[2]-pot_elec[1] # E1
      ECG[2]= pot_elec[3]-pot_elec[1] # E2
      ECG[3]= pot_elec[3]-pot_elec[2] # E3
      ECG[4]= pot_elec[2]-0.5*pot_elec[1]-0.5*pot_elec[3] # aVL
      ECG[5]= pot_elec[1]-0.5*pot_elec[2]-0.5*pot_elec[3] # aVR
      ECG[6]= pot_elec[3]-0.5*pot_elec[1]-0.5*pot_elec[2] # aVF
      for i = 7:12 # V1 V2 V3 V4 V5 V6
            ECG[i] = pot_elec[i-3]-WCT
      end
      return ECG
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

@return     (Array{Int64,2}) `ECG_calcs`
            Array with dimensions (amountCalcs,numberOfLeads), filled with zeros.
"""
function initializeECG(t::String,amountCalcs::Int64)
      z = zeros(amountCalcs)
      ECG_calcs_init = Dict{String,Array{Float64,1}}()

      if t == "beam"
            name_elecs = ["1","2"]
      elseif t == "heart3"
            name_elecs = ["VR","VL","VF"]
      elseif t == "heart12"
            name_elecs = ["VR","VL","VF","V1","V2","V3","V4","V5","V6"]
      end

      for name in name_elecs
            ECG_calcs_init[name] = z
      end

      return ECG((type=t,ECG_calcs=ECG_calcs_init,n_calcs=amountCalcs)...)
end


"""
"""
function updateECG!(ECGstruct::ECG,celAutom::CellulaireAutomaat)
      if ECGstruct.type == "beam"
            return ecg_beam!(ECGstruct,celAutom)
      elseif ECGstruct.type == "heart3"
            # doe iets anders
      elseif ECGstruct.type == "heart12"
            # doe nog iets anders
      end

end
