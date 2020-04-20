# Using functions from 'RandomPerm.gi' and 'Orbit.gi'
Read("RandomPerm.gi"); Read("Orbit.gi");

#######################################################################
# Description: A function to calculate the pointwise image of a list. #
# Input: - list of integers B = [b_1,...,b_k]                         #   
#        - permutation x as cycle                                     #
# Output: - list of integers [b_1^x,...,b_k^x]                        #
#######################################################################

BaseImage := function(B, x)
    local i, image;
    image := [];
    for i in [1..Length(B)] do
        Add(image, B[i]^x);
    od;
    return image;
end;

#######################################################################
#######################################################################

#######################################################################
# Description: A function to setup partial generating sets of         #
#              stabilizeres.                                          #
# Input: - list of integers B = [b_1,...,b_k]                         #   
#        - generating system S as list of permutations                #
# Output: - list of generating system (lists of permutations),        #
#           one generating system for each stabilizer in chain        #
#######################################################################

GeneratingSets := function(B, S)
    local s, x, i, partialBase, sets;
    sets := [S];   
    for i in [1..Length(B)-1] do
        Add(sets, []);
    od;
    for x in S do
        for i in [1..Length(B)-1] do
            partialBase := B{[1..i]};
            if BaseImage(partialBase, x) = partialBase then
                Add(sets[i+1], x);
            else
                break;
            fi;
        od;
    od;
    return sets;
end;

#######################################################################
#######################################################################

#######################################################################
# Description: A function to strip a permutation w.r.t. a base,       #
#              described in section 3.2.                              #    
# Input: - permutation g as cycle                                     #
#        - list of integers B = [b_1,...,b_k]                         #   
#        - generating system S as list of permutations                #
#        - data structure returned by 'GeneratingSets' function       #                                                
# Output: - residue and dropout level of g                            #
#######################################################################

Strip := function(g, B, orbits, sgs)
    local m, i, u, b;
    m := g;
    for i in [1..Length(B)] do
        b := B[i]^m;
        if orbits[i][2][b] = 0 then 
            return [m, i]; 
        fi;
        u := Transversal(b, orbits[i][2], sgs[i]);  m := m*(u^(-1));                
    od;
    return [m, Length(B)+1];
end;

#######################################################################
#######################################################################
	
#######################################################################
# Description: The standard deterministic Schreier-Sims-Algorithm     #
#              described in section 3.3.                              #
# Input: - list B of base points (might be empty list [] first)       #        
#        - list S of given generators of group G                      #
#        - integer n specifying the degree of symmetric group the     #
#          generated group G is embedded in                           #
# Output: - list 'BSGS' with three entries including                  #
#           -- the calculated base B as list in BSGS[1]               #
#           -- a strong generating set S as list of partial           #
#              generating sets for each stabilizer in BSGS[2]         #
#           -- a list of the sizes of the orbits belonging to B       #
#              and S in BSGS[3]                                       #
#######################################################################

SchreierSims := function(B, S, n)
    local x, u, y, h, b, i, j, k, t, orbits, sgs, sizes, stripped;
        
    # extending partial base so that no point of S  
    # fixes all points
    for x in S do
        if BaseImage(B, x) = B then
            Add(B, MovedPoints(x)[1]);
        fi;
    od; 
    # calculating generating sets for each stabilizer 
    sgs := GeneratingSets(B, S);
        
    # calculate orbits and schreier vectors of base 
    # points under action of generated stabilizers        
    orbits := [];
    for i in [1..Length(B)] do
        Add(orbits, OrbitSV(B[i], sgs[i], n));
    od;
    i := Length(B);
    while i >= 1 do
        for b in orbits[i][1] do
            t := Transversal(b, orbits[i][2], sgs[i]);            
            for x in sgs[i] do
                u := Transversal(b^x, orbits[i][2], sgs[i]);  y := true;
                if not t*x = u then
                    # check whether Schreier generator 
                    # belongs to H^(i+1)
                    stripped := Strip(t*x*(u^(-1)), B, orbits, sgs);
                    j := stripped[2];  h := stripped[1];
                    
                    # in case the Schreier generator does  
                    # not belong to H^(i+1) do ...
                    if j <= Length(B) then
                        y := false;
                    elif not h = () then                        
                        Add(B, MovedPoints(h)[1]);
                        Add(sgs, []);  y := false;
                    fi;
                    if y = false then
                        # update generating sets, orbits 
                        # and Schreier vectors
                        for k in [i+1..j] do
                            Add(sgs[k], h);
                            orbits[k] := OrbitSV(B[k], sgs[k], n);                           
                        od;                
                        i := j;  break;
                    fi;
                fi;
            od;
            if y = false then
                break;
            fi;
        od;  
        if y = false then
            continue;
        fi;               
        i := i - 1;     
    od; 
    sizes := [];
    for i in [1..Length(B)] do
       Add(sizes, Length(orbits[i][1]));
    od;
    return [B, sgs, sizes];
end;

#######################################################################
#######################################################################

#######################################################################
# Description: The randomized Schreier-Sims-Algorithm,                #
#              described in section 4.2.                              #      
# Input: - list B of base points (might be empty list [] first)       #
#        - list S of given generators of group G                      #
#        - an integer 'bound' capping the number of elements          #
#          getting stripped                                           #
#          with trivial residue in a row                              #
#        - integer n specifying the degree of symmetric group the     #
#          generated group G is embedded in                           #
# Output: - list 'BSGS' with three entries including                  #    
#           -- the calculated base B as list in BSGS[1]               #
#           -- a strong generating set S as list of partial           #
#              generating sets for each stabilizer in BSGS[2]         #
#           -- a list of the sizes of the orbits belonging to B       #
#              and S in BSGS[3]                                       #
#######################################################################

RandomSchreierSims := function(B, S, bound, n)
    local x, g, i, h, X, c, j, k, y, sgs, orbits, sizes, SGS, stripped;
    
    for x in S do
        if BaseImage(B, x) = B then
            Add(B, MovedPoints(x)[1]);
        fi;
    od;
    
    sgs := GeneratingSets(B, S);  orbits := [];  SGS := ShallowCopy(S);
   
    for i in [1..Length(B)] do
        Add(orbits, OrbitSV(B[i], sgs[i], n));
    od;
    
   X := InitGenerator(S);  c := 0;
   while c < bound do
       g := RandomPerm(S, X);
       stripped := Strip(g, B, orbits, sgs);
       j := stripped[2];  h := stripped[1];  y := true;
       if j <= Length(B) then
           y := false;
       elif not h = () then
           y := false;  Add(B, MovedPoints(h)[1]);  Add(sgs, []);
       fi;
       if y = false then
           Add(SGS, h);
           for k in [2..j] do
               Add(sgs[k], h);
               orbits[k] := OrbitSV(B[k], sgs[k], n);      
           od;      
           c := 0;                
       fi;
       c := c + 1;    
   od;
   
   sizes := [];  sgs[1] := SGS;
   for i in [1..Length(B)] do
       Add(sizes, Length(orbits[i][1]));
   od;
   Add(sgs, [()]);
   return [B, sgs, sizes];
end;

#######################################################################
#######################################################################

#######################################################################
# Description: A function to remove redundant generators from         #       
#              a strong generating set, described in section 4.2.     #                             
# Input: - BSGS data structure returned by RandomSchreierSims         # 
#        - Generating set which defined the group. NOT the strong     #
#          generating set returned by RandomSchreierSims!             #
#        - integer n specifying the size of the set the group         #
#          operates on
#######################################################################

RemoveRedundantGens := function(BSGS, initialGenSet, n)
    local orbits, i, g, newGenSys, removeGen;
    orbits := [];
    
    removeGen := function(SGS, gen, level)
        local j;
        for j in [1..level] do
            RemoveSet(SGS[j], gen);
        od;
    end;    
    for i in [1..Length(BSGS[1])] do
        BSGS[2][i] := Set(BSGS[2][i]);
        Add(orbits, Set(OrbitSimple(BSGS[1][i], BSGS[2][i]))); 
    od;
    for i in Reversed([1..Length(BSGS[1])-1]) do
        for g in BSGS[2][i] do
            if not g in BSGS[2][i+1] 
               and not g in initialGenSet then
                newGenSys := Set(BSGS[2][i]); 
                RemoveSet(newGenSys, g);
                if Set(OrbitSimple(BSGS[1][i], newGenSys)) = orbits[i] 
                	then removeGen(BSGS[2], g, i);
                fi;
            fi;
        od;
    od;
    orbits := []; 
    for i in [1..Length(BSGS[1])] do
        Add(orbits, OrbitSV(BSGS[1][i], BSGS[2][i], n));
    od;
    Add(BSGS, orbits);  MakeImmutable(BSGS);
end;

#######################################################################
#######################################################################

#######################################################################
# Description: The sequential deterministic verification,             #       
#              described in section 4.3                               #                             
# Input: - BSGS data structure returned by RandomSchreierSims         # 
#        - integer m specifying which layer H^{(m)}_{b_m} = H^{(m+1)} #
#          has to be verified                                         #
# Output: - boolean, true if H^{(m)}_{b_m} = H^{(m+1)} holds,         #
#           else false                                                #
#######################################################################

SimpleVerify := function(BSGS, m)
    local x, u, b, t, orbits, gen, stripped, sgs, B;
          
    B := BSGS[1]; sgs := BSGS[2]; orbits := BSGS[4];
    
    for b in orbits[m][1] do
         t := Transversal(b, orbits[m][2], sgs[m]);
         for x in sgs[m] do
            if not x in sgs[m+1] or orbits[m+1][2][b] = 0 then
                 u := Transversal(b^x, orbits[m][2], sgs[m]);
                 if not t*x = u then 
                    gen := t*x*(u^(-1));
                    # check whether Schreier generator 'gen' 
                    # belongs to H^(i+1)
                    stripped := Strip(gen, B, orbits, sgs);                                                             
                    # if not, return                                                            
                    if stripped[2] <= Length(B) then
                        return false;
                    elif not stripped[1] = () then
                        return false;
                    fi;                    
                 fi;
              fi;
         od;
    od;                   
    return true;
end;

#######################################################################
#######################################################################

#######################################################################
# Description: The deterministic verification, parallelized w.r.t.    #       
#              the length of the basic orbits, described as method 2  #
#              in section 5.1.                                        #                             
# Input: - BSGS data structure returned by RandomSchreierSims         # 
#        - integer m specifying which layer H^{(m)}_{b_m} = H^{(m+1)} #
#          has to be verified                                         #
# Output: - boolean, true if H^{(m)}_{b_m} = H^{(m+1)} holds,         #
#           else false                                                #
#######################################################################

ParallelVerify := function(BSGS, m)
    local b, i, orbits, sgs, B, Verify, tasks, results, numOfTasks, 
          remaining, partition, pointsPerTask;
    
    # the well-known deterministic verify routine
    Verify := function(B, sgs, orbits, m, partition)
        local x, u, i, stripped, gen, t;
        for b in partition do
            t := Transversal(b, orbits[m][2], sgs[m]);
                 for x in sgs[m] do
                    if not x in sgs[m+1] or orbits[m+1][2][b] = 0 then
                         u := Transversal(b^x, orbits[m][2], sgs[m]);
                         if not t*x = u then 
                            gen := t*x*(u^(-1));
                            stripped := Strip(gen, B, orbits, sgs); 
                                                                               
                            if stripped[2] <= Length(B) then
                                return false;
                            elif not stripped[1] = () then
                                return false;
                            fi;                    
                         fi;
                      fi;
                 od;
        od;
        return true;
    end;

    B := BSGS[1]; sgs := BSGS[2]; 
    orbits := BSGS[4]; partition := [];
        
    pointsPerTask := Int(Length(orbits[m][1])*0.05) + 1;
    numOfTasks := Int(Length(orbits[m][1])/pointsPerTask) + 1;
    remaining := Length(orbits[m][1]) mod pointsPerTask;
    
    tasks := [1..numOfTasks]; results := [1..numOfTasks];
    
    # slice orbit into 'numOfTask' many nearly equally sized subsets
    for i in [1..numOfTasks-1] do
        Add(partition, orbits[m][1]{[(i-1)*pointsPerTask+1..
                                         i*pointsPerTask]});
    od;
    
    Add(partition, orbits[m][1]{[(numOfTasks-1)*pointsPerTask+1..
                    ((numOfTasks-1)*pointsPerTask+remaining)]});
    
    # start several verifications ...
    for i in [1..numOfTasks] do
        tasks[i] := RunTask(Verify, B, sgs, orbits, m, partition[i]);
    od;
    
    # ... and evaluate their results afterwards
    for i in [1..numOfTasks] do
        results[i] := TaskResult(tasks[i]);
    od;
    if false in results then
        return false;
        
    else
        return true;
    fi;   
end;

