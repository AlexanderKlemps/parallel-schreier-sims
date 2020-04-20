#######################################################################
# Description: A function to initialize a generator for random        #
#              permutations, described in section 4.1.                #    
# Input: - list of generators S                                       #
# Output: - list of permutations X including repetitions of generators#
#           from S                                                    #
#######################################################################

InitGenerator := function(S)
    local i, X, k;
    X := List(S); k := Length(X);
    for i in [k+1..Maximum(11, Length(S))] do
        Add(X, X[i-k]);
    od;
    if IsReadOnlyGlobal("X") = true then 
        MakeReadWriteGlobal("X");
    fi;
    return X;
end;

#######################################################################
#######################################################################

#######################################################################
# Description: A function to calculate a random permutation,          #
#              described in section 4.1.                              #    
# Input: - list of generators S, generator X returned by InitGenerator#
# Output: - random permutation a as cycle                             #
#######################################################################

RandomPerm := function(S, X)
    local Randomize, a, r, i;
    
    r := Maximum(11, Length(S));    
    Randomize := function()
        local s, t, e, range;
        
        # indices s and t specify which permutations X[s], X[t]
        # will be multiplied with each other
        # e specifies whether X[t] or its inverse will be used
        range := [1..r];  s := Random(range);  e := Random([-1,1]);
        t := Random(range);  Remove(range, s);  
                
        # alternate the order of multiplication 
        # of X[t] and X[s] randomly
        # update generator X
        
        if Random([1,2]) = 1 then
           X[s] := X[s]*X[t]^e;
           a := a*X[s];
        else
           X[s] := X[t]^e*X[s];
           a := X[s]*a;
        fi;
        return a;
    end;    
    
    # calculate random permutation a as
    # product of random permutations from X
    a := ();
    for i in [1..50] do
        Randomize();
    od;
    return a;
end;
