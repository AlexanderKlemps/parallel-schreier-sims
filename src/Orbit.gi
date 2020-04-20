#######################################################################
# Description: A function to calculate orbits, described in           #
#              section 3.1.                                           #                             
# Input: - integer 0 < a <= n + 1                                     #   
#        - list of generators S                                       #
# Output: - list of integers 'orbit'                                  #
#######################################################################

OrbitSimple := function(a, S)
    local orbit, b, x;
    orbit := [a];
    for b in orbit do
        for x in S do
            if not b^x in orbit then
                Add(orbit, b^x);
            fi;
        od;
    od;
    return orbit;
end;

#######################################################################
#######################################################################

#######################################################################
# Description: A function to calculate orbits and belonging Schreier- #
#              vectors, described in section 3.1.                     #    
# Input: - integer 0 < a <= n                                         #   
#        - list of generators S                                       #
# Output: - list of integers 'orbit'                                  #
#         - belonging Schreier-vector 'sv' as list of integers        #
#######################################################################

OrbitSV := function(a, S, n)
    local orbit, b, x, i, sv;
    orbit := [a];  sv := [1..n]*0;  sv[a] := -1;
    for b in orbit do
        for i in [1..Length(S)] do
            if not b^S[i] in orbit then
                Add(orbit, b^S[i]);  sv[b^S[i]] := i;
            fi;
        od;
    od;
    return [orbit, sv];
end;

#######################################################################
#######################################################################

#######################################################################
# Description: A function to calculate transversal elements, described#
#              in section 3.1.                                        #    
# Input: - orbit element b of type integer 0 < b <= n                 #   
#        - Schreier-vector 'sv' as list of integers in range [1..n]   #
#        - list of generators S                                       #
# Output: - transversal element of b w.r.t. a and <S>                 # 
#######################################################################

Transversal := function(b, sv, S)
    local u, k;
    if sv[b] = 0 then
        return false;
    fi;
    u := ();  k := sv[b];
    while not k = -1 do
        u := S[k]*u;  b := b^(S[k]^-1);  k := sv[b];
    od;
    return u;
end; 
