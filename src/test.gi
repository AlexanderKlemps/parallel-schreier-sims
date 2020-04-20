Read("SchreierSims.gi");
LoadPackage( "AtlasRep" );

MicroSeconds := function()
  local t;
  t := IO_gettimeofday();
  return t.tv_sec * 1000000 + t.tv_usec;
end;

# choose which version of Schreier-Sims-Algorithm will be applied
# randomized = true will apply the randomized Schreier-Sims,
# randomited = false will apply the deterministic algorithm
randomized := true;

# group is the short cut for the group that will be investigated
# is used by AtlasGroups to load group information 
group := "M11";

# define your group here:
G := AtlasGroup(group);
#G := SymmetricGroup(100);

# Rubik group 
#G:= Group(
#( 1, 3, 8, 6)( 2, 5, 7, 4)( 9,33,25,17)(10,34,26,18)(11,35,27,19),
#( 9,11,16,14)(10,13,15,12)( 1,17,41,40)( 4,20,44,37)( 6,22,46,35),
#(17,19,24,22)(18,21,23,20)( 6,25,43,16)( 7,28,42,13)( 8,30,41,11),
#(25,27,32,30)(26,29,31,28)( 3,38,43,19)( 5,36,45,21)( 8,33,48,24),
#(33,35,40,38)(34,37,39,36)( 3, 9,46,32)( 2,12,47,29)( 1,14,48,27),
#(41,43,48,46)(42,45,47,44)(14,22,30,38)(15,23,31,39)(16,24,32,40) );

S := Set(GeneratorsOfGroup(G));
n := LargestMovedPoint(G);
B := [];

Print("*********************************************************\n");
if randomized then
    Print("Apply randomized Schreier-Sims-Algorithm to initial sets...\n\n");
else
   Print("Apply deterministic Schreier-Sims-Algorithm to initial sets...\n\n"); 
fi;
Print("Info:\n");
Print(" - Group: "); Print(group); Print(" on "); Print(n); Print(" points\n");
Print(" - Size: "); Print(Size(G)); Print("\n");

t1 := MicroSeconds();
if randomized then
    BSGS := RandomSchreierSims(B, S, 2, n);
else
    BSGS := SchreierSims(B, S, n);
fi;
Add(BSGS[2], [()]);
Print("\nProcedure done. Possible BSGS calculated. Time needed: "); Print((MicroSeconds()-t1) * 1.0 / 1000000); Print("\n");

Print("Base: "); Print(BSGS[1]); Print("\n");
Print("Size of orbits: "); Print(BSGS[3]); Print("\n");

if randomized then
    Print("Length of SGS: "); Print(Length(BSGS[2][1])); Print("\n");
    Print("Remove redundant generators from SGS ...\n");
    RemoveRedundantGens(BSGS, S, n);
    Print("SGS reduced. New length of SGS: "); Print(Length(BSGS[2][1])); Print("\n");
fi;


#dummy := []; 
#for i in [1..Length(BSGS[1])] do
#    Add(dummy, OrbitSV(BSGS[1][i], BSGS[2][i], n));
#od;
#Add(BSGS, dummy);

#MakeImmutable(BSGS);
