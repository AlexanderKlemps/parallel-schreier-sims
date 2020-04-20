Read("SchreierSims.gi");
LoadPackage( "AtlasRep" );

MicroSeconds := function()
  local t;
  t := IO_gettimeofday();
  return t.tv_sec * 1000000 + t.tv_usec;
end;

n:=11;

Groups := ["M11", "M23", "M24", "J1", "J2", "J3", "Suz", "Co2", "Fi22"];
#Groups := [120];

for group in Groups do
    counter := 0;
    time := [[],[],[],[]];  
    Print("Group: ");  Print(group); 
    while counter < 5 do
        G := AtlasGroup(group);
        #G := SymmetricGroup(group);
        S := Set(GeneratorsOfGroup(G)); n := LargestMovedPoint(G); B := [];
        t1 := MicroSeconds();
        BSGS := RandomSchreierSims(B, S, 10, n);
        t2 := (MicroSeconds() - t1) * 1.0 / 1000000;
        RemoveRedundantGens(BSGS, S, n);
        upper_bound := Length(BSGS[1]);
        results := [1..Length(BSGS[1])];
        
        ## parallelization w.r.t. base
        t1 := MicroSeconds();
        tasks := [1..Length(BSGS[1])];
        for l in [1..upper_bound] do
            tasks[l] := RunTask(SimpleVerify, BSGS, l);
        od;
        for l in [1..upper_bound] do
            results[l] := TaskResult(tasks[l]);
        od;
        if false in results then
            continue; 
        fi;
        Add(time[4], t2);
        counter := counter + 1;
        Add(time[1], (MicroSeconds() - t1) * 1.0 / 1000000);
        
        results := [1..Length(BSGS[1])];
        
        ## parallelization w.r.t. base
        t1 := MicroSeconds();
        for l in Reversed([1..upper_bound]) do
            results[l] := ParallelVerify(BSGS, l);            
        od;
        Add(time[2], (MicroSeconds() - t1) * 1.0 / 1000000);
        
        results := [1..Length(BSGS[1])];

        ## no parallelization
        t1 := MicroSeconds();
        for l in Reversed([1..upper_bound]) do
            results[l] := SimpleVerify(BSGS, l);     
        od;
        Add(time[3], (MicroSeconds() - t1) * 1.0 / 1000000);
    od;
    average := [0,0,0,0];
    for i in [1..5] do
        average[1] := average[1] + time[1][i]/5;
        average[2] := average[2] + time[2][i]/5;
        average[3] := average[3] + time[3][i]/5;
        average[4] := average[4] + time[4][i]/5;
    od;
    Print(", calculation time: "); Print(average[4]); Print(", times: "); Print(average{[1..3]}); Print("\n");
od;
