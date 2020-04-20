Read("SchreierSims.gi");

#######################################################################
# Description: The deterministic verification, parallelized w.r.t.    #       
#              base points, described as method 1 in section 5.1.     #                             
# Input: - BSGS data structure returned by RandomSchreierSims         # 
# Output: - boolean, true if BSGS is correct, else false              #
#######################################################################

Method1 := function(BSGS)
    local l, tasks, results;
    
    Print("Start verification.\n");
    Print("Apply parallel computing method w.r.t. base points.\n");
    
    tasks := [1..Length(BSGS[1])];  results := [1..Length(BSGS[1])];
    
    for l in [1..Length(BSGS[1])] do
        tasks[l] := RunTask(SimpleVerify, BSGS, l);
    od;
    for l in [1..Length(BSGS[1])] do
        results[l] := TaskResult(tasks[l]);
    od;
    
    Print("\nVerification done.");
    
    if false in results then
        Print("\nBSGS is wrong!\n");
    else
        Print("\nCorrectness of BSGS verified!\n");
    fi;
    Print("*********************************************************\n");
end;

#######################################################################
#######################################################################

#######################################################################
# Description: The deterministic verification, parallelized w.r.t.    #       
#              orbit points, described as method 2 in section 5.1.    #                             
# Input: - BSGS data structure returned by RandomSchreierSims         # 
# Output: - boolean, true if BSGS is correct, else false              #
#######################################################################

Method2 := function(BSGS)
    local l, tasks, results;
    
    Print("Start verification.\n");
    Print("Apply parallel computing method w.r.t. orbit points.\n");
    
    tasks := [1..Length(BSGS[1])];  results := [1..Length(BSGS[1])];
        
    for l in Reversed([1..Length(BSGS[1])]) do
        results[l] := ParallelVerify(BSGS, l);
        if false in results then
            break;
        fi;
    od;
    
    Print("\nVerification done.");
    
    if false in results then
        Print("\nBSGS is wrong!\n");
    else
        Print("\nCorrectness of BSGS verified!\n");
    fi;
    Print("*********************************************************\n");
end;
