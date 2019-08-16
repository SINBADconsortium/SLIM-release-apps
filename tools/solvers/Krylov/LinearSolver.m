classdef LinearSolver < FuncObj   
    
    methods
        function op = LinearSolver(lin_system,solve_opts)
            op = op@FuncObj(@linearsolve,{lin_system,[],[],[],solve_opts},{'A','b','x','mode','solve_opts'});
        end                
    end
end