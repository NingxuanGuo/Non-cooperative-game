function myGA()
%Select the solution with maximum social welfare for
%(beta,alpha)-parameterization
%The calculated result x is the start point for (beta,alpha)-parameterization

ObjectiveFunction = @mydisDatacenter11;
ObjectiveFunction = @Original;
nvars = 312;%Number of variables
LB = [0.6*ones(72,1);zeros(72,1)];%lower bounds
UB = [5 * ones(72,1);3 * ones(72,1)];%upper bounds
options = gaoptimset('PopulationSize',20,'Generations',100); 
[x,fval] = ga(ObjectiveFunction,nvars,[],[],[],[],LB,UB,[],options)
end

