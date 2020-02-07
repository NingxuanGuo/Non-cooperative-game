function dsre = myrevert(customer_id,ds)
%% Reduction model

global shiftload

j = customer_id;
k = find(shiftload(:,5)==j);
Ns = size(k,1);

yalmip('clear')

% Define variables
x = binvar(Ns,24,'full');

% Define constraints 
F = sum(x,2)==shiftload(k,2);
%Time window
% for b = 1 : Ns
%     F = [F,sum(x(b,1:shiftload(b,3)-1))==0];
%     F = [F,sum(x(b,shiftload(b,4)+1:24))==0];
% end

% Define an objective
f = (shiftload(k,1)' * x - ds')*(shiftload(k,1)' * x - ds')';


% Set some options for YALMIP and solver
%ops = sdpsettings('debug','1');

% Solve the problem
sol = optimize(F,f);

shift = value(x);%Shiftable loads scheduling timg-slot
dsre = shiftload(k,1)'*shift;%shiftable loads power
end

