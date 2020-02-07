clear;

%Load all global variables

global shiftload; %shiftable load data
%The content of one row in shiftload is: power,duration,minimum time
%limit,maximum time limit,customer_id,the original start supply time
global satisfaction_1; %linear coefficient for satisfaction functions of customers,randomly selected from [4,10]
global fixload;% fixed load data
%Each row in fixload is each customer's fixed load data at each time-slot
global oriD;% The original demand load curve
global cost_a1;%linear coefficient for cost functions of firms, randomly selected from [0.1,0.3]
global cost_a2;%2*quadratic coefficient for cost functions of firms,randomly selected from [0.09,0.11]

%shiftload = cell2mat(struct2cell(load('NREL_shift.mat')));
satisfaction_1 = cell2mat(struct2cell(load('satisfaction_1.mat')));
%fixload = cell2mat(struct2cell(load('new_fix.mat')));
shiftload = cell2mat(struct2cell(load('shift0101.mat')));
fixload = cell2mat(struct2cell(load('fffix.mat')));
%cost_1 = cell2mat(struct2cell(load('cost_1.mat')));
cost_a1 = cell2mat(struct2cell(load('cost_a1.mat')));
cost_a2 = cell2mat(struct2cell(load('cost_a2.mat')));
oriD = cell2mat(struct2cell(load('loadfile.mat')));

% Constraint matrix parameters for the gradient projection method
global E;
global b;
global Q;
global A;
load('Projection_parameters_no_timewindow.mat','E','b','Q','A'); 


