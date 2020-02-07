function f = mydisDatacenter11()
%(beta,alpha)-parameterization


global shiftload;
global fixload;
global satisfaction_1;
global cost_a1;
global cost_a2;
global oriD;
%load ga_x.mat
load originalba.mat

Nf =3;
Nc = 10;
% Initialization
e = 10^-6;
h = 0.5 * 0.5 * ones(24,1);

%D0 = sum(oriD,2)'*Nc/10;
%D0 = oriD(:,5)'*Nc;

% Companies parameters initialization
%  b0 = reshape(ga_x(1:72),[3,24]);
%  a0 = reshape(ga_x(73:144),[3,24]);
%  cost_a2 = 0.09+0.02*rand(Nf,24);
%  cost_a1 = 0.1+0.2*rand(Nf,24);
%  b0 = 1./(2*cost_a2);
%  a0 = cost_a1;

% price calculation
%r0 = (D0 + diag(b0' * a0)')./ sum(b0);
r0 = price;


%calculate customers' initial bills and payoff
d0 = oriD';
%bill0 = d0 * r0';
for j = 1 : Nc
    pay0(1,j) = satisfaction_1(j,:) * d0(j,:)' - 0.5 * 0.5 * d0(j,:) * d0(j,:)' - d0(j,:) * r0';
end

% ds0 = zeros(Nc,24);
% for j = 1 : Nc
% %     k = find(shiftload(:,5)==j);
%     o = mod(j,10);
%     if o==0
%         o=10;
%     end
%     k = find(shiftload(:,5)==o);
%     S = size(k,1);
%     for a = 1 : S
%         ds0(j,shiftload(k(a),6):shiftload(k(a),6)+shiftload(k(a),2)-1) = ds0(j,shiftload(k(a),6):shiftload(k(a),6)+shiftload(k(a),2)-1)+shiftload(k(a),1);
%     end
% end

%calculate companies' initial profit
% pro0 = zeros(1,3);
% cost0 = zeros(1,3);
% for i = 1 : 3
%     s = b0(i,:) .* (r0 - a0(i,:));
%     cost0(1,i) = cost_a2(i,:) .* s * s' + cost_a1(i,:) * s';
%     pro0(1,i) = s * r0' - cost0(1,i);
% end

% Optimization
price = r0;
m = 1;
r(m,:) = price;
B(:,:,m) = b0;
A(:,:,m) = a0;
D(m,:) = D0;
pro(m,:) = pro0;
payoff(m,:) = pay0;
cost(m,:) = cost0;
x(:,:,1) = ds0';
c = 100;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   

tic
while c > e
    m = m + 1;
    % company
    for i = 1 : Nf 
%s=bt(p-at),beta,alpha-parameterization
         s0 = B(i,:,m-1).*(price-A(i,:,m-1));
         gg = price - s0./(sum(B(:,:,m-1))-B(i,:,m-1))-(2*cost_a2(i,:).*s0+cost_a1(i,:));
         gb = (price - A(i,:,m-1)).*(sum(B(:,:,m-1))-B(i,:,m-1))./sum(B(:,:,m-1)).*gg;
         ga = (B(i,:,m-1) - sum(B(:,:,m-1))).* B(i,:,m-1)./sum(B(:,:,m-1)).*gg;
         B(i,:,m) = B(i,:,m-1) + 0.1 * gb;
         A(i,:,m) = A(i,:,m-1) + 0.08 * ga;
%          s = B(i,:,m) .* (price - A(i,:,m));
%          cost(m,i) = cost_a2(i,:) .* s * s' + cost_a1(i,:) * s';%cost_1=[0.09,0.11]
%         pro(m,i) = s * price' - cost(m,i); 
    end
%customer
    D(m,:) = sum(fixload)*Nc/10;
    for j = 1 : Nc
        x(:,j,m) = myProjection2(j,price,x(:,j,m-1),B(:,:,m-1));
        D(m,:) = D(m,:) + x(:,j,m)';
%         [d,pay] = customer2(j,price,x(:,j,m-1),B(:,:,m-1),A(:,:,m-1),D(m-1,:));
%         x(:,j,m) = d-fixload(j,:)';
%         D(m,:) = D(m,:) +  x(:,j,m)';
    end
    r(m,:) = (D(m,:) + diag(B(:,:,m)' * A(:,:,m))')./ sum(B(:,:,m));
     c1 = max(max(abs(x(:,:,m)-x(:,:,m-1))));
     c2 = max(max(abs(B(:,:,m)-B(:,:,m-1))));
     c3 = max(abs(r(m,:)-r(m-1,:)));
     c = max([c1,c2,c3]);   
    price = r(m,:);
    s = B(:,:,m) .* (ones(3,1) * price - A(:,:,m));
    cost(m,:) = sum(cost_a2 .* s.^2 + cost_a1 .* s,2)';
    pro(m,:) = price * s' - cost(m,:); 
    payoff(m,:) = (diag((fixload + x(:,:,m)') * satisfaction_1') - (fixload + x(:,:,m)').^2 * h - (fixload + x(:,:,m)') * price')';
end
toc   
    %f = -SW;
    SW1 = sum(diag((fixload + x(:,:,m)') * satisfaction_1') - (fixload + x(:,:,m)').^2 * h)-sum(sum(1./(2*B(:,:,m)).*s.^2+A(:,:,m).*s));
    equivalent_ds = reshape(x(:,:,m),[24,Nc]);
    dsre = zeros(Nc,24);
    Dre = sum(fixload);
    for j = 1:Nc
        dsre(j,:) = myrevert(j, equivalent_ds(:,j));
        Dre = Dre + dsre(j,:);
    end
    price = (Dre + diag(B(:,:,m)' * A(:,:,m))')./ sum(B(:,:,m));
    s = B(:,:,m) .* (ones(3,1) * price - A(:,:,m));
    costfval = sum(cost_a2 .* s.^2 + cost_a1 .* s,2)';
    profval = price * s' - costfval; 
    payofffval = (diag((fixload + x(:,:,m)') * satisfaction_1') - (fixload + x(:,:,m)').^2 * h - (fixload + x(:,:,m)') * price')';
    SW = sum(diag((fixload + x(:,:,m)') * satisfaction_1') - (fixload + x(:,:,m)').^2 * h)-sum(sum(1./(2*B(:,:,m)).*s.^2+A(:,:,m).*s));
end

