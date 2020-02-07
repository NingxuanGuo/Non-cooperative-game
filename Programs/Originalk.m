function Originalk()
%Without DR case for k-parameterization
%The m'th iteration is the results of without DR case

global shiftload;
global fixload;
global satisfaction_1;
global cost_a1;
global cost_a2;
global oriD;

Nc = size(fixload,1);
% Initialization
e = 10^-6;
h = 0.5 * 0.5 * ones(24,1);

D0 = sum(oriD,2)';

% Companies parameters initialization
b = 1./(2*cost_a2);
a = cost_a1;
k = 1*ones(3,24);

% price calculation
r0 = (D0 + sum(b.*a))./ sum(b./k);


%calculate customers' initial bills and payoff
% d0 = oriD';
% for j = 1 : Nc
%     pay0(1,j) = satisfaction_1(j,:) * d0(j,:)' - 0.5 * 0.5 * d0(j,:) * d0(j,:)' - d0(j,:) * r0';
% end

ds0 = zeros(Nc,24);
for j = 1 : Nc
    y = find(shiftload(:,5)==j);
    S = size(y,1);
    for o = 1 : S
        ds0(j,shiftload(y(o),6):shiftload(y(o),6)+shiftload(y(o),2)-1) = ds0(j,shiftload(y(o),6):shiftload(y(o),6)+shiftload(y(o),2)-1)+shiftload(y(o),1);
    end
end

%calculate companies' initial profit
% pro0 = zeros(1,3);
% cost0 = zeros(1,3);
% for i = 1 : 3
%     s = b(i,:)./k(i,:).* (r0 - k(i,:).*a(i,:));
%     cost0(1,i) = cost_a2(i,:) .* s * s' + cost_a1(i,:) * s';
%     pro0(1,i) = s * r0' - cost0(1,i);
% end

% Optimization
price = r0;
m = 1;
r(m,:) = price;
K(:,:,m) = k;
c = 100;

tic
while c > e
    m = m + 1;
    % company
    for i = 1 : 3
         s0 = b(i,:)./K(i,:,m-1).*(price-K(i,:,m-1).*a(i,:));
         gg = price - s0./(sum(b./K(:,:,m-1))-b(i,:)./K(i,:,m-1))-(2*cost_a2(i,:).*s0+cost_a1(i,:));
         gk = -price.*b(i,:)./(K(i,:,m-1).^2).*(sum(b./K(:,:,m-1))-b(i,:)./K(i,:,m-1))./sum(b./K(:,:,m-1)).*gg;
         K(i,:,m) = K(i,:,m-1) + 0.05 * gk;
    end
    r(m,:) = (D0+ sum(b.*a))./ sum(b./K(:,:,m));
    c2 = max(max(abs(K(:,:,m)-K(:,:,m-1))));
    c4 = max(abs(r(m,:)-r(m-1,:)));
    c = max([c2,c4]);
    price = r(m,:);
    s = b./K(:,:,m) .* (ones(3,1) * price - a.*K(:,:,m));
    cost0 = sum(cost_a2 .* s.^2 + cost_a1 .* s,2)';
    pro0 = price * s' - cost0; 
end
toc
    SW0 = sum(diag((fixload + ds0)* satisfaction_1') - (fixload + ds0).^2 * h)-sum(sum(K(:,:,m)./(2*b).*s.^2+K(:,:,m).*a.*s));
    a0=a;
    b0=b;
    k0=K(:,:,m);
end

