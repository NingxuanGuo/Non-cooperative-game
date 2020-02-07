function f = Original()
%Without DR case for beta-parameterization, alpha-parameterization, 
%and (beta,alpha)-parameterization
%The m'th iteration is the results of without DR case

global shiftload;
global fixload;
global satisfaction_1;
global cost_a1;
global cost_a2;
global oriD;

Nc = size(fixload,1);
% Initialization
e = 10^-9;
h = 0.5 * 0.5 * ones(24,1);

D0 = sum(oriD,2)';

% Companies parameters initialization
b = 1./(2*cost_a2);
a = cost_a1;
% b = reshape(ga_x(1:72),[3,24]);
% a = reshape(ga_x(73:144),[3,24]);

% price calculation
r = (D0 + diag(b' * a)')./ sum(b);

ds0 = zeros(Nc,24);
for j = 1 : Nc
    k = find(shiftload(:,5)==j);
    S = size(k,1);
    for o = 1 : S
        ds0(j,shiftload(k(o),6):shiftload(k(o),6)+shiftload(k(o),2)-1) = ds0(j,shiftload(k(o),6):shiftload(k(o),6)+shiftload(k(o),2)-1)+shiftload(k(o),1);
    end
end

% Optimization
price = r;
m = 1;
r(m,:) = price;
c = 100;  
%% b-para
B(:,:,m) = b;

tic
while c > e
    m = m + 1;
    % company
    for i = 1 : 3
         s0 = B(i,:,m-1).*(price-a(i,:));
         gg = price - s0./(sum(B(:,:,m-1))-B(i,:,m-1))-(2*cost_a2(i,:).*s0+cost_a1(i,:));
         gb = (price - a(i,:)).*(sum(B(:,:,m-1))-B(i,:,m-1))./sum(B(:,:,m-1)).*gg;
         B(i,:,m) = B(i,:,m-1) + 0.1 * gb;
    end
    r(m,:) = (D0 + diag(B(:,:,m)' * a)')./ sum(B(:,:,m));
    c2 = max(max(abs(B(:,:,m)-B(:,:,m-1))));
    c4 = max(abs(r(m,:)-r(m-1,:)));
    c = max([c2,c4]); 
    price = r(m,:);
    s = B(:,:,m) .* (ones(3,1) * price - a);
    cost0 = sum(cost_a2 .* s.^2 + cost_a1 .* s,2)';
    pro0 = price * s' - cost0; 
end
toc
    SW0 = sum(diag((fixload + ds0) * satisfaction_1') - (fixload + ds0).^2 * h)-sum(sum(1./(2*B(:,:,m)).*s.^2+a.*s));
    b0 = B(:,:,m);
    a0 = a;
    
 %% a-para
% A(:,:,m) = a;
% 
% tic
% while c > e
%     m = m + 1;
%     % company
%     for i = 1 : 3
%          s0 = b(i,:).*(price-A(i,:,m-1));
%          gg = price - s0./(sum(b)-b(i,:))-(2*cost_a2(i,:).*s0+cost_a1(i,:));
%          ga = (b(i,:) - sum(b)).* b(i,:)./sum(b).*gg;
%          A(i,:,m) = A(i,:,m-1) + 0.1 * ga;
%     end
%     r(m,:) = (D0 + diag(b' * A(:,:,m))')./ sum(b);
%     c2 = max(max(abs(A(:,:,m)-A(:,:,m-1))));
%     c4 = max(abs(r(m,:)-r(m-1,:)));
%     c = max([c2,c4]); 
%     price = r(m,:);
%     s = b .* (ones(3,1) * price - A(:,:,m));
%     cost0 = sum(cost_a2 .* s.^2 + cost_a1 .* s,2)';
%     pro0 = price * s' - cost0; 
% end
% toc
%     SW0 = sum(diag((fixload + ds0) * satisfaction_1') - (fixload + ds0).^2 * h)-sum(sum(1./(2*b).*s.^2+A(:,:,m).*s));
%     b0 = b;
%     a0 = A(:,:,m);
    
%% ba-para
% B(:,:,m) = b;
% A(:,:,m) = a;
% 
% tic
% while c > e
%     m = m + 1;
%     % company
%     for i = 1 : 3
%          s0 = B(i,:,m-1).*(price-A(i,:,m-1));
%          gg = price - s0./(sum(B(:,:,m-1))-B(i,:,m-1))-(2*cost_a2(i,:).*s0+cost_a1(i,:));
%          gb = (price - A(i,:,m-1)).*(sum(B(:,:,m-1))-B(i,:,m-1))./sum(B(:,:,m-1)).*gg;
%          ga = (B(i,:,m-1) - sum(B(:,:,m-1))).* B(i,:,m-1)./sum(B(:,:,m-1)).*gg;
%          B(i,:,m) = B(i,:,m-1) + 0.1 * gb;
%          A(i,:,m) = A(i,:,m-1) + 0.05 * ga;
%     end
%     r(m,:) = (D0 + diag(B(:,:,m)' * A(:,:,m))')./ sum(B(:,:,m));
%     c2 = max(max(abs(B(:,:,m)-B(:,:,m-1))));
%     c3 = max(max(abs(A(:,:,m)-A(:,:,m-1))));
%     c4 = max(abs(r(m,:)-r(m-1,:)));
%     c = max([c2,c3,c4]); 
%     price = r(m,:);
%     s = B(:,:,m) .* (ones(3,1) * price - A(:,:,m));
%     cost0 = sum(cost_a2 .* s.^2 + cost_a1 .* s,2)';
%     pro0 = price * s' - cost0; 
% end
% toc
%     SW0 = sum(diag((fixload + ds0) * satisfaction_1') - (fixload + ds0).^2 * h)-sum(sum(1./(2*B(:,:,m)).*s.^2+A(:,:,m).*s));
%     b0 = B(:,:,m);
%     a0 = A(:,:,m);
%     f = -SW0;
end




