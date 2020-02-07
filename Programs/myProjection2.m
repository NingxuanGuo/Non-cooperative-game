function f = myProjection2(customer_id,price,x,B)
% Gradient projection method for each customer
% x is the shiftable load demand

global shiftload
global fixload
global satisfaction_1
global E
global Q
global b
global A
% A = eye(24,24);
% b = zeros(24,10);

%j = customer_id;
j = mod(customer_id,10);
if j==0
    j=10;
end
e = 10^-3;
%step 2
o = find(abs(A * x - b(:,j))<=e);
A1 = A(o,:);
A2 = A;
A2(o,:) = [];
b1 = b(o,j);
b2 = b(:,j);
b2(o,:) = [];
bb = b2 - A2 * x;

G = diag(0.5+2./sum(B));
c = (-satisfaction_1(j,:)-(fixload(j,:)+x')./sum(B)+price)';
% G = -diag(0.5+2./sum(B));
% c = (satisfaction_1(j,:) - (D-fixload(j,:)-x')./sum(B) - diag(B'*a)'./sum(B))';
g = G * (x + fixload(j,:)') + c;
%step3
while 1
L = size(A1,1);
M = [A1;E];
if rank(M) == rank(A1)
    M = A1;
end
%M = A1;
P = eye(24,24) - M'* (M*M')^(-1) *M;
d = - P * g;
if norm(d) <= e
    W = (M*M')^(-1) * M * g;
    if W(1:L,:) >= 0
    %if W(1:L,:) <= 0
        f = x;
        break;
    else
        if ~isempty(A1)
            %k = find(W(1:L,:)<0);
            k = find(W(1:L,:)==min(W(1:L,:)));
            Nk = size(k,1);
            A1(k(Nk),:) = [];
        else
            f = x;
            break;
        end        
    end
else
    dd = A2 * d;
    if dd >= 0
        rmax = inf;     
    else
        i = find(dd<0);
        rmax = min(bb(i,:)./dd(i,:));
    end
    rr = 1;
    %rr = - (0.5 *x'* G * d+ 0.5*d'*G*x+fixload(j,:)* G*d+c'*d)/(d'*G*d);
    if rr > rmax
        r = rmax;
    elseif rr < 0
        r = 0;
    else
        %r = 0.6;
         r = rr;     
    end
    for h =0:100
        l=0.5^(h+1)*rr*g'*d;
        m=0.5*(x + 0.5^h*rr*d+fixload(j,:)')'*G*(x + 0.5^h*rr*d+fixload(j,:)')+c'*(x + 0.5^h*rr*d+fixload(j,:)');
        n=0.5*(x+fixload(j,:)')'*G*(x+fixload(j,:)')+c'*(x+fixload(j,:)');
        y = l-m-n;
        if y>=0
        	break;
        end
    end
     f = x + 0.5^h * r * d;
    %f = x + r * d;
    break;
end
end
 
end


