function Calculate_parameters()
%Calculate parameters for the gradient projection method
%Calculate shiftable load related constraint matrix
%The calculated results is saved in Projection_parameters_no_timewindow.mat
global shiftload
global fixload

Nc = size(fixload,1);
b = zeros(48,Nc);
Q = zeros(1,Nc);
for j = 1 : Nc
    k = find(shiftload(:,5)==j);
    Ns = size(k,1);
    
%     tmin = min(shiftload(k,3));
%     tmax = max(shiftload(k,4));
%     Nt = tmax - tmin + 1;
%     t = (tmin : tmax)';
%     dsmin = zeros(Nt,1);
%     dsmax = zeros(Nt,1);

    t = (1 : 24)';
    dsmin = zeros(24,1);
    dsmax = zeros(24,1);
    
%     for a = 1 : Ns
%         dsmax = dsmax + sign(max((t-shiftload(k(a),3)+1).*(shiftload(k(a),4)-t+1),0)) .* shiftload(k(a),1);
%     end

    for a = 1 : Ns
        dsmax = dsmax + sign(max((t-1+1).*(24-t+1),0)) .* shiftload(k(a),1);
    end
    Q(1,j) = shiftload(k,1)'*shiftload(k,2);
    b(:,j) = [dsmin;-dsmax];
    
end
    A = [eye(24,24);-1 * eye(24,24)];
    E = ones(1,24);
end

