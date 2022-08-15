function z_gate= gate_meas_gms(z,gamma,model,m,P)

valid_idx = [];
zlength = size(z,2); if zlength==0, z_gate= []; return; end
plength = size(m,2);

for j=1:plength

    [H_ekf, ~]= ekf_update_mat(model, m(:, j));                 % user specified function for application

    Sj= model.R + H_ekf * P(:,:,j) * H_ekf';
    Vs= chol(Sj);
    %det_Sj= prod(diag(Vs))^2;
    inv_sqrt_Sj= inv(Vs);
    %iSj= inv_sqrt_Sj*inv_sqrt_Sj';
    nu= z - H_ekf * repmat(m(:,j), [1 zlength]);

    % _/!\_ We must check that the difference between two angles is between -pi and pi _/!\_
    is_gg_pi = nu(1,:) > pi;
    nu(1, is_gg_pi) = nu(1, is_gg_pi) - 2*pi;
    is_ll_pi = nu(1,:) < -pi;
    nu(1, is_ll_pi) = nu(1, is_ll_pi) + 2*pi;

    dist= sum((inv_sqrt_Sj'*nu).^2);
    valid_idx= union(valid_idx,find( dist < gamma ));
end
z_gate = z(:,valid_idx);
