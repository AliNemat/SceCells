function [frhs] = Dpp_signaling(Dpp_mat,A_mat,L_mat,NI_mat,Dpp_source)
% Dpp_vec = [Dpp, Tkv, Dpp-Tkv, pMad]

bio_parameters;

Dpp = Dpp_mat(:,1);
Tkv = Dpp_mat(:,2);
DT  = Dpp_mat(:,3);
pMad = Dpp_mat(:,4);

mesh_size = length(Dpp);

Diff_Dpp = D_dpp*sum((repmat(Dpp',mesh_size,1)-repmat(Dpp,1,mesh_size)).*A_mat./L_mat.*NI_mat,2);

Dpp_rhs = Diff_Dpp - kon*Dpp.*Tkv + koff*DT + Prod_dpp*Dpp_source - d_dpp*Dpp;
Tkv_rhs = - kon*Dpp.*Tkv + koff*DT + Prod_tkv./(1+(pMad/kp).^n) - d_tkv*Tkv;
DT_rhs  = kon*Dpp.*Tkv - koff*DT - d_dt*DT;
pMad_rhs = Prod_pmad./(1+(DT/kdt).^(-n)) - d_pmad*pMad;

frhs = [Dpp_rhs, Tkv_rhs, DT_rhs, pMad_rhs];