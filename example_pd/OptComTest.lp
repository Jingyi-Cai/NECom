\ Model matlab
\ LP format - for model browsing. Use MPS format to capture full model detail.
Maximize
  v_1_Biomass_sp1 + v_1_Biomass_sp2
Subject To
 SV=0_S[e]_sp1_lb: - v_1_EX_S_sp1 - v_1_T_S_sp1 - v_1_EX_S_2_sp1 >= 0
 SV=0_A[e]_sp1_lb: - v_1_EX_A_sp1 - v_1_T_A_sp1 - v_1_EX_A_2_sp1 >= 0
 SV=0_B[e]_sp1_lb: - v_1_EX_B_sp1 - v_1_T_B_sp1 - v_1_EX_B_2_sp1 >= 0
 SV=0_P[e]_sp1_lb: - v_1_EX_P_sp1 - v_1_T_P_sp1 - v_1_EX_P_2_sp1 >= 0
 SV=0_h[e]_sp1_lb: - v_1_EX_h_sp1 - v_1_T_h_sp1 - v_1_EX_h_2_sp1 >= 0
 SV=0_pi[e]_sp1_lb: - v_1_EX_pi_sp1 - v_1_T_pi_sp1 - v_1_EX_pi_2_sp1 >= 0
 SV=0_h2o[e]_sp1_lb: - v_1_EX_h2o_sp1 - v_1_T_h2o_sp1 - v_1_EX_h2o_2_sp1
   >= 0
 SV=0_S[c]_sp1_lb: v_1_T_S_sp1 - v_1_R_1_sp1 >= 0
 SV=0_A[c]_sp1_lb: v_1_T_A_sp1 + v_1_R_3_sp1 - 3 v_1_OBJ_sp1 >= 0
 SV=0_B[c]_sp1_lb: v_1_T_B_sp1 + v_1_R_2_sp1 - 3 v_1_OBJ_sp1 >= 0
 SV=0_h[c]_sp1_lb: v_1_T_h_sp1 - 2 v_1_R_1_sp1 + v_1_R_2_sp1
   + 3 v_1_R_3_sp1 + v_1_R_4_sp1 + 5 v_1_OBJ_sp1 >= 0
 SV=0_P[c]_sp1_lb: v_1_T_P_sp1 + v_1_R_1_sp1 - v_1_R_2_sp1 - v_1_R_3_sp1
   >= 0
 SV=0_pi[c]_sp1_lb: v_1_T_pi_sp1 + v_1_R_4_sp1 >= 0
 SV=0_h2o[c]_sp1_lb: v_1_T_h2o_sp1 - v_1_R_4_sp1 >= 0
 SV=0_adp[c]_sp1_lb: - 2 v_1_R_1_sp1 + v_1_R_2_sp1 + 3 v_1_R_3_sp1
   + v_1_R_4_sp1 + 5 v_1_OBJ_sp1 >= 0
 SV=0_atp[c]_sp1_lb: 2 v_1_R_1_sp1 - v_1_R_2_sp1 - 3 v_1_R_3_sp1
   - v_1_R_4_sp1 - 5 v_1_OBJ_sp1 >= 0
 SV=0_biomass[c]_sp1_lb: v_1_OBJ_sp1 - v_1_Biomass_sp1 >= 0
 SV=0_S[e]_sp2_lb: - v_1_EX_S_sp2 - v_1_T_S_sp2 - v_1_EX_S_2_sp2 >= 0
 SV=0_A[e]_sp2_lb: - v_1_EX_A_sp2 - v_1_T_A_sp2 - v_1_EX_A_2_sp2 >= 0
 SV=0_B[e]_sp2_lb: - v_1_EX_B_sp2 - v_1_T_B_sp2 - v_1_EX_B_2_sp2 >= 0
 SV=0_P[e]_sp2_lb: - v_1_EX_P_sp2 - v_1_T_P_sp2 - v_1_EX_P_2_sp2 >= 0
 SV=0_h[e]_sp2_lb: - v_1_EX_h_sp2 - v_1_T_h_sp2 - v_1_EX_h_2_sp2 >= 0
 SV=0_pi[e]_sp2_lb: - v_1_EX_pi_sp2 - v_1_T_pi_sp2 - v_1_EX_pi_2_sp2 >= 0
 SV=0_h2o[e]_sp2_lb: - v_1_EX_h2o_sp2 - v_1_T_h2o_sp2 - v_1_EX_h2o_2_sp2
   >= 0
 SV=0_S[c]_sp2_lb: v_1_T_S_sp2 - v_1_R_1_sp2 >= 0
 SV=0_A[c]_sp2_lb: v_1_T_A_sp2 + v_1_R_3_sp2 - 3 v_1_OBJ_sp2 >= 0
 SV=0_B[c]_sp2_lb: v_1_T_B_sp2 + v_1_R_2_sp2 - 3 v_1_OBJ_sp2 >= 0
 SV=0_h[c]_sp2_lb: v_1_T_h_sp2 - 2 v_1_R_1_sp2 + 3 v_1_R_2_sp2
   + v_1_R_3_sp2 + v_1_R_4_sp2 + 5 v_1_OBJ_sp2 >= 0
 SV=0_P[c]_sp2_lb: v_1_T_P_sp2 + v_1_R_1_sp2 - v_1_R_2_sp2 - v_1_R_3_sp2
   >= 0
 SV=0_pi[c]_sp2_lb: v_1_T_pi_sp2 + v_1_R_4_sp2 >= 0
 SV=0_h2o[c]_sp2_lb: v_1_T_h2o_sp2 - v_1_R_4_sp2 >= 0
 SV=0_adp[c]_sp2_lb: - 2 v_1_R_1_sp2 + 3 v_1_R_2_sp2 + v_1_R_3_sp2
   + v_1_R_4_sp2 + 5 v_1_OBJ_sp2 >= 0
 SV=0_atp[c]_sp2_lb: 2 v_1_R_1_sp2 - 3 v_1_R_2_sp2 - v_1_R_3_sp2
   - v_1_R_4_sp2 - 5 v_1_OBJ_sp2 >= 0
 SV=0_biomass[c]_sp2_lb: v_1_OBJ_sp2 - v_1_Biomass_sp2 >= 0
 V-V_ut+V_ex=0_EX_A_2_sp1_lb: v_1_EX_A_2_sp1 + vut_EX_A_2_sp1
   - vex_EX_A_2_sp1 >= 0
 V-V_ut+V_ex=0_EX_B_2_sp1_lb: v_1_EX_B_2_sp1 + vut_EX_B_2_sp1
   - vex_EX_B_2_sp1 >= 0
 V-V_ut+V_ex=0_EX_A_2_sp2_lb: v_1_EX_A_2_sp2 + vut_EX_A_2_sp2
   - vex_EX_A_2_sp2 >= 0
 V-V_ut+V_ex=0_EX_B_2_sp2_lb: v_1_EX_B_2_sp2 + vut_EX_B_2_sp2
   - vex_EX_B_2_sp2 >= 0
 Xn*W+psi=0_EX_A_2_sp1_lb: psi_EX_A_2_sp1 + 0.5 omega_EX_A_2_sp1 >= 0
 Xn*W+psi=0_EX_B_2_sp1_lb: psi_EX_B_2_sp1 + 0.5 omega_EX_B_2_sp1 >= 0
 Xn*W+psi=0_EX_A_2_sp2_lb: psi_EX_A_2_sp2 + 0.5 omega_EX_A_2_sp2 >= 0
 Xn*W+psi=0_EX_B_2_sp2_lb: psi_EX_B_2_sp2 + 0.5 omega_EX_B_2_sp2 >= 0
 S*lambda_1+mu_UB-mu_LB+psi=cEX_S_sp1_lb: - mu_LB_EX_S_sp1 + mu_UB_EX_S_sp1
   - lambda_S[e]_sp1_1 >= 0
 S*lambda_1+mu_UB-mu_LB+psi=cEX_A_sp1_lb: - mu_LB_EX_A_sp1 + mu_UB_EX_A_sp1
   - lambda_A[e]_sp1_1 >= 0
 S*lambda_1+mu_UB-mu_LB+psi=cEX_B_sp1_lb: - mu_LB_EX_B_sp1 + mu_UB_EX_B_sp1
   - lambda_B[e]_sp1_1 >= 0
 S*lambda_1+mu_UB-mu_LB+psi=cEX_P_sp1_lb: - mu_LB_EX_P_sp1 + mu_UB_EX_P_sp1
   - lambda_P[e]_sp1_1 >= 0
 S*lambda_1+mu_UB-mu_LB+psi=cEX_h_sp1_lb: - mu_LB_EX_h_sp1 + mu_UB_EX_h_sp1
   - lambda_h[e]_sp1_1 >= 0
 S*lambda_1+mu_UB-mu_LB+psi=cEX_pi_sp1_lb: - mu_LB_EX_pi_sp1
   + mu_UB_EX_pi_sp1 - lambda_pi[e]_sp1_1 >= 0
 S*lambda_1+mu_UB-mu_LB+psi=cEX_h2o_sp1_lb: - mu_LB_EX_h2o_sp1
   + mu_UB_EX_h2o_sp1 - lambda_h2o[e]_sp1_1 >= 0
 S*lambda_1+mu_UB-mu_LB+psi=cT_S_sp1_lb: - mu_LB_T_S_sp1 + mu_UB_T_S_sp1
   - lambda_S[e]_sp1_1 + lambda_S[c]_sp1_1 >= 0
 S*lambda_1+mu_UB-mu_LB+psi=cT_A_sp1_lb: - mu_LB_T_A_sp1 + mu_UB_T_A_sp1
   - lambda_A[e]_sp1_1 + lambda_A[c]_sp1_1 >= 0
 S*lambda_1+mu_UB-mu_LB+psi=cT_B_sp1_lb: - mu_LB_T_B_sp1 + mu_UB_T_B_sp1
   - lambda_B[e]_sp1_1 + lambda_B[c]_sp1_1 >= 0
 S*lambda_1+mu_UB-mu_LB+psi=cT_h_sp1_lb: - mu_LB_T_h_sp1 + mu_UB_T_h_sp1
   - lambda_h[e]_sp1_1 + lambda_h[c]_sp1_1 >= 0
 S*lambda_1+mu_UB-mu_LB+psi=cT_P_sp1_lb: - mu_LB_T_P_sp1 + mu_UB_T_P_sp1
   - lambda_P[e]_sp1_1 + lambda_P[c]_sp1_1 >= 0
 S*lambda_1+mu_UB-mu_LB+psi=cT_pi_sp1_lb: - mu_LB_T_pi_sp1 + mu_UB_T_pi_sp1
   - lambda_pi[e]_sp1_1 + lambda_pi[c]_sp1_1 >= 0
 S*lambda_1+mu_UB-mu_LB+psi=cT_h2o_sp1_lb: - mu_LB_T_h2o_sp1
   + mu_UB_T_h2o_sp1 - lambda_h2o[e]_sp1_1 + lambda_h2o[c]_sp1_1 >= 0
 S*lambda_1+mu_UB-mu_LB+psi=cR_1_sp1_lb: - mu_LB_R_1_sp1 + mu_UB_R_1_sp1
   - lambda_S[c]_sp1_1 - 2 lambda_h[c]_sp1_1 + lambda_P[c]_sp1_1
   - 2 lambda_adp[c]_sp1_1 + 2 lambda_atp[c]_sp1_1 >= 0
 S*lambda_1+mu_UB-mu_LB+psi=cR_2_sp1_lb: - mu_LB_R_2_sp1 + mu_UB_R_2_sp1
   + lambda_B[c]_sp1_1 + lambda_h[c]_sp1_1 - lambda_P[c]_sp1_1
   + lambda_adp[c]_sp1_1 - lambda_atp[c]_sp1_1 >= 0
 S*lambda_1+mu_UB-mu_LB+psi=cR_3_sp1_lb: - mu_LB_R_3_sp1 + mu_UB_R_3_sp1
   + lambda_A[c]_sp1_1 + 3 lambda_h[c]_sp1_1 - lambda_P[c]_sp1_1
   + 3 lambda_adp[c]_sp1_1 - 3 lambda_atp[c]_sp1_1 >= 0
 S*lambda_1+mu_UB-mu_LB+psi=cR_4_sp1_lb: - mu_LB_R_4_sp1 + mu_UB_R_4_sp1
   + lambda_h[c]_sp1_1 + lambda_pi[c]_sp1_1 - lambda_h2o[c]_sp1_1
   + lambda_adp[c]_sp1_1 - lambda_atp[c]_sp1_1 >= 0
 S*lambda_1+mu_UB-mu_LB+psi=cOBJ_sp1_lb: - mu_LB_OBJ_sp1 + mu_UB_OBJ_sp1
   - 3 lambda_A[c]_sp1_1 - 3 lambda_B[c]_sp1_1 + 5 lambda_h[c]_sp1_1
   + 5 lambda_adp[c]_sp1_1 - 5 lambda_atp[c]_sp1_1
   + lambda_biomass[c]_sp1_1 >= 0
 S*lambda_1+mu_UB-mu_LB+psi=cBiomass_sp1_lb: - mu_LB_Biomass_sp1
   + mu_UB_Biomass_sp1 - lambda_biomass[c]_sp1_1 >= 1
 S*lambda_1+mu_UB-mu_LB+psi=cEX_S_sp2_lb: - mu_LB_EX_S_sp2 + mu_UB_EX_S_sp2
   - lambda_S[e]_sp2_1 >= 0
 S*lambda_1+mu_UB-mu_LB+psi=cEX_A_sp2_lb: - mu_LB_EX_A_sp2 + mu_UB_EX_A_sp2
   - lambda_A[e]_sp2_1 >= 0
 S*lambda_1+mu_UB-mu_LB+psi=cEX_B_sp2_lb: - mu_LB_EX_B_sp2 + mu_UB_EX_B_sp2
   - lambda_B[e]_sp2_1 >= 0
 S*lambda_1+mu_UB-mu_LB+psi=cEX_P_sp2_lb: - mu_LB_EX_P_sp2 + mu_UB_EX_P_sp2
   - lambda_P[e]_sp2_1 >= 0
 S*lambda_1+mu_UB-mu_LB+psi=cEX_h_sp2_lb: - mu_LB_EX_h_sp2 + mu_UB_EX_h_sp2
   - lambda_h[e]_sp2_1 >= 0
 S*lambda_1+mu_UB-mu_LB+psi=cEX_pi_sp2_lb: - mu_LB_EX_pi_sp2
   + mu_UB_EX_pi_sp2 - lambda_pi[e]_sp2_1 >= 0
 S*lambda_1+mu_UB-mu_LB+psi=cEX_h2o_sp2_lb: - mu_LB_EX_h2o_sp2
   + mu_UB_EX_h2o_sp2 - lambda_h2o[e]_sp2_1 >= 0
 S*lambda_1+mu_UB-mu_LB+psi=cT_S_sp2_lb: - mu_LB_T_S_sp2 + mu_UB_T_S_sp2
   - lambda_S[e]_sp2_1 + lambda_S[c]_sp2_1 >= 0
 S*lambda_1+mu_UB-mu_LB+psi=cT_A_sp2_lb: - mu_LB_T_A_sp2 + mu_UB_T_A_sp2
   - lambda_A[e]_sp2_1 + lambda_A[c]_sp2_1 >= 0
 S*lambda_1+mu_UB-mu_LB+psi=cT_B_sp2_lb: - mu_LB_T_B_sp2 + mu_UB_T_B_sp2
   - lambda_B[e]_sp2_1 + lambda_B[c]_sp2_1 >= 0
 S*lambda_1+mu_UB-mu_LB+psi=cT_h_sp2_lb: - mu_LB_T_h_sp2 + mu_UB_T_h_sp2
   - lambda_h[e]_sp2_1 + lambda_h[c]_sp2_1 >= 0
 S*lambda_1+mu_UB-mu_LB+psi=cT_P_sp2_lb: - mu_LB_T_P_sp2 + mu_UB_T_P_sp2
   - lambda_P[e]_sp2_1 + lambda_P[c]_sp2_1 >= 0
 S*lambda_1+mu_UB-mu_LB+psi=cT_pi_sp2_lb: - mu_LB_T_pi_sp2 + mu_UB_T_pi_sp2
   - lambda_pi[e]_sp2_1 + lambda_pi[c]_sp2_1 >= 0
 S*lambda_1+mu_UB-mu_LB+psi=cT_h2o_sp2_lb: - mu_LB_T_h2o_sp2
   + mu_UB_T_h2o_sp2 - lambda_h2o[e]_sp2_1 + lambda_h2o[c]_sp2_1 >= 0
 S*lambda_1+mu_UB-mu_LB+psi=cR_1_sp2_lb: - mu_LB_R_1_sp2 + mu_UB_R_1_sp2
   - lambda_S[c]_sp2_1 - 2 lambda_h[c]_sp2_1 + lambda_P[c]_sp2_1
   - 2 lambda_adp[c]_sp2_1 + 2 lambda_atp[c]_sp2_1 >= 0
 S*lambda_1+mu_UB-mu_LB+psi=cR_2_sp2_lb: - mu_LB_R_2_sp2 + mu_UB_R_2_sp2
   + lambda_B[c]_sp2_1 + 3 lambda_h[c]_sp2_1 - lambda_P[c]_sp2_1
   + 3 lambda_adp[c]_sp2_1 - 3 lambda_atp[c]_sp2_1 >= 0
 S*lambda_1+mu_UB-mu_LB+psi=cR_3_sp2_lb: - mu_LB_R_3_sp2 + mu_UB_R_3_sp2
   + lambda_A[c]_sp2_1 + lambda_h[c]_sp2_1 - lambda_P[c]_sp2_1
   + lambda_adp[c]_sp2_1 - lambda_atp[c]_sp2_1 >= 0
 S*lambda_1+mu_UB-mu_LB+psi=cR_4_sp2_lb: - mu_LB_R_4_sp2 + mu_UB_R_4_sp2
   + lambda_h[c]_sp2_1 + lambda_pi[c]_sp2_1 - lambda_h2o[c]_sp2_1
   + lambda_adp[c]_sp2_1 - lambda_atp[c]_sp2_1 >= 0
 S*lambda_1+mu_UB-mu_LB+psi=cOBJ_sp2_lb: - mu_LB_OBJ_sp2 + mu_UB_OBJ_sp2
   - 3 lambda_A[c]_sp2_1 - 3 lambda_B[c]_sp2_1 + 5 lambda_h[c]_sp2_1
   + 5 lambda_adp[c]_sp2_1 - 5 lambda_atp[c]_sp2_1
   + lambda_biomass[c]_sp2_1 >= 0
 S*lambda_1+mu_UB-mu_LB+psi=cBiomass_sp2_lb: - mu_LB_Biomass_sp2
   + mu_UB_Biomass_sp2 - lambda_biomass[c]_sp2_1 >= 1
 S*lambda_1+mu_UB-mu_LB+psi=cEX_A_2_sp1_lb: - mu_LB_EX_A_2_sp1
   + mu_UB_EX_A_2_sp1 - lambda_A[e]_sp1_1 + psi_EX_A_2_sp1 >= 0
 S*lambda_1+mu_UB-mu_LB+psi=cEX_B_2_sp1_lb: - mu_LB_EX_B_2_sp1
   + mu_UB_EX_B_2_sp1 - lambda_B[e]_sp1_1 + psi_EX_B_2_sp1 >= 0
 S*lambda_1+mu_UB-mu_LB+psi=cEX_P_2_sp1_lb: - mu_LB_EX_P_2_sp1
   + mu_UB_EX_P_2_sp1 - lambda_P[e]_sp1_1 >= 0
 S*lambda_1+mu_UB-mu_LB+psi=cEX_S_2_sp1_lb: - mu_LB_EX_S_2_sp1
   + mu_UB_EX_S_2_sp1 - lambda_S[e]_sp1_1 >= 0
 S*lambda_1+mu_UB-mu_LB+psi=cEX_h2o_2_sp1_lb: - mu_LB_EX_h2o_2_sp1
   + mu_UB_EX_h2o_2_sp1 - lambda_h2o[e]_sp1_1 >= 0
 S*lambda_1+mu_UB-mu_LB+psi=cEX_h_2_sp1_lb: - mu_LB_EX_h_2_sp1
   + mu_UB_EX_h_2_sp1 - lambda_h[e]_sp1_1 >= 0
 S*lambda_1+mu_UB-mu_LB+psi=cEX_pi_2_sp1_lb: - mu_LB_EX_pi_2_sp1
   + mu_UB_EX_pi_2_sp1 - lambda_pi[e]_sp1_1 >= 0
 S*lambda_1+mu_UB-mu_LB+psi=cEX_A_2_sp2_lb: - mu_LB_EX_A_2_sp2
   + mu_UB_EX_A_2_sp2 - lambda_A[e]_sp2_1 + psi_EX_A_2_sp2 >= 0
 S*lambda_1+mu_UB-mu_LB+psi=cEX_B_2_sp2_lb: - mu_LB_EX_B_2_sp2
   + mu_UB_EX_B_2_sp2 - lambda_B[e]_sp2_1 + psi_EX_B_2_sp2 >= 0
 S*lambda_1+mu_UB-mu_LB+psi=cEX_P_2_sp2_lb: - mu_LB_EX_P_2_sp2
   + mu_UB_EX_P_2_sp2 - lambda_P[e]_sp2_1 >= 0
 S*lambda_1+mu_UB-mu_LB+psi=cEX_S_2_sp2_lb: - mu_LB_EX_S_2_sp2
   + mu_UB_EX_S_2_sp2 - lambda_S[e]_sp2_1 >= 0
 S*lambda_1+mu_UB-mu_LB+psi=cEX_h2o_2_sp2_lb: - mu_LB_EX_h2o_2_sp2
   + mu_UB_EX_h2o_2_sp2 - lambda_h2o[e]_sp2_1 >= 0
 S*lambda_1+mu_UB-mu_LB+psi=cEX_h_2_sp2_lb: - mu_LB_EX_h_2_sp2
   + mu_UB_EX_h_2_sp2 - lambda_h[e]_sp2_1 >= 0
 S*lambda_1+mu_UB-mu_LB+psi=cEX_pi_2_sp2_lb: - mu_LB_EX_pi_2_sp2
   + mu_UB_EX_pi_2_sp2 - lambda_pi[e]_sp2_1 >= 0
 Xn*V_ut-ub<=0_EX_A_2_sp1_lb: 0.5 vut_EX_A_2_sp1 - ub_EX_A_2_sp1
   >= -1e+100
 Xn*V_ut-ub<=0_EX_B_2_sp1_lb: 0.5 vut_EX_B_2_sp1 - ub_EX_B_2_sp1
   >= -1e+100
 Xn*V_ut-ub<=0_EX_A_2_sp2_lb: 0.5 vut_EX_A_2_sp2 - ub_EX_A_2_sp2
   >= -1e+100
 Xn*V_ut-ub<=0_EX_B_2_sp2_lb: 0.5 vut_EX_B_2_sp2 - ub_EX_B_2_sp2
   >= -1e+100
 v_ut=max(0,v_oth)_1_EX_A_2_sp1_lb: - 0.5 v_1_EX_A_2_sp2
   - 1000 delta_EX_A_2_sp1 + ub_EX_A_2_sp1 >= -1e+100
 v_ut=max(0,v_oth)_1_EX_B_2_sp1_lb: - 0.5 v_1_EX_B_2_sp2
   - 1000 delta_EX_B_2_sp1 + ub_EX_B_2_sp1 >= -1e+100
 v_ut=max(0,v_oth)_1_EX_A_2_sp2_lb: - 0.5 v_1_EX_A_2_sp1
   - 1000 delta_EX_A_2_sp2 + ub_EX_A_2_sp2 >= -1e+100
 v_ut=max(0,v_oth)_1_EX_B_2_sp2_lb: - 0.5 v_1_EX_B_2_sp1
   - 1000 delta_EX_B_2_sp2 + ub_EX_B_2_sp2 >= -1e+100
 v_ut=max(0,v_oth)_2_EX_A_2_sp1_lb: 1000 delta_EX_A_2_sp1 + ub_EX_A_2_sp1
   >= -1e+100
 v_ut=max(0,v_oth)_2_EX_B_2_sp1_lb: 1000 delta_EX_B_2_sp1 + ub_EX_B_2_sp1
   >= -1e+100
 v_ut=max(0,v_oth)_2_EX_A_2_sp2_lb: 1000 delta_EX_A_2_sp2 + ub_EX_A_2_sp2
   >= -1e+100
 v_ut=max(0,v_oth)_2_EX_B_2_sp2_lb: 1000 delta_EX_B_2_sp2 + ub_EX_B_2_sp2
   >= -1e+100
 v_ut=max(0,v_oth)_3_EX_A_2_sp1_lb: - 0.5 v_1_EX_A_2_sp2 + ub_EX_A_2_sp1
   >= 0
 v_ut=max(0,v_oth)_3_EX_B_2_sp1_lb: - 0.5 v_1_EX_B_2_sp2 + ub_EX_B_2_sp1
   >= 0
 v_ut=max(0,v_oth)_3_EX_A_2_sp2_lb: - 0.5 v_1_EX_A_2_sp1 + ub_EX_A_2_sp2
   >= 0
 v_ut=max(0,v_oth)_3_EX_B_2_sp2_lb: - 0.5 v_1_EX_B_2_sp1 + ub_EX_B_2_sp2
   >= 0
 SV=0_S[e]_sp1_ub: - v_1_EX_S_sp1 - v_1_T_S_sp1 - v_1_EX_S_2_sp1 <= 0
 SV=0_A[e]_sp1_ub: - v_1_EX_A_sp1 - v_1_T_A_sp1 - v_1_EX_A_2_sp1 <= 0
 SV=0_B[e]_sp1_ub: - v_1_EX_B_sp1 - v_1_T_B_sp1 - v_1_EX_B_2_sp1 <= 0
 SV=0_P[e]_sp1_ub: - v_1_EX_P_sp1 - v_1_T_P_sp1 - v_1_EX_P_2_sp1 <= 0
 SV=0_h[e]_sp1_ub: - v_1_EX_h_sp1 - v_1_T_h_sp1 - v_1_EX_h_2_sp1 <= 0
 SV=0_pi[e]_sp1_ub: - v_1_EX_pi_sp1 - v_1_T_pi_sp1 - v_1_EX_pi_2_sp1 <= 0
 SV=0_h2o[e]_sp1_ub: - v_1_EX_h2o_sp1 - v_1_T_h2o_sp1 - v_1_EX_h2o_2_sp1
   <= 0
 SV=0_S[c]_sp1_ub: v_1_T_S_sp1 - v_1_R_1_sp1 <= 0
 SV=0_A[c]_sp1_ub: v_1_T_A_sp1 + v_1_R_3_sp1 - 3 v_1_OBJ_sp1 <= 0
 SV=0_B[c]_sp1_ub: v_1_T_B_sp1 + v_1_R_2_sp1 - 3 v_1_OBJ_sp1 <= 0
 SV=0_h[c]_sp1_ub: v_1_T_h_sp1 - 2 v_1_R_1_sp1 + v_1_R_2_sp1
   + 3 v_1_R_3_sp1 + v_1_R_4_sp1 + 5 v_1_OBJ_sp1 <= 0
 SV=0_P[c]_sp1_ub: v_1_T_P_sp1 + v_1_R_1_sp1 - v_1_R_2_sp1 - v_1_R_3_sp1
   <= 0
 SV=0_pi[c]_sp1_ub: v_1_T_pi_sp1 + v_1_R_4_sp1 <= 0
 SV=0_h2o[c]_sp1_ub: v_1_T_h2o_sp1 - v_1_R_4_sp1 <= 0
 SV=0_adp[c]_sp1_ub: - 2 v_1_R_1_sp1 + v_1_R_2_sp1 + 3 v_1_R_3_sp1
   + v_1_R_4_sp1 + 5 v_1_OBJ_sp1 <= 0
 SV=0_atp[c]_sp1_ub: 2 v_1_R_1_sp1 - v_1_R_2_sp1 - 3 v_1_R_3_sp1
   - v_1_R_4_sp1 - 5 v_1_OBJ_sp1 <= 0
 SV=0_biomass[c]_sp1_ub: v_1_OBJ_sp1 - v_1_Biomass_sp1 <= 0
 SV=0_S[e]_sp2_ub: - v_1_EX_S_sp2 - v_1_T_S_sp2 - v_1_EX_S_2_sp2 <= 0
 SV=0_A[e]_sp2_ub: - v_1_EX_A_sp2 - v_1_T_A_sp2 - v_1_EX_A_2_sp2 <= 0
 SV=0_B[e]_sp2_ub: - v_1_EX_B_sp2 - v_1_T_B_sp2 - v_1_EX_B_2_sp2 <= 0
 SV=0_P[e]_sp2_ub: - v_1_EX_P_sp2 - v_1_T_P_sp2 - v_1_EX_P_2_sp2 <= 0
 SV=0_h[e]_sp2_ub: - v_1_EX_h_sp2 - v_1_T_h_sp2 - v_1_EX_h_2_sp2 <= 0
 SV=0_pi[e]_sp2_ub: - v_1_EX_pi_sp2 - v_1_T_pi_sp2 - v_1_EX_pi_2_sp2 <= 0
 SV=0_h2o[e]_sp2_ub: - v_1_EX_h2o_sp2 - v_1_T_h2o_sp2 - v_1_EX_h2o_2_sp2
   <= 0
 SV=0_S[c]_sp2_ub: v_1_T_S_sp2 - v_1_R_1_sp2 <= 0
 SV=0_A[c]_sp2_ub: v_1_T_A_sp2 + v_1_R_3_sp2 - 3 v_1_OBJ_sp2 <= 0
 SV=0_B[c]_sp2_ub: v_1_T_B_sp2 + v_1_R_2_sp2 - 3 v_1_OBJ_sp2 <= 0
 SV=0_h[c]_sp2_ub: v_1_T_h_sp2 - 2 v_1_R_1_sp2 + 3 v_1_R_2_sp2
   + v_1_R_3_sp2 + v_1_R_4_sp2 + 5 v_1_OBJ_sp2 <= 0
 SV=0_P[c]_sp2_ub: v_1_T_P_sp2 + v_1_R_1_sp2 - v_1_R_2_sp2 - v_1_R_3_sp2
   <= 0
 SV=0_pi[c]_sp2_ub: v_1_T_pi_sp2 + v_1_R_4_sp2 <= 0
 SV=0_h2o[c]_sp2_ub: v_1_T_h2o_sp2 - v_1_R_4_sp2 <= 0
 SV=0_adp[c]_sp2_ub: - 2 v_1_R_1_sp2 + 3 v_1_R_2_sp2 + v_1_R_3_sp2
   + v_1_R_4_sp2 + 5 v_1_OBJ_sp2 <= 0
 SV=0_atp[c]_sp2_ub: 2 v_1_R_1_sp2 - 3 v_1_R_2_sp2 - v_1_R_3_sp2
   - v_1_R_4_sp2 - 5 v_1_OBJ_sp2 <= 0
 SV=0_biomass[c]_sp2_ub: v_1_OBJ_sp2 - v_1_Biomass_sp2 <= 0
 V-V_ut+V_ex=0_EX_A_2_sp1_ub: v_1_EX_A_2_sp1 + vut_EX_A_2_sp1
   - vex_EX_A_2_sp1 <= 0
 V-V_ut+V_ex=0_EX_B_2_sp1_ub: v_1_EX_B_2_sp1 + vut_EX_B_2_sp1
   - vex_EX_B_2_sp1 <= 0
 V-V_ut+V_ex=0_EX_A_2_sp2_ub: v_1_EX_A_2_sp2 + vut_EX_A_2_sp2
   - vex_EX_A_2_sp2 <= 0
 V-V_ut+V_ex=0_EX_B_2_sp2_ub: v_1_EX_B_2_sp2 + vut_EX_B_2_sp2
   - vex_EX_B_2_sp2 <= 0
 Xn*W+psi=0_EX_A_2_sp1_ub: psi_EX_A_2_sp1 + 0.5 omega_EX_A_2_sp1 <= 1e+100
 Xn*W+psi=0_EX_B_2_sp1_ub: psi_EX_B_2_sp1 + 0.5 omega_EX_B_2_sp1 <= 1e+100
 Xn*W+psi=0_EX_A_2_sp2_ub: psi_EX_A_2_sp2 + 0.5 omega_EX_A_2_sp2 <= 1e+100
 Xn*W+psi=0_EX_B_2_sp2_ub: psi_EX_B_2_sp2 + 0.5 omega_EX_B_2_sp2 <= 1e+100
 S*lambda_1+mu_UB-mu_LB+psi=cEX_S_sp1_ub: - mu_LB_EX_S_sp1 + mu_UB_EX_S_sp1
   - lambda_S[e]_sp1_1 <= 0
 S*lambda_1+mu_UB-mu_LB+psi=cEX_A_sp1_ub: - mu_LB_EX_A_sp1 + mu_UB_EX_A_sp1
   - lambda_A[e]_sp1_1 <= 0
 S*lambda_1+mu_UB-mu_LB+psi=cEX_B_sp1_ub: - mu_LB_EX_B_sp1 + mu_UB_EX_B_sp1
   - lambda_B[e]_sp1_1 <= 0
 S*lambda_1+mu_UB-mu_LB+psi=cEX_P_sp1_ub: - mu_LB_EX_P_sp1 + mu_UB_EX_P_sp1
   - lambda_P[e]_sp1_1 <= 0
 S*lambda_1+mu_UB-mu_LB+psi=cEX_h_sp1_ub: - mu_LB_EX_h_sp1 + mu_UB_EX_h_sp1
   - lambda_h[e]_sp1_1 <= 0
 S*lambda_1+mu_UB-mu_LB+psi=cEX_pi_sp1_ub: - mu_LB_EX_pi_sp1
   + mu_UB_EX_pi_sp1 - lambda_pi[e]_sp1_1 <= 0
 S*lambda_1+mu_UB-mu_LB+psi=cEX_h2o_sp1_ub: - mu_LB_EX_h2o_sp1
   + mu_UB_EX_h2o_sp1 - lambda_h2o[e]_sp1_1 <= 0
 S*lambda_1+mu_UB-mu_LB+psi=cT_S_sp1_ub: - mu_LB_T_S_sp1 + mu_UB_T_S_sp1
   - lambda_S[e]_sp1_1 + lambda_S[c]_sp1_1 <= 0
 S*lambda_1+mu_UB-mu_LB+psi=cT_A_sp1_ub: - mu_LB_T_A_sp1 + mu_UB_T_A_sp1
   - lambda_A[e]_sp1_1 + lambda_A[c]_sp1_1 <= 0
 S*lambda_1+mu_UB-mu_LB+psi=cT_B_sp1_ub: - mu_LB_T_B_sp1 + mu_UB_T_B_sp1
   - lambda_B[e]_sp1_1 + lambda_B[c]_sp1_1 <= 0
 S*lambda_1+mu_UB-mu_LB+psi=cT_h_sp1_ub: - mu_LB_T_h_sp1 + mu_UB_T_h_sp1
   - lambda_h[e]_sp1_1 + lambda_h[c]_sp1_1 <= 0
 S*lambda_1+mu_UB-mu_LB+psi=cT_P_sp1_ub: - mu_LB_T_P_sp1 + mu_UB_T_P_sp1
   - lambda_P[e]_sp1_1 + lambda_P[c]_sp1_1 <= 0
 S*lambda_1+mu_UB-mu_LB+psi=cT_pi_sp1_ub: - mu_LB_T_pi_sp1 + mu_UB_T_pi_sp1
   - lambda_pi[e]_sp1_1 + lambda_pi[c]_sp1_1 <= 0
 S*lambda_1+mu_UB-mu_LB+psi=cT_h2o_sp1_ub: - mu_LB_T_h2o_sp1
   + mu_UB_T_h2o_sp1 - lambda_h2o[e]_sp1_1 + lambda_h2o[c]_sp1_1 <= 0
 S*lambda_1+mu_UB-mu_LB+psi=cR_1_sp1_ub: - mu_LB_R_1_sp1 + mu_UB_R_1_sp1
   - lambda_S[c]_sp1_1 - 2 lambda_h[c]_sp1_1 + lambda_P[c]_sp1_1
   - 2 lambda_adp[c]_sp1_1 + 2 lambda_atp[c]_sp1_1 <= 0
 S*lambda_1+mu_UB-mu_LB+psi=cR_2_sp1_ub: - mu_LB_R_2_sp1 + mu_UB_R_2_sp1
   + lambda_B[c]_sp1_1 + lambda_h[c]_sp1_1 - lambda_P[c]_sp1_1
   + lambda_adp[c]_sp1_1 - lambda_atp[c]_sp1_1 <= 0
 S*lambda_1+mu_UB-mu_LB+psi=cR_3_sp1_ub: - mu_LB_R_3_sp1 + mu_UB_R_3_sp1
   + lambda_A[c]_sp1_1 + 3 lambda_h[c]_sp1_1 - lambda_P[c]_sp1_1
   + 3 lambda_adp[c]_sp1_1 - 3 lambda_atp[c]_sp1_1 <= 0
 S*lambda_1+mu_UB-mu_LB+psi=cR_4_sp1_ub: - mu_LB_R_4_sp1 + mu_UB_R_4_sp1
   + lambda_h[c]_sp1_1 + lambda_pi[c]_sp1_1 - lambda_h2o[c]_sp1_1
   + lambda_adp[c]_sp1_1 - lambda_atp[c]_sp1_1 <= 0
 S*lambda_1+mu_UB-mu_LB+psi=cOBJ_sp1_ub: - mu_LB_OBJ_sp1 + mu_UB_OBJ_sp1
   - 3 lambda_A[c]_sp1_1 - 3 lambda_B[c]_sp1_1 + 5 lambda_h[c]_sp1_1
   + 5 lambda_adp[c]_sp1_1 - 5 lambda_atp[c]_sp1_1
   + lambda_biomass[c]_sp1_1 <= 0
 S*lambda_1+mu_UB-mu_LB+psi=cBiomass_sp1_ub: - mu_LB_Biomass_sp1
   + mu_UB_Biomass_sp1 - lambda_biomass[c]_sp1_1 <= 1
 S*lambda_1+mu_UB-mu_LB+psi=cEX_S_sp2_ub: - mu_LB_EX_S_sp2 + mu_UB_EX_S_sp2
   - lambda_S[e]_sp2_1 <= 0
 S*lambda_1+mu_UB-mu_LB+psi=cEX_A_sp2_ub: - mu_LB_EX_A_sp2 + mu_UB_EX_A_sp2
   - lambda_A[e]_sp2_1 <= 0
 S*lambda_1+mu_UB-mu_LB+psi=cEX_B_sp2_ub: - mu_LB_EX_B_sp2 + mu_UB_EX_B_sp2
   - lambda_B[e]_sp2_1 <= 0
 S*lambda_1+mu_UB-mu_LB+psi=cEX_P_sp2_ub: - mu_LB_EX_P_sp2 + mu_UB_EX_P_sp2
   - lambda_P[e]_sp2_1 <= 0
 S*lambda_1+mu_UB-mu_LB+psi=cEX_h_sp2_ub: - mu_LB_EX_h_sp2 + mu_UB_EX_h_sp2
   - lambda_h[e]_sp2_1 <= 0
 S*lambda_1+mu_UB-mu_LB+psi=cEX_pi_sp2_ub: - mu_LB_EX_pi_sp2
   + mu_UB_EX_pi_sp2 - lambda_pi[e]_sp2_1 <= 0
 S*lambda_1+mu_UB-mu_LB+psi=cEX_h2o_sp2_ub: - mu_LB_EX_h2o_sp2
   + mu_UB_EX_h2o_sp2 - lambda_h2o[e]_sp2_1 <= 0
 S*lambda_1+mu_UB-mu_LB+psi=cT_S_sp2_ub: - mu_LB_T_S_sp2 + mu_UB_T_S_sp2
   - lambda_S[e]_sp2_1 + lambda_S[c]_sp2_1 <= 0
 S*lambda_1+mu_UB-mu_LB+psi=cT_A_sp2_ub: - mu_LB_T_A_sp2 + mu_UB_T_A_sp2
   - lambda_A[e]_sp2_1 + lambda_A[c]_sp2_1 <= 0
 S*lambda_1+mu_UB-mu_LB+psi=cT_B_sp2_ub: - mu_LB_T_B_sp2 + mu_UB_T_B_sp2
   - lambda_B[e]_sp2_1 + lambda_B[c]_sp2_1 <= 0
 S*lambda_1+mu_UB-mu_LB+psi=cT_h_sp2_ub: - mu_LB_T_h_sp2 + mu_UB_T_h_sp2
   - lambda_h[e]_sp2_1 + lambda_h[c]_sp2_1 <= 0
 S*lambda_1+mu_UB-mu_LB+psi=cT_P_sp2_ub: - mu_LB_T_P_sp2 + mu_UB_T_P_sp2
   - lambda_P[e]_sp2_1 + lambda_P[c]_sp2_1 <= 0
 S*lambda_1+mu_UB-mu_LB+psi=cT_pi_sp2_ub: - mu_LB_T_pi_sp2 + mu_UB_T_pi_sp2
   - lambda_pi[e]_sp2_1 + lambda_pi[c]_sp2_1 <= 0
 S*lambda_1+mu_UB-mu_LB+psi=cT_h2o_sp2_ub: - mu_LB_T_h2o_sp2
   + mu_UB_T_h2o_sp2 - lambda_h2o[e]_sp2_1 + lambda_h2o[c]_sp2_1 <= 0
 S*lambda_1+mu_UB-mu_LB+psi=cR_1_sp2_ub: - mu_LB_R_1_sp2 + mu_UB_R_1_sp2
   - lambda_S[c]_sp2_1 - 2 lambda_h[c]_sp2_1 + lambda_P[c]_sp2_1
   - 2 lambda_adp[c]_sp2_1 + 2 lambda_atp[c]_sp2_1 <= 0
 S*lambda_1+mu_UB-mu_LB+psi=cR_2_sp2_ub: - mu_LB_R_2_sp2 + mu_UB_R_2_sp2
   + lambda_B[c]_sp2_1 + 3 lambda_h[c]_sp2_1 - lambda_P[c]_sp2_1
   + 3 lambda_adp[c]_sp2_1 - 3 lambda_atp[c]_sp2_1 <= 0
 S*lambda_1+mu_UB-mu_LB+psi=cR_3_sp2_ub: - mu_LB_R_3_sp2 + mu_UB_R_3_sp2
   + lambda_A[c]_sp2_1 + lambda_h[c]_sp2_1 - lambda_P[c]_sp2_1
   + lambda_adp[c]_sp2_1 - lambda_atp[c]_sp2_1 <= 0
 S*lambda_1+mu_UB-mu_LB+psi=cR_4_sp2_ub: - mu_LB_R_4_sp2 + mu_UB_R_4_sp2
   + lambda_h[c]_sp2_1 + lambda_pi[c]_sp2_1 - lambda_h2o[c]_sp2_1
   + lambda_adp[c]_sp2_1 - lambda_atp[c]_sp2_1 <= 0
 S*lambda_1+mu_UB-mu_LB+psi=cOBJ_sp2_ub: - mu_LB_OBJ_sp2 + mu_UB_OBJ_sp2
   - 3 lambda_A[c]_sp2_1 - 3 lambda_B[c]_sp2_1 + 5 lambda_h[c]_sp2_1
   + 5 lambda_adp[c]_sp2_1 - 5 lambda_atp[c]_sp2_1
   + lambda_biomass[c]_sp2_1 <= 0
 S*lambda_1+mu_UB-mu_LB+psi=cBiomass_sp2_ub: - mu_LB_Biomass_sp2
   + mu_UB_Biomass_sp2 - lambda_biomass[c]_sp2_1 <= 1
 S*lambda_1+mu_UB-mu_LB+psi=cEX_A_2_sp1_ub: - mu_LB_EX_A_2_sp1
   + mu_UB_EX_A_2_sp1 - lambda_A[e]_sp1_1 + psi_EX_A_2_sp1 <= 0
 S*lambda_1+mu_UB-mu_LB+psi=cEX_B_2_sp1_ub: - mu_LB_EX_B_2_sp1
   + mu_UB_EX_B_2_sp1 - lambda_B[e]_sp1_1 + psi_EX_B_2_sp1 <= 0
 S*lambda_1+mu_UB-mu_LB+psi=cEX_P_2_sp1_ub: - mu_LB_EX_P_2_sp1
   + mu_UB_EX_P_2_sp1 - lambda_P[e]_sp1_1 <= 0
 S*lambda_1+mu_UB-mu_LB+psi=cEX_S_2_sp1_ub: - mu_LB_EX_S_2_sp1
   + mu_UB_EX_S_2_sp1 - lambda_S[e]_sp1_1 <= 0
 S*lambda_1+mu_UB-mu_LB+psi=cEX_h2o_2_sp1_ub: - mu_LB_EX_h2o_2_sp1
   + mu_UB_EX_h2o_2_sp1 - lambda_h2o[e]_sp1_1 <= 0
 S*lambda_1+mu_UB-mu_LB+psi=cEX_h_2_sp1_ub: - mu_LB_EX_h_2_sp1
   + mu_UB_EX_h_2_sp1 - lambda_h[e]_sp1_1 <= 0
 S*lambda_1+mu_UB-mu_LB+psi=cEX_pi_2_sp1_ub: - mu_LB_EX_pi_2_sp1
   + mu_UB_EX_pi_2_sp1 - lambda_pi[e]_sp1_1 <= 0
 S*lambda_1+mu_UB-mu_LB+psi=cEX_A_2_sp2_ub: - mu_LB_EX_A_2_sp2
   + mu_UB_EX_A_2_sp2 - lambda_A[e]_sp2_1 + psi_EX_A_2_sp2 <= 0
 S*lambda_1+mu_UB-mu_LB+psi=cEX_B_2_sp2_ub: - mu_LB_EX_B_2_sp2
   + mu_UB_EX_B_2_sp2 - lambda_B[e]_sp2_1 + psi_EX_B_2_sp2 <= 0
 S*lambda_1+mu_UB-mu_LB+psi=cEX_P_2_sp2_ub: - mu_LB_EX_P_2_sp2
   + mu_UB_EX_P_2_sp2 - lambda_P[e]_sp2_1 <= 0
 S*lambda_1+mu_UB-mu_LB+psi=cEX_S_2_sp2_ub: - mu_LB_EX_S_2_sp2
   + mu_UB_EX_S_2_sp2 - lambda_S[e]_sp2_1 <= 0
 S*lambda_1+mu_UB-mu_LB+psi=cEX_h2o_2_sp2_ub: - mu_LB_EX_h2o_2_sp2
   + mu_UB_EX_h2o_2_sp2 - lambda_h2o[e]_sp2_1 <= 0
 S*lambda_1+mu_UB-mu_LB+psi=cEX_h_2_sp2_ub: - mu_LB_EX_h_2_sp2
   + mu_UB_EX_h_2_sp2 - lambda_h[e]_sp2_1 <= 0
 S*lambda_1+mu_UB-mu_LB+psi=cEX_pi_2_sp2_ub: - mu_LB_EX_pi_2_sp2
   + mu_UB_EX_pi_2_sp2 - lambda_pi[e]_sp2_1 <= 0
 Xn*V_ut-ub<=0_EX_A_2_sp1_ub: 0.5 vut_EX_A_2_sp1 - ub_EX_A_2_sp1 <= 0
 Xn*V_ut-ub<=0_EX_B_2_sp1_ub: 0.5 vut_EX_B_2_sp1 - ub_EX_B_2_sp1 <= 0
 Xn*V_ut-ub<=0_EX_A_2_sp2_ub: 0.5 vut_EX_A_2_sp2 - ub_EX_A_2_sp2 <= 0
 Xn*V_ut-ub<=0_EX_B_2_sp2_ub: 0.5 vut_EX_B_2_sp2 - ub_EX_B_2_sp2 <= 0
 v_ut=max(0,v_oth)_1_EX_A_2_sp1_ub: - 0.5 v_1_EX_A_2_sp2
   - 1000 delta_EX_A_2_sp1 + ub_EX_A_2_sp1 <= 0
 v_ut=max(0,v_oth)_1_EX_B_2_sp1_ub: - 0.5 v_1_EX_B_2_sp2
   - 1000 delta_EX_B_2_sp1 + ub_EX_B_2_sp1 <= 0
 v_ut=max(0,v_oth)_1_EX_A_2_sp2_ub: - 0.5 v_1_EX_A_2_sp1
   - 1000 delta_EX_A_2_sp2 + ub_EX_A_2_sp2 <= 0
 v_ut=max(0,v_oth)_1_EX_B_2_sp2_ub: - 0.5 v_1_EX_B_2_sp1
   - 1000 delta_EX_B_2_sp2 + ub_EX_B_2_sp2 <= 0
 v_ut=max(0,v_oth)_2_EX_A_2_sp1_ub: 1000 delta_EX_A_2_sp1 + ub_EX_A_2_sp1
   <= 1000
 v_ut=max(0,v_oth)_2_EX_B_2_sp1_ub: 1000 delta_EX_B_2_sp1 + ub_EX_B_2_sp1
   <= 1000
 v_ut=max(0,v_oth)_2_EX_A_2_sp2_ub: 1000 delta_EX_A_2_sp2 + ub_EX_A_2_sp2
   <= 1000
 v_ut=max(0,v_oth)_2_EX_B_2_sp2_ub: 1000 delta_EX_B_2_sp2 + ub_EX_B_2_sp2
   <= 1000
 v_ut=max(0,v_oth)_3_EX_A_2_sp1_ub: - 0.5 v_1_EX_A_2_sp2 + ub_EX_A_2_sp1
   <= 1e+100
 v_ut=max(0,v_oth)_3_EX_B_2_sp1_ub: - 0.5 v_1_EX_B_2_sp2 + ub_EX_B_2_sp1
   <= 1e+100
 v_ut=max(0,v_oth)_3_EX_A_2_sp2_ub: - 0.5 v_1_EX_A_2_sp1 + ub_EX_A_2_sp2
   <= 1e+100
 v_ut=max(0,v_oth)_3_EX_B_2_sp2_ub: - 0.5 v_1_EX_B_2_sp1 + ub_EX_B_2_sp2
   <= 1e+100
 bilinear1: v_1_Biomass_sp1 - 10 mu_LB_EX_S_sp1 - 5 mu_LB_EX_A_2_sp1
   - 5 mu_LB_EX_B_2_sp1 + [ - omega_EX_A_2_sp1 * ub_EX_A_2_sp1
   - omega_EX_B_2_sp1 * ub_EX_B_2_sp1 ] = 0
 bilinear2: v_1_Biomass_sp2 - 10 mu_LB_EX_S_sp2 - 5 mu_LB_EX_A_2_sp2
   - 5 mu_LB_EX_B_2_sp2 + [ - omega_EX_A_2_sp2 * ub_EX_A_2_sp2
   - omega_EX_B_2_sp2 * ub_EX_B_2_sp2 ] = 0
Bounds
 v_1_EX_S_sp1 = -10
 v_1_EX_A_sp1 = 0
 v_1_EX_B_sp1 = 0
 v_1_EX_P_sp1 = 0
 v_1_EX_h_sp1 = -0
 v_1_EX_pi_sp1 = -998.3333333333334
 v_1_EX_h2o_sp1 = -1.666666666666666
 v_1_T_S_sp1 = 10
 v_1_T_A_sp1 = 5
 v_1_T_B_sp1 = -5
 v_1_T_h_sp1 = 0
 v_1_T_P_sp1 = 0
 v_1_T_pi_sp1 = -1.666666666666666
 v_1_T_h2o_sp1 = 1.666666666666666
 v_1_R_1_sp1 = 10
 v_1_R_2_sp1 = 10
 v_1_R_3_sp1 = 0
 v_1_R_4_sp1 = 1.666666666666666
 v_1_OBJ_sp1 = 1.666666666666667
 v_1_Biomass_sp1 = 1.666666666666667
 v_1_EX_S_sp2 = -10
 v_1_EX_A_sp2 = 0
 v_1_EX_B_sp2 = 0
 v_1_EX_P_sp2 = 0
 v_1_EX_h_sp2 = 0
 v_1_EX_pi_sp2 = -998.3333333333334
 v_1_EX_h2o_sp2 = -1.666666666666666
 v_1_T_S_sp2 = 10
 v_1_T_A_sp2 = -5
 v_1_T_B_sp2 = 5
 v_1_T_h_sp2 = -0
 v_1_T_P_sp2 = 0
 v_1_T_pi_sp2 = -1.666666666666666
 v_1_T_h2o_sp2 = 1.666666666666666
 v_1_R_1_sp2 = 10
 v_1_R_2_sp2 = 0
 v_1_R_3_sp2 = 10
 v_1_R_4_sp2 = 1.666666666666666
 v_1_OBJ_sp2 = 1.666666666666667
 v_1_Biomass_sp2 = 1.666666666666667
 v_1_EX_A_2_sp1 = -5
 v_1_EX_B_2_sp1 = 5
 v_1_EX_P_2_sp1 = -0
 v_1_EX_S_2_sp1 = 0
 v_1_EX_h2o_2_sp1 = 0
 v_1_EX_h_2_sp1 = 0
 v_1_EX_pi_2_sp1 = 1000
 v_1_EX_A_2_sp2 = 5
 v_1_EX_B_2_sp2 = -5
 v_1_EX_P_2_sp2 = -0
 v_1_EX_S_2_sp2 = 0
 v_1_EX_h2o_2_sp2 = 0
 v_1_EX_h_2_sp2 = 0
 v_1_EX_pi_2_sp2 = 1000
 mu_LB_EX_h_sp1 = 0
 mu_LB_EX_pi_sp1 = 0
 mu_LB_EX_h2o_sp1 = 0
 mu_LB_T_S_sp1 = 0
 mu_LB_T_A_sp1 = 0
 mu_LB_T_B_sp1 = 0
 mu_LB_T_h_sp1 = 0
 mu_LB_T_P_sp1 = 0
 mu_LB_T_pi_sp1 = 0
 mu_LB_T_h2o_sp1 = 0
 mu_LB_EX_h_sp2 = 0
 mu_LB_EX_pi_sp2 = 0
 mu_LB_EX_h2o_sp2 = 0
 mu_LB_T_S_sp2 = 0
 mu_LB_T_A_sp2 = 0
 mu_LB_T_B_sp2 = 0
 mu_LB_T_h_sp2 = 0
 mu_LB_T_P_sp2 = 0
 mu_LB_T_pi_sp2 = 0
 mu_LB_T_h2o_sp2 = 0
 mu_UB_EX_h_sp1 = 0
 mu_UB_EX_pi_sp1 = 0
 mu_UB_EX_h2o_sp1 = 0
 mu_UB_T_S_sp1 = 0
 mu_UB_T_A_sp1 = 0
 mu_UB_T_B_sp1 = 0
 mu_UB_T_h_sp1 = 0
 mu_UB_T_P_sp1 = 0
 mu_UB_T_pi_sp1 = 0
 mu_UB_T_h2o_sp1 = 0
 mu_UB_R_1_sp1 = 0
 mu_UB_R_2_sp1 = 0
 mu_UB_R_3_sp1 = 0
 mu_UB_R_4_sp1 = 0
 mu_UB_OBJ_sp1 = 0
 mu_UB_Biomass_sp1 = 0
 mu_UB_EX_h_sp2 = 0
 mu_UB_EX_pi_sp2 = 0
 mu_UB_EX_h2o_sp2 = 0
 mu_UB_T_S_sp2 = 0
 mu_UB_T_A_sp2 = 0
 mu_UB_T_B_sp2 = 0
 mu_UB_T_h_sp2 = 0
 mu_UB_T_P_sp2 = 0
 mu_UB_T_pi_sp2 = 0
 mu_UB_T_h2o_sp2 = 0
 mu_UB_R_1_sp2 = 0
 mu_UB_R_2_sp2 = 0
 mu_UB_R_3_sp2 = 0
 mu_UB_R_4_sp2 = 0
 mu_UB_OBJ_sp2 = 0
 mu_UB_Biomass_sp2 = 0
 mu_UB_EX_A_2_sp1 = 0
 mu_UB_EX_B_2_sp1 = 0
 mu_UB_EX_P_2_sp1 = 0
 mu_UB_EX_h2o_2_sp1 = 0
 mu_UB_EX_h_2_sp1 = 0
 mu_UB_EX_pi_2_sp1 = 0
 mu_UB_EX_A_2_sp2 = 0
 mu_UB_EX_B_2_sp2 = 0
 mu_UB_EX_P_2_sp2 = 0
 mu_UB_EX_h2o_2_sp2 = 0
 mu_UB_EX_h_2_sp2 = 0
 mu_UB_EX_pi_2_sp2 = 0
 lambda_S[e]_sp1_1 free
 lambda_A[e]_sp1_1 free
 lambda_B[e]_sp1_1 free
 lambda_P[e]_sp1_1 free
 lambda_h[e]_sp1_1 free
 lambda_pi[e]_sp1_1 free
 lambda_h2o[e]_sp1_1 free
 lambda_S[c]_sp1_1 free
 lambda_A[c]_sp1_1 free
 lambda_B[c]_sp1_1 free
 lambda_h[c]_sp1_1 free
 lambda_P[c]_sp1_1 free
 lambda_pi[c]_sp1_1 free
 lambda_h2o[c]_sp1_1 free
 lambda_adp[c]_sp1_1 free
 lambda_atp[c]_sp1_1 free
 lambda_biomass[c]_sp1_1 free
 lambda_S[e]_sp2_1 free
 lambda_A[e]_sp2_1 free
 lambda_B[e]_sp2_1 free
 lambda_P[e]_sp2_1 free
 lambda_h[e]_sp2_1 free
 lambda_pi[e]_sp2_1 free
 lambda_h2o[e]_sp2_1 free
 lambda_S[c]_sp2_1 free
 lambda_A[c]_sp2_1 free
 lambda_B[c]_sp2_1 free
 lambda_h[c]_sp2_1 free
 lambda_P[c]_sp2_1 free
 lambda_pi[c]_sp2_1 free
 lambda_h2o[c]_sp2_1 free
 lambda_adp[c]_sp2_1 free
 lambda_atp[c]_sp2_1 free
 lambda_biomass[c]_sp2_1 free
 -infinity <= psi_EX_A_2_sp1 <= 0
 -infinity <= psi_EX_B_2_sp1 <= 0
 -infinity <= psi_EX_A_2_sp2 <= 0
 -infinity <= psi_EX_B_2_sp2 <= 0
Binaries
 delta_EX_A_2_sp1 delta_EX_B_2_sp1 delta_EX_A_2_sp2 delta_EX_B_2_sp2
End
