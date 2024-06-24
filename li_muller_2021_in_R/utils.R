logLik_GLR <- function(betahat_long, betahat_short, sig_biv, r_d, CV_struct) {
  
  # Extract CV
  CV_vals <- CV_struct$CV.vals
  angle_vec <- t(as.matrix(CV_struct$angle.vec))
  lograd_vec <- t(as.matrix(CV_struct$lograd.vec))
  
  # Compute sig
  sig_L <- sqrt(sig_biv$sigL2)
  sig_LS <- sig_biv$sigLS
  sig_S <- sqrt(sig_biv$sigS2)
  
  rho_LS <- sig_LS / (sig_S * sig_L)
  
  # 
  tauhat_b <- betahat_long / sig_L
  a1 <- (rho_LS^(-2) - 1)^(-1/2) * (1 / sig_L)
  a2 <- -(rho_LS^(-2) - 1)^(-1/2) * (1 / (sig_S * rho_LS))
  tauhat_theta <- a1 * betahat_long + a2 * betahat_short
  psi_b <- sig_L * (rho_LS^(-2) - 1)^(-1/2) * (1 / sig_L - 1 / (rho_LS * sig_S))
  max_tau_g <- abs( (rho_LS^(-2) - 1)^(-1/2) * (1 / (rho_LS * sig_S)) * r_d )
  
  # Compute angle and lograd
  angle <- atan(max_tau_g / abs(psi_b)) / (0.5 * pi)
  lograd <- 0.5 * ( 0.5 * log( abs(psi_b)^2 + max_tau_g^2 )/log(200) + 1 )
  
  # Find nearest idx
  temp <- findnearest(angle, angle_vec)
  CV_ind_angle1 <- temp$outlower
  CV_ind_angle2 <- temp$outupper
  temp <- findnearest(lograd, lograd_vec)
  CV_ind_lograd1 <- temp$outlower
  CV_ind_lograd2 <- temp$outupper
  
  # 
  out_logLik_GLR <- matrix(0, nrow = length(betahat_long), ncol = 1)
  
  # Calculate log likelihood
  ii <- 1
  
  temp <- logLik_GLR_norm(tauhat_b[ii], tauhat_theta[ii], psi_b[ii], max_tau_g[ii])
  wgt_angle <- ( angle[ii] - angle_vec[CV_ind_angle2[ii]] ) / ( angle_vec[CV_ind_angle1[ii]] - angle_vec[CV_ind_angle2[ii]] )
  wgt_lograd <- ( lograd[ii] - lograd_vec[CV_ind_lograd2[ii]] ) / ( lograd_vec[CV_ind_lograd1[ii]] - lograd_vec[CV_ind_lograd2[ii]] )
  wgt_angle <- min(abs(wgt_angle), 1)
  wgt_lograd <- min(abs(wgt_lograd), 1)
  
  CV_temp <- wgt_angle * wgt_lograd * CV_vals[CV_ind_lograd1[ii], CV_ind_angle1[ii]] +
    (1 - wgt_angle) * wgt_lograd * CV_vals[CV_ind_lograd1[ii], CV_ind_angle2[ii]] +
    wgt_angle * (1 - wgt_lograd) * CV_vals[CV_ind_lograd2[ii], CV_ind_angle1[ii]] +
    (1 - wgt_angle) * (1 - wgt_lograd) * CV_vals[CV_ind_lograd2[ii], CV_ind_angle2[ii]]
  
  out_logLik_GLR[ii] <- (temp > CV_temp)
  
  return(out_logLik_GLR)
}


findnearest <- function(A, B) {
  
  nA <- length(A) # Coz scalar, should be 1
  nB <- length(B) # Coz vector, can use length
  
  AB <- rbind(A, B)
  indAB <- rbind(t(1:nA), matrix(-1, nB, 1))
  indsortAB <- order(AB)
  indAB <- indAB[indsortAB]
  modAB <- indAB
  modAB[modAB > 0] <- 0
  
  modABlower <- -cumsum(modAB)
  modABlower <- pmax(modABlower, 1)
  
  outlower <- modABlower[indAB > 0]
  indsortB <- order(B)
  outlower <- indsortB[outlower]
  indsortA <- order(A)
  indsortA_inv <- order(indsortA)
  
  outlower <- outlower[indsortA_inv]
  modABupper <- -cumsum(modAB) + 1
  modABupper <- pmin(modABupper, nB)
  outupper <- modABupper[indAB > 0]
  outupper <- indsortB[outupper]
  
  outupper <- outupper[indsortA_inv]
  
  return(list(outlower = outlower, outupper = outupper))
}



logLik_GLR_norm <- function(tauhat_b, tauhat_theta, psi_b, max_tau_g) {
  
  v_tau <- rbind(1, psi_b)
  P_tau <- v_tau %*% solve(t(v_tau) %*% v_tau) %*% t(v_tau)
  M_tau <- diag(2) - P_tau
  tauhat <- rbind(tauhat_b, tauhat_theta)
  
  # 
  right <- rbind(0, max_tau_g)
  left <- rbind(0, -max_tau_g)
  right_resid <- M_tau %*% ( tauhat - right[, rep(1:ncol(right), ncol(tauhat))] )
  left_resid <- M_tau %*% ( tauhat - left[, rep(1:ncol(left), ncol(tauhat))] )
  
  GlogLik_1right <- colSums(right_resid^2)
  GlogLik_1left <- colSums(left_resid^2)
  
  GlogLik_1 <- -(colSums(left_resid * right_resid) >= 0) * min(GlogLik_1right, GlogLik_1left)
  right_diff <- rbind(tauhat_b, tauhat_theta) - right
  left_diff <- rbind(tauhat_b, tauhat_theta) - left
  GlogLik_0 <- -(tauhat_theta > max_tau_g) * colSums(right_diff^2) - 
    (tauhat_theta < -max_tau_g) * colSums(left_diff^2) - 
    ( (-max_tau_g <= tauhat_theta) & (tauhat_theta <= max_tau_g) ) * (tauhat_b^2)
  
  out_logLik_GLR_norm <- GlogLik_1 - GlogLik_0
  
  return(out_logLik_GLR_norm)
}


