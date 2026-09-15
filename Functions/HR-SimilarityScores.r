ms_plot <- function(zA,rA,title=NULL){
  
  min_mz = min(zA)-5
  max_mz = max(zA)+5 
  
  plot(zA,rA,type="h",col="black",
       xlim=c(min_mz,max_mz),
       ylim=c(0,1),
       main=title,
       xlab="mass-to-charge ratio",
       ylab="relative intensity",
       yaxt="n")
  axis(2,at=seq(0,1,0.5),labels=abs(seq(0,1,0.5)),las=2)
}

HT_plot <- function(zA,rA,zB,rB, m1 = NULL, m2=NULL, title=NULL,score=NULL){
  min_mz = min(zA,zB)
  max_mz = max(zA,zB)
  
  if(is.null(score)==FALSE){
    title = paste0(title,"\n(",score,")")
  }
  
  cA = rep("black",length(zA))
  cB = rep("black",length(zB))
  
  if(is.null(m1)==FALSE){
    i = which(zA<m1)
    if(length(i)>0) cA[i] = "red"
    i = which(zB<m1)
    if(length(i)>0) cB[i] = "red"
  }
  
  if(is.null(m2)==FALSE){
    i = which(zA>m2)
    if(length(i)>0) cA[i] = "red"
    i = which(zB>m2)
    if(length(i)>0) cB[i] = "red"
  }
  
  plot(zA,rA,type="h",col=cA,
       xlim=c(min_mz,max_mz),
       ylim=c(-1,1),
       main=title,
       xlab="mass-to-charge ratio",
       ylab="relative intensity",
       yaxt="n")
  segments(zB,rep(0,length(zB)),zB,-rB,col=cB)  
  axis(2,at=seq(-1,1,0.5),labels=abs(seq(-1,1,0.5)),las=2)

}

spec_processing <-function(zA,rA,tau,epsilon){
  
  relevant_A = which(rA>=tau)
  
  bar_zA = zA[relevant_A]
  bar_rA = rA[relevant_A]
  
  hat_zA = round(bar_zA,(epsilon+1)) 
  
  star_zA = unique(hat_zA)
  star_rA = numeric(length(star_zA))
  
  for(i in 1:length(star_zA)){
    lambda = which(hat_zA==star_zA[i])
    star_rA[i] = sum(bar_rA[lambda])
  }
  
  return(list(mz = star_zA,
              intensity = star_rA))
  
}

computeSim <- function(zA,rA,zB,rB,epsilon,m1=NULL,m2=NULL,wL=0.33,wM=0.33,wH=0.33){
  
  mz_tol = 10^(-epsilon)
  u_z = sort(unique(c(zA,zB)))
  a = which(c(u_z,0)-c(0,u_z)<=mz_tol)
  a = a[-length(a)]
  if(length(a)>0) {
    u_z = u_z[-a]
  }
  
  lambda = length(u_z)
  
  # Generate Relevant Vectors
  rAlpha = numeric(lambda)
  rBeta = numeric(lambda)
  
  for(i in 1:lambda){
    
    beta = which(abs(zB - u_z[i])<=mz_tol)
    if(length(beta)==0){ 
      rBeta[i] = 0
    } else {
      rBeta[i] = sum(rB[beta])    
    }
    
    alpha = which(abs(zA - u_z[i])<=mz_tol)
    if(length(alpha)==0){ 
      rAlpha[i] = 0
    } else {
      rAlpha[i] = sum(rA[alpha])    
    }
    
  }
  
  ## Compute cosine similarity between vectors
  if(is.null(m1)==FALSE){
    if(is.null(m2==FALSE)){
      ## Three regions
      csl = 0
      csm = 0
      csh = 0
      i = which(u_z<m1)
      if(length(i)>0){
        csl = wL * cosSim(rAlpha[i],rBeta[i])
      } 
      
      i = which(u_z>m2)
      if(length(i)>0){
        csh = wH * cosSim(rAlpha[i],rBeta[i])
      } 
      
      i = intersect(which(u_z>=m1),which(u_z<=m2))
      if(length(i)>0){
        csm = wM * cosSim(rAlpha[i],rBeta[i])
      }
      
      sim = csl+csm+csh
      
    } else {
      ## Two regions (low or mid/high)
      csl = 0
      csmh = 0

      i = which(u_z<m1)
      if(length(i)>0){
        csl = wL * cosSim(rAlpha[i],rBeta[i])
      } 
      
      i = which(u_z>=m1)
      if(length(i)>0){
        csmh = (wM+wH) * cosSim(rAlpha[i],rBeta[i])
      } 

      
      sim = csl+csmh
    }
  } else {
    if(is.null(m2==FALSE)){
      ## Two regions (low/mid or high)
      cslm = 0
      csh = 0
      
      i = which(u_z<m2)
      if(length(i)>0){
        cslm = (wL+wM) * cosSim(rAlpha[i],rBeta[i])
      } 
      
      i = which(u_z>=m2)
      if(length(i)>0){
        csh = wJ * cosSim(rAlpha[i],rBeta[i])
      } 
      sim = cslm+csh
    } else {
      ## One region
      sim = cosSim(rAlpha,rBeta)
    }
  }
  
return(sim)
  
  
}

cosSim <- function(rAlpha,rBeta){
  N = sum(rAlpha*rBeta)
  DA = sqrt(sum(rAlpha*rAlpha)) 
  DB = sqrt(sum(rBeta*rBeta)) 
  cs = N/(DA*DB + 1e-8)
  
  return(cs)
  
}
