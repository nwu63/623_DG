! TODO: add fail flag
! TODO: add applyMatRes subroutine
! TODO: check iteration counting in print out
! TODO: try dgemm instead of matmul
subroutine timeIntegration(q,p,I2E,B2E,In,Bn,qnrm,Jinv,Jinv2,detJ,detJ2,Minv,xy,w,gphi,w1,rBC,resids,&
    phiL,phiR,xyL,xyR,qlist,resnorm,gamma,Rgas,CFL,convtol,min_iter,max_iter,nelem,niface,nbface,nqelem,Ng,Ng1)
    ! -----------------------------------------------------------------------
    ! Purpose: use forward Euler to timestep the governing equations
    ! 
    ! Inputs:
    !   q[nelem,4] = the initial state
    !   I2E[niface,4]
    !   B2E[nbface,3]
    !   In[niface,3] = normal vector for interior faces, and length
    !   Bn[nbface,3] = normal vector for boundary faces, and length
    !   rBC[5] = information on freestream/boundary values
    !   gamma = heat ratio
    !   Rgas = specific gas constant
    !   dt = the delta t for each timestep
    !   ndt = number of time steps
    ! 
    ! Outs:
    !   q[nelem,4] = the state after ndt timesteps
    !   resids[nelem,4] = the residual after ndt timesteps
    ! 
    ! -----------------------------------------------------------------------
    implicit none
    integer, intent(in) :: nelem,niface,nbface,min_iter,max_iter,p,nqelem,Ng,Ng1
    real(8), intent(inout), dimension(nelem,(p+1)*(p+2)/2,4) :: q
    integer, intent(in), dimension(niface,4) :: I2E
    integer, intent(in), dimension(nbface,3) :: B2E
    real(8), intent(in), dimension(niface,3) :: In ! includes length as 3rd component
    real(8), intent(in), dimension(nbface,3) :: Bn ! includes length as 3rd component
    real(8), intent(in), dimension(nqelem,Ng1,2) :: qnrm ! nrm for high-q elems
    real(8), intent(in), dimension(nelem) :: detJ
    real(8), intent(in), dimension(nqelem,Ng) :: detJ2
    real(8), intent(in), dimension(nelem,2,2) :: Jinv
    real(8), intent(in), dimension(nqelem,Ng,2,2) :: Jinv2
    real(8), intent(in), dimension(Ng,(p+1)*(p+2)/2,2) :: gphi
    real(8), intent(in), dimension(Ng,2) :: xy
    real(8), intent(in), dimension(Ng) :: w
    real(8), intent(in), dimension(Ng1) :: w1
    real(8), intent(in), dimension(3,Ng1,(p+1)*(p+2)/2) :: phiL,phiR
    real(8), intent(in), dimension(3,Ng1,2) :: xyL,xyR
    integer, intent(in), dimension(nqelem) :: qlist
    real(8), intent(in), dimension(nelem,(p+1)*(p+2)/2,(p+1)*(p+2)/2) :: Minv
    real(8), intent(out),dimension(nelem,(p+1)*(p+2)/2,4) :: resids
    real(8), intent(in), dimension(5) :: rBC
    real(8), intent(in) :: gamma, Rgas, convtol, CFL
    real(8), intent(out), dimension(max_iter) :: resnorm
!f2py intent(in) node,E2N,I2E,B2E,In,Bn,gamma,rBC,Jinv,Minv
!f2py intent(out) resids, resnorm
!f2py intent(in,out) q
    integer :: ielem, iter
    integer, dimension(3) :: loc
    real(8), dimension(nelem) :: wavespeed
    real(8), dimension(nelem,(p+1)*(p+2)/2,4) :: q1,q2
    real(8), dimension((p+1)*(p+2)/2,4) :: new_resids
    real(8) :: dt
    resnorm(:) = -1 ! set to high value to allow first pass of while loop
    
    do iter = 1,max_iter
        call getResidual(q,p,I2E,B2E,In,Bn,qnrm,rBC,resids,Jinv,Jinv2,detJ,detJ2,xy,w,gphi,w1,&
            phiL,phiR,xyL,xyR,qlist,wavespeed,gamma,Rgas,nelem,niface,nbface,nqelem,Ng,Ng1)
        resnorm(iter) = maxval(resids)
        loc = maxloc(resids)
        if (mod(iter,100) == 0) then
            print*, iter, resnorm(iter), loc(1), loc(2), loc(3)
        endif
        
        do ielem = 1,nelem
            call applyMatRes(Minv(ielem,:,:),resids(ielem,:,:),new_resids,p)
            dt = CFL*detJ(ielem)/wavespeed(ielem) ! we let detJ = 2*area
            q1(ielem,:,:) = q(ielem,:,:) - dt*new_resids ! the area term cancels with the dt expression
        enddo
        if (p == 0) then 
            q = q1
        elseif (p == 1) then
            call getResidual(q1,p,I2E,B2E,In,Bn,qnrm,rBC,resids,Jinv,Jinv2,detJ,detJ2,xy,w,gphi,w1,&
                phiL,phiR,xyL,xyR,qlist,wavespeed,gamma,Rgas,nelem,niface,nbface,nqelem,Ng,Ng1)
            do ielem = 1,nelem
                call applyMatRes(Minv(ielem,:,:),resids(ielem,:,:),new_resids,p)
                dt = CFL*detJ(ielem)/wavespeed(ielem)
                q(ielem,:,:) = 0.5d0*(q(ielem,:,:) + q1(ielem,:,:) - dt*new_resids)
            enddo
        elseif (p == 2) then
            call getResidual(q1,p,I2E,B2E,In,Bn,qnrm,rBC,resids,Jinv,Jinv2,detJ,detJ2,xy,w,gphi,w1,&
                phiL,phiR,xyL,xyR,qlist,wavespeed,gamma,Rgas,nelem,niface,nbface,nqelem,Ng,Ng1)
            do ielem = 1,nelem
                call applyMatRes(Minv(ielem,:,:),resids(ielem,:,:),new_resids,p)
                dt = CFL*detJ(ielem)/wavespeed(ielem)
                q2(ielem,:,:) = 0.25d0*(3.d0*q(ielem,:,:) + q1(ielem,:,:) - dt*new_resids)
            enddo
            call getResidual(q2,p,I2E,B2E,In,Bn,qnrm,rBC,resids,Jinv,Jinv2,detJ,detJ2,xy,w,gphi,w1,&
                phiL,phiR,xyL,xyR,qlist,wavespeed,gamma,Rgas,nelem,niface,nbface,nqelem,Ng,Ng1)
            
            do ielem = 1,nelem
                call applyMatRes(Minv(ielem,:,:),resids(ielem,:,:),new_resids,p)
                dt = CFL*detJ(ielem)/wavespeed(ielem)
                q(ielem,:,:) = 1.d0/3.d0*(q(ielem,:,:) + 2.d0*q2(ielem,:,:) - 2.d0*dt*new_resids)
            enddo
        endif
        if ((iter > min_iter .and. resnorm(iter) < convtol) .or. resnorm(iter) /= resnorm(iter)) then
            exit
        endif
    enddo
    print*, "Converged! Took ", iter-1, " iterations to reduce residual to ",resnorm(iter)
end subroutine timeIntegration

subroutine applyMatRes(Minv,resids,new_resids,p)
    implicit none
    integer, intent(in) :: p
    real(8), intent(in), dimension((p+1)*(p+2)/2,(p+1)*(p+2)/2) :: Minv
    real(8), intent(in),dimension((p+1)*(p+2)/2,4) :: resids
    real(8), intent(out),dimension((p+1)*(p+2)/2,4) :: new_resids
    integer :: Nb
    Nb = (p+1)*(p+2)/2
    call dgemm('N','N',Nb,4,Nb,1.d0,Minv,Nb,resids,Nb,0.d0,new_resids,Nb)
    ! new_resids = matmul(Minv,resids)
end subroutine applyMatRes