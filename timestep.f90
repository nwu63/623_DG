subroutine timeIntegration(q,p,I2E,B2E,In,Bn,Jinv,detJ,Minv,rBC,resids,resnorm,&
    gamma,Rgas,CFL,convtol,min_iter,max_iter,nelem,niface,nbface)
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
    integer, intent(in) :: nelem,niface,nbface,min_iter,max_iter,p
    real(8), intent(inout), dimension(nelem,(p+1)*(p+2)/2,4) :: q
    integer, intent(in), dimension(niface,4) :: I2E
    integer, intent(in), dimension(nbface,3) :: B2E
    real(8), intent(in), dimension(niface,3) :: In ! includes length as 3rd component
    real(8), intent(in), dimension(nbface,3) :: Bn ! includes length as 3rd component
    real(8), intent(in), dimension(nelem) :: detJ
    real(8), intent(in), dimension(nelem,2,2) :: Jinv
    real(8), intent(in), dimension(nelem,(p+1)*(p+2)/2,(p+1)*(p+2)/2) :: Minv
    real(8), intent(out),dimension(nelem,(p+1)*(p+2)/2,4) :: resids
    real(8), intent(in), dimension(5) :: rBC
    real(8), intent(in) :: gamma, Rgas, convtol, CFL
    real(8), intent(out), dimension(max_iter) :: resnorm
!f2py intent(in) node,E2N,I2E,B2E,In,Bn,gamma,rBC,Jinv,Minv
!f2py intent(out) resids, resnorm
!f2py intent(in,out) q
    integer :: ielem, iter, iq
    integer, dimension(3) :: loc
    real(8), dimension(nelem) :: wavespeed, temp_wavespeed ! temp is used b/c we want to use same wavespeed from FE step
    real(8), dimension(nelem,(p+1)*(p+2)/2,4) :: qFE
    real(8) :: dt
    resnorm(:) = -1 ! set to high value to allow first pass of while loop
    
    do iter = 1,max_iter
        call getResidual(q,p,I2E,B2E,In,Bn,rBC,resids,Jinv,detJ,wavespeed,gamma,Rgas,nelem,niface,nbface)
        resnorm(iter) = maxval(resids)
        loc = maxloc(resids)
        if (mod(iter,100) == 0) then
            print*, iter, resnorm(iter), loc(1), loc(2), loc(3)
        endif
        
        do ielem = 1,nelem
            dt = 2*CFL/wavespeed(ielem)
            do iq = 1,4
                qFE(ielem,:,iq) = q(ielem,:,iq) - dt*matmul(Minv(ielem,:,:),resids(ielem,:,iq)) ! the area term cancels with the dt expression
            enddo
        enddo
        if (p >= 0) then ! for now do FE for everything
            q = qFE
        elseif (p >= 1) then
            call getResidual(q,p,I2E,B2E,In,Bn,rBC,resids,Jinv,detJ,wavespeed,gamma,Rgas,nelem,niface,nbface)
            do ielem = 1,nelem
                dt = 2*CFL/wavespeed(ielem)
                q(ielem,:,:) = 0.5d0*(q(ielem,:,:) + qFE(ielem,:,:) - dt*resids(ielem,:,:))
            enddo
        endif
        if ((iter > max_iter) .or. (iter > min_iter .and. resnorm(iter) < convtol) .or. resnorm(iter) /= resnorm(iter)) then
            exit
        endif
    enddo
    print*, "Converged! Took ", iter, " iterations to reduce residual to ",resnorm(iter)
end subroutine timeIntegration