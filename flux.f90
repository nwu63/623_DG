subroutine eulerFluxVec(q,F,gamma)
    ! -----------------------------------------------------------------------
    ! Purpose: Calculates the analytical Euler flux given a state and normal
    ! 
    ! Inputs:
    !   q[4] = the state
    !   nrm[2] = the normal vector
    !   gamma = heat ratio
    ! 
    ! Outs:
    !   F[4] = analytical Euler flux
    ! 
    ! -----------------------------------------------------------------------
    implicit none
    real(8), intent(in), dimension(0:3) :: q
    real(8), intent(in) :: gamma
    real(8), intent(out), dimension(0:3,0:1) :: F
!f2py intent(in) q, nrm,gamma
!f2py intent(out) F
    real(8) :: velnorm,p,H

    velnorm = sqrt((q(1)/q(0))**2.d0 + (q(2)/q(0))**2.d0) ! this is the norm of velocity vector
    p = (gamma-1)*(q(3) - 0.5d0*q(0)*velnorm**2.d0)
    H = q(3)/q(0) + p/q(0)
    ! ---------- construct the entire flux vector
    F(0,0) = q(1)
    F(1,0) = q(1)**2.d0/q(0) + p
    F(2,0) = q(1)*q(2)/q(0)
    F(3,0) = q(1)*H

    F(0,1) = q(2)
    F(1,1) = q(1)*q(2)/q(0)
    F(2,1) = q(2)**2.d0/q(0) + p
    F(3,1) = q(2)*H
end subroutine eulerFluxVec


subroutine eulerFlux(q,F,vec,gamma,smax)
    ! -----------------------------------------------------------------------
    ! Purpose: Calculates the analytical Euler flux given a state and normal
    ! 
    ! Inputs:
    !   q[4] = the state
    !   nrm[2] = the normal vector
    !   gamma = heat ratio
    ! 
    ! Outs:
    !   F[4] = analytical Euler flux
    ! 
    ! -----------------------------------------------------------------------
    implicit none
    real(8), intent(in), dimension(0:3) :: q
    real(8), intent(in), dimension(0:1) :: vec
    real(8), intent(in) :: gamma
    real(8), intent(out), dimension(0:3) :: F
    real(8), intent(out) :: smax
!f2py intent(in) q, nrm,gamma
!f2py intent(out) F
    real(8) :: velnorm,p,H,unorm,c
    real(8), dimension(0:3) :: F1, F2
    real(8), dimension(0:2) :: lambda
    real(8), dimension(0:1) :: nrm
    nrm = vec / norm2(vec) ! we always normalize nrm
    unorm = q(1)/q(0)*nrm(0) + q(2)/q(0)*nrm(1)
    velnorm = sqrt((q(1)/q(0))**2.d0 + (q(2)/q(0))**2.d0) ! this is the norm of velocity vector
    p = (gamma-1)*(q(3) - 0.5d0*q(0)*velnorm**2.d0)
    H = q(3)/q(0) + p/q(0)
    c = sqrt(gamma*p/q(0))
    lambda(0) = unorm+c
    lambda(1) = unorm-c
    lambda(2) = unorm
    smax = maxval(abs(lambda))
    ! ---------- F1 and F2 are the fluxes in the x and y directions
    F1(0) = q(1)
    F1(1) = q(1)**2.d0/q(0) + p
    F1(2) = q(1)*q(2)/q(0)
    F1(3) = q(1)*H

    F2(0) = q(2)
    F2(1) = q(1)*q(2)/q(0)
    F2(2) = q(2)**2.d0/q(0) + p
    F2(3) = q(2)*H

    F = F1(:)*nrm(0) + F2(:)*nrm(1) ! Here we project it to the normal direction
    if (norm2(vec) <= 5.d-16) then
        F = 0.d0
    endif
end subroutine eulerFlux


subroutine roeAvgState(qL,qR,qA,gamma)
    ! -----------------------------------------------------------------------
    ! Purpose: Calculates the Roe-averaged state qA from left and right states
    ! 
    ! Inputs:
    !   qL[4] = the left state
    !   qR[4] = the right state
    !   gamma = heat ratio
    ! 
    ! Outs:
    !   qA[3] = Roe-averaged states [u,v,H]
    ! 
    ! -----------------------------------------------------------------------
    implicit none
    real(8), intent(in), dimension(0:3) :: qL, qR
    real(8), intent(in) :: gamma
    real(8), intent(out), dimension(0:2) :: qA
!f2py intent(in) qL, qR, gamma
!f2py intent(out) qA
    real(8) :: velnormL, velnormR, pL, pR, HL, HR, RL, RR
    ! _____________________
    ! begin main execution 
    velnormL = sqrt((qL(1)/qL(0))**2.d0 + (qL(2)/qL(0))**2.d0)
    velnormR = sqrt((qR(1)/qR(0))**2.d0 + (qR(2)/qR(0))**2.d0)
    pL = (gamma-1.d0)*(qL(3) - 0.5d0*qL(0)*velnormL**2.d0)
    pR = (gamma-1.d0)*(qR(3) - 0.5d0*qR(0)*velnormR**2.d0)
    HL = qL(3)/qL(0) + pL/qL(0)
    HR = qR(3)/qR(0) + pR/qR(0)

    rL = sqrt(qL(0)) 
    rR = sqrt(qR(0))
    qA(0) = (rL*qL(1)/qL(0) + rR*qR(1)/qR(0))/(rL + rR)
    qA(1) = (rL*qL(2)/qL(0) + rR*qR(2)/qR(0))/(rL + rR)
    qA(2) = (rL*HL + rR*HR)/(rL + rR)
end subroutine roeAvgState


subroutine roeFlux(qL,qR,F,vec,gamma,smax)
    ! -----------------------------------------------------------------------
    ! Purpose: computes the Roe flux given left/right states and the normal
    ! 
    ! Inputs:
    !   qL[4] = the left state
    !   qL[4] = the left state
    !   nrm[2] = the normal vector, pointing from L to R
    !   gamma = heat ratio
    ! 
    ! Outs:
    !   F[4] = Roe flux vector
    ! 
    ! -----------------------------------------------------------------------
    implicit none
    real(8), intent(in), dimension(0:3) :: qL, qR
    real(8), intent(in), dimension(0:1) :: vec
    real(8), intent(in) :: gamma
    real(8), intent(out), dimension(0:3) :: F
    real(8), intent(out) :: smax
!f2py intent(in) qL, qR, nrm, gamma
!f2py intent(out) F,smax
    real(8), dimension(0:3) :: FL, FR, dq, R
    real(8), dimension(0:2) :: qA, lambda
    real(8) :: unorm, velnorm, H, c, M, G1, G2, s1, s2, C1, C2, epsilon
    real(8), dimension(0:1) :: nrm
    ! _____________________
    ! begin main execution 
    nrm = vec/norm2(vec) ! always normalize
    call eulerFlux(qL,FL,nrm,gamma,smax) ! these smax values are overwritten later on
    call eulerFlux(qR,FR,nrm,gamma,smax)

    dq = qR - qL ! difference in two states
    call roeAvgState(qL,qR,qA,gamma) ! remember that qA is not scaled by density! [rho,u,v,E]
    unorm = qA(0)*nrm(0) + qA(1)*nrm(1) ! this is velocity projected to the nrm dir
    velnorm = sqrt(qA(0)**2.d0 + qA(1)**2.d0) ! this is the norm of the velocity vector
    H = qA(2)
    c = sqrt((gamma-1.d0)*(H - 0.5d0*velnorm**2.d0))
    M = velnorm/c
    ! ------------------ see Kfid 623 notes for definitions/approach
    lambda(0) = unorm+c ! these are the eigenvalues of flux Jacobian
    lambda(1) = unorm-c
    lambda(2) = unorm
    epsilon = 0.1*c
    if (abs(lambda(0)) < epsilon) then
        lambda(0) = (epsilon**2 + lambda(0)**2)/(2*epsilon)
    endif
    if (abs(lambda(1)) < epsilon) then
        lambda(1) = (epsilon**2 + lambda(1)**2)/(2*epsilon)
    endif
    if (abs(lambda(2)) < epsilon) then
        lambda(2) = (epsilon**2 + lambda(2)**2)/(2*epsilon)
    endif
    smax = maxval(abs(lambda))
    G1 = (gamma-1)*(velnorm**2.d0/2.d0*dq(0) - qA(0)*dq(1) - qA(1)*dq(2) + dq(3))
    G2 = -unorm*dq(0) + dq(1)*nrm(0) + dq(2)*nrm(1)
    s1 = 0.5d0*(abs(lambda(0))+abs(lambda(1)))
    s2 = 0.5d0*(abs(lambda(0))-abs(lambda(1)))
    C1 = G1/c**2.d0*(s1-abs(lambda(2))) + G2/c*s2
    C2 = G1/c*s2 + (s1-abs(lambda(2)))*G2

    ! ----------------- R vector is the second vector subtracted from avg flux
    R(0) = abs(lambda(2))*dq(0) + C1
    R(1) = abs(lambda(2))*dq(1) + C1*qA(0) + C2*nrm(0)
    R(2) = abs(lambda(2))*dq(2) + C1*qA(1) + C2*nrm(1)
    R(3) = abs(lambda(2))*dq(3) + C1*H + C2*unorm

    F = 0.5d0*(FL + FR) - 0.5d0*R
end subroutine roeFlux