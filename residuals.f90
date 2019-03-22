subroutine getResidual(q,p,I2E,B2E,In,Bn,rBC,resids,wavespeed,gamma,Rgas,nelem,niface,nbface)
    ! -----------------------------------------------------------------------
    ! Purpose: Calculates the residual over the entire domain
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
    ! 
    ! Outs:
    !   resids[nelem,4] = residual
    !   wavespeed[nelem] = maximum wavespeed, edge-weighted already
    ! 
    ! -----------------------------------------------------------------------
    implicit none
    integer, intent(in) :: nelem,niface,nbface,p
    real(8), intent(in), dimension(0:nelem-1,0:3) :: q
    integer, intent(in), dimension(0:niface-1,0:3) :: I2E
    integer, intent(in), dimension(0:nbface-1,0:2) :: B2E
    real(8), intent(in), dimension(0:niface-1,0:2) :: In ! includes length as 3rd component
    real(8), intent(in), dimension(0:nbface-1,0:2) :: Bn ! includes length as 3rd component
    real(8), intent(in), dimension(0:4) :: rBC
    real(8), intent(in) :: gamma, Rgas
    real(8), intent(out), dimension(0:nelem-1,0:3) :: resids
    real(8), intent(out), dimension(0:nelem-1) :: wavespeed ! this is in fact the sum of length*wavespeed for each cell
!f2py intent(in) q,I2E,B2E,In,Bn,rBC,gamma, Rgas,area,centroidvec
!f2py intent(out) resids,wavespeed
    integer :: iface, elemL,elemR,faceL,faceR,elem,btype,face
    real(8), dimension(0:1) :: nrm
    real(8), dimension(0:3) :: F,qL,qR,qI,qBC
    real(8) :: length,pinf,Mb,rhoinf,Tt,pt,alpha,Tb,pb,cb,Splus,Jplus,uI,vI,unb,unplus,cplus,pplus,ub,vb,dn,a,b,c,det,smax
    logical :: dirichlet ! sets all BC to Dirichlet equal to qBC
    dirichlet = .true.
    resids(:,:) = 0 ! reset resids to zero
    wavespeed(:) = 0

    ! ------------------- interior flux
    do iface = 0, niface-1
        nrm = In(iface,0:1)
        length = In(iface,2)
        elemL = I2E(iface,0)
        faceL = I2E(iface,1)
        elemR = I2E(iface,2)
        faceR = I2E(iface,3)
        qL = q(elemL,:)
        qR = q(elemR,:)
        call roeFlux(qL,qR,F,nrm,gamma,smax)
        resids(elemL,:) = resids(elemL,:) + F(:)*length
        resids(elemR,:) = resids(elemR,:) - F(:)*length
        wavespeed(elemL) = wavespeed(elemL) + smax*length
        wavespeed(elemR) = wavespeed(elemR) + smax*length
    enddo
    ! -------------------- unpack rBC
    pinf = rBC(0)
    rhoinf = rBC(1)
    Tt = rBC(2)
    pt = rBC(3)
    alpha = rBC(4)
    ! -------------------- boundary faces
    do iface = 0, nbface-1
        nrm = Bn(iface,0:1)
        length = Bn(iface,2)
        elem = B2E(iface,0)
        btype = B2E(iface,2)
        face = B2E(iface,1)
        qI = q(elem,:)
        uI = qI(1)/qI(0) ! internal velocity components
        vI = qI(2)/qI(0)
        unplus = uI*nrm(0) + vI*nrm(1) ! internal projected velocity
        if (dirichlet) then ! DIRICHLET BC
            call getIC(rBC,qBC,gamma)
            call roeFlux(qI,qBC,F,nrm,gamma,smax)
        else if (btype == 3 .or. btype == 4) then ! WALL BC
            ub = uI - unplus*nrm(0)
            vb = vI - unplus*nrm(1)
            pb = (gamma-1.d0)*(qI(3) - 0.5d0*qI(0)*(ub**2+vb**2))
            ! ------- directly construct flux
            F(0) = 0.d0
            F(1) = pb*nrm(0)
            F(2) = pb*nrm(1)
            F(3) = 0.d0
            ! should you zero out smax here?
            ! smax = 0.d0
        else if (btype == 1) then ! INLET BC
            dn = cos(alpha)*nrm(0) + sin(alpha)*nrm(1)
            pplus = (gamma-1.d0)*(qI(3) - 0.5d0*qI(0)*(uI**2+vI**2))
            cplus = sqrt(gamma*pplus/qI(0))
            Jplus = unplus + 2*cplus/(gamma-1)
            ! --------------- solve for quadratic in Mb
            a = gamma*Rgas*Tt*dn**2 - (gamma-1)/2*Jplus**2
            b = 4*gamma*Rgas*Tt*dn/(gamma-1)
            c = 4*gamma*Rgas*Tt/(gamma-1)**2 - Jplus**2
            det = sqrt(b**2-4*a*c)
            if (-b/(2*a)-det/(2*a) < 0) then
                Mb = (-b+det)/(2*a)
            else
                Mb = (-b-det)/(2*a)
            endif
            Tb = Tt/(1d0 + 0.5d0*(gamma-1d0)*Mb**2)
            pb = pt*(Tb/Tt)**(gamma/(gamma-1d0))
            qBC(0) = pb/(Rgas*Tb)
            cb = sqrt(gamma*pb/qBC(0))
            qBC(1) = Mb*cb*cos(alpha)
            qBC(2) = Mb*cb*sin(alpha)
            qBC(3) = pb/(gamma-1.d0)+0.5d0*qBC(0)*(qBC(1)**2+qBC(2)**2)
            qBC(1) = qBC(1)*qBC(0)
            qBC(2) = qBC(2)*qBC(0)
            call eulerFlux(qBC,F,nrm,gamma,smax)
        else if (btype == 2) then ! OUTLET BC
            pplus = (gamma-1.d0)*(qI(3) - 0.5d0*qI(0)*(uI**2+vI**2))
            Splus = pplus/(qI(0)**gamma)
            qBC(0) = (pinf/Splus)**(1.d0/gamma)
            cplus = sqrt(gamma*pplus/qI(0))
            Jplus = unplus + 2*cplus/(gamma-1.d0)
            cb = sqrt(gamma*pinf/qBC(0))
            unb = Jplus - 2*cb/(gamma-1.d0)
            qBC(1) = uI + (unb-unplus)*nrm(0)
            qBC(2) = vI + (unb-unplus)*nrm(1)
            qBC(3) = pinf/(gamma-1.d0) + 0.5d0 * qBC(0)*(qBC(1)**2 + qBC(2)**2)
            qBC(1) = qBC(1) * qBC(0)
            qBC(2) = qBC(2) * qBC(0)
            call eulerFlux(qBC,F,nrm,gamma,smax)
        end if
        resids(elem,:) = resids(elem,:) + F(:)*length
        wavespeed(elem) = wavespeed(elem) + smax*length
    enddo
end subroutine getResidual

subroutine getIC(rBC,q,gamma)
    ! -----------------------------------------------------------------------
    ! Purpose: Calculates the initial condition, similar to the python routine
    !          currently only used for qBC in Dirichlet BC
    ! 
    ! Inputs:
    !   rBC[5] = information on freestream/boundary values
    !   gamma = heat ratio
    ! 
    ! Outs:
    !   q[4] = the state
    ! 
    ! -----------------------------------------------------------------------
    implicit none
    real(8), intent(out), dimension(0:3) :: q
    real(8), intent(in), dimension(0:4) :: rBC
    real(8), intent(in) :: gamma
!f2py intent(in) gamma
!f2py intent(out) q
    real(8) :: pinf,rhoinf,Minf,Tt,pt,alpha,c,u,v
    pinf = rBC(0)
    rhoinf = rBC(1)
    Tt = rBC(2)
    pt = rBC(3)
    alpha = rBC(4)
    c = sqrt(gamma*pinf/rhoinf)
    Minf = sqrt((Tt - 1)*2/(gamma-1))
    u = Minf*c*cos(alpha)
    v = Minf*c*sin(alpha)
    q(0) = rhoinf
    q(1) = rhoinf*u
    q(2) = rhoinf*v
    q(3) = pinf/(gamma-1) + 0.5*rhoinf*(u**2+v**2)
end subroutine getIC