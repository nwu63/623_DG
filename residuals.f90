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
    real(8), intent(in), dimension(nelem,4) :: q
    integer, intent(in), dimension(niface,4) :: I2E
    integer, intent(in), dimension(nbface,3) :: B2E
    real(8), intent(in), dimension(niface,3) :: In ! includes length as 3rd component
    real(8), intent(in), dimension(nbface,3) :: Bn ! includes length as 3rd component
    real(8), intent(in), dimension(5) :: rBC
    real(8), intent(in) :: gamma, Rgas
    real(8), intent(out), dimension(nelem,4) :: resids
    real(8), intent(out), dimension(nelem) :: wavespeed ! this is in fact the sum of length*wavespeed for each cell
!f2py intent(in) q,I2E,B2E,In,Bn,rBC,gamma, Rgas,area,centroidvec
!f2py intent(out) resids,wavespeed
    integer :: iface, elemL,elemR,faceL,faceR,elem,btype,face
    real(8), dimension(2) :: nrm
    real(8), dimension(4) :: F,qL,qR,qI,qBC
    real(8) :: length,pinf,Mb,rhoinf,Tt,pt,alpha,Tb,pb,cb,Splus,Jplus,uI,vI,unb,unplus,cplus,pplus,ub,vb,dn,a,b,c,det,smax
    logical :: dirichlet ! sets all BC to Dirichlet equal to qBC
    dirichlet = .true.
    resids(:,:) = 0 ! reset resids to zero
    wavespeed(:) = 0

    ! ------------------- interior flux
    do iface = 1, niface
        nrm = In(iface,1:2)
        length = In(iface,3)
        elemL = I2E(iface,1)
        faceL = I2E(iface,2)
        elemR = I2E(iface,3)
        faceR = I2E(iface,4)
        qL = q(elemL,:)
        qR = q(elemR,:)
        call roeFlux(qL,qR,F,nrm,gamma,smax)
        resids(elemL,:) = resids(elemL,:) + F(:)*length
        resids(elemR,:) = resids(elemR,:) - F(:)*length
        wavespeed(elemL) = wavespeed(elemL) + smax*length
        wavespeed(elemR) = wavespeed(elemR) + smax*length
    enddo
    ! -------------------- unpack rBC
    pinf = rBC(1)
    rhoinf = rBC(2)
    Tt = rBC(3)
    pt = rBC(4)
    alpha = rBC(5)
    ! -------------------- boundary faces
    do iface = 1, nbface
        nrm = Bn(iface,1:2)
        length = Bn(iface,3)
        elem = B2E(iface,1)
        btype = B2E(iface,3)
        face = B2E(iface,2)
        qI = q(elem,:)
        uI = qI(2)/qI(1) ! internal velocity components
        vI = qI(3)/qI(1)
        unplus = uI*nrm(1) + vI*nrm(2) ! internal projected velocity
        if (dirichlet) then ! DIRICHLET BC
            call getIC(rBC,qBC,gamma)
            call roeFlux(qI,qBC,F,nrm,gamma,smax)
        else if (btype == 3 .or. btype == 4) then ! WALL BC
            ub = uI - unplus*nrm(1)
            vb = vI - unplus*nrm(2)
            pb = (gamma-1.d0)*(qI(4) - 0.5d0*qI(1)*(ub**2+vb**2))
            ! ------- directly construct flux
            F(1) = 0.d0
            F(2) = pb*nrm(1)
            F(3) = pb*nrm(2)
            F(4) = 0.d0
            ! should you zero out smax here?
            ! smax = 0.d0
        else if (btype == 1) then ! INLET BC
            dn = cos(alpha)*nrm(1) + sin(alpha)*nrm(2)
            pplus = (gamma-1.d0)*(qI(4) - 0.5d0*qI(1)*(uI**2+vI**2))
            cplus = sqrt(gamma*pplus/qI(1))
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
            qBC(1) = pb/(Rgas*Tb)
            cb = sqrt(gamma*pb/qBC(1))
            qBC(2) = Mb*cb*cos(alpha)
            qBC(3) = Mb*cb*sin(alpha)
            qBC(4) = pb/(gamma-1.d0)+0.5d0*qBC(1)*(qBC(2)**2+qBC(3)**2)
            qBC(2) = qBC(2)*qBC(1)
            qBC(3) = qBC(3)*qBC(1)
            call eulerFlux(qBC,F,nrm,gamma,smax)
        else if (btype == 2) then ! OUTLET BC
            pplus = (gamma-1.d0)*(qI(3) - 0.5d0*qI(1)*(uI**2+vI**2))
            Splus = pplus/(qI(1)**gamma)
            qBC(1) = (pinf/Splus)**(1.d0/gamma)
            cplus = sqrt(gamma*pplus/qI(1))
            Jplus = unplus + 2*cplus/(gamma-1.d0)
            cb = sqrt(gamma*pinf/qBC(1))
            unb = Jplus - 2*cb/(gamma-1.d0)
            qBC(2) = uI + (unb-unplus)*nrm(1)
            qBC(3) = vI + (unb-unplus)*nrm(2)
            qBC(4) = pinf/(gamma-1.d0) + 0.5d0 * qBC(1)*(qBC(2)**2 + qBC(3)**2)
            qBC(2) = qBC(2) * qBC(1)
            qBC(3) = qBC(3) * qBC(1)
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
    real(8), intent(out), dimension(4) :: q
    real(8), intent(in), dimension(5) :: rBC
    real(8), intent(in) :: gamma
!f2py intent(in) gamma
!f2py intent(out) q
    real(8) :: pinf,rhoinf,Minf,Tt,pt,alpha,c,u,v
    pinf = rBC(1)
    rhoinf = rBC(2)
    Tt = rBC(3)
    pt = rBC(4)
    alpha = rBC(5)
    c = sqrt(gamma*pinf/rhoinf)
    Minf = sqrt((Tt - 1)*2/(gamma-1))
    u = Minf*c*cos(alpha)
    v = Minf*c*sin(alpha)
    q(1) = rhoinf
    q(2) = rhoinf*u
    q(3) = rhoinf*v
    q(4) = pinf/(gamma-1) + 0.5*rhoinf*(u**2+v**2)
end subroutine getIC