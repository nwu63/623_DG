! TODO: have precomputation/evaluation of quadratures in DG
subroutine getResidual(q,p,I2E,B2E,In,Bn,qnrm,rBC,resids,Jinv,Jinv2,detJ,detJ2,xy,w2,gphi,&
    w1,phiL,phiR,xyL,xyR,qlist,wavespeed,gamma,Rgas,nelem,niface,nbface,nqelem,Ng,Ng1)
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
    integer, intent(in) :: nelem,niface,nbface,p,nqelem,Ng,Ng1
    real(8), intent(in), dimension(nelem,(p+1)*(p+2)/2,4) :: q
    integer, intent(in), dimension(niface,4) :: I2E
    integer, intent(in), dimension(nbface,3) :: B2E
    real(8), intent(in), dimension(niface,3) :: In ! includes length as 3rd component
    real(8), intent(in), dimension(nbface,3) :: Bn ! includes length as 3rd component
    real(8), intent(in), dimension(nqelem,Ng1,2) :: qnrm ! nrm for high-q elems
    real(8), intent(in), dimension(nelem,2,2) :: Jinv
    real(8), intent(in), dimension(nqelem,Ng,2,2) :: Jinv2
    real(8), intent(in), dimension(nelem) :: detJ
    real(8), intent(in), dimension(nqelem,Ng) :: detJ2
    real(8), intent(in), dimension(Ng,(p+1)*(p+2)/2,2) :: gphi
    real(8), intent(in), dimension(Ng,2) :: xy
    real(8), intent(in), dimension(Ng) :: w2
    real(8), intent(in), dimension(Ng1) :: w1
    real(8), intent(in), dimension(3,Ng1,(p+1)*(p+2)/2) :: phiL,phiR
    real(8), intent(in), dimension(3,Ng1,2) :: xyL,xyR
    integer, intent(in), dimension(nqelem) :: qlist
    real(8), intent(in), dimension(5) :: rBC
    real(8), intent(in) :: gamma, Rgas
    real(8), intent(out), dimension(nelem,(p+1)*(p+2)/2,4) :: resids
    real(8), intent(out), dimension(nelem) :: wavespeed ! this is in fact the sum of length*wavespeed for each cell
!f2py intent(in) q,I2E,B2E,In,Bn,rBC,gamma, Rgas, Jinv, detJ,xy,w2,gphi,xyL,xyR,phiL,phiR
!f2py intent(out) resids,wavespeed
    integer :: iface, elemL,elemR,faceL,faceR,elem,btype,face,Nb,ig,ib,g1,g,idx
    real(8), dimension(2) :: nrm
    real(8), dimension(4) :: F,qBC,qI
    real(8) :: length,pinf,Mb,rhoinf,Tt,pt,alpha,Tb,pb,cb,Splus,Jplus,uI,vI,unb,unplus,cplus,pplus,ub,vb,dn,a,b,c,det,smax
    logical :: dirichlet ! sets all BC to Dirichlet equal to qBC
    logical :: qelem     ! if the current element is high-q or low-q
    real(8), dimension(:,:), allocatable :: vec, qState, qL, qR
    
    dirichlet = .false.
    qelem = .false.
    resids(:,:,:) = 0.d0 ! reset resids to zero
    wavespeed(:) = 0.d0
    ! -----------------------------------
    ! Pre-compute quantities
    ! -----------------------------------
    Nb = (p+1)*(p+2)/2
    allocate(vec(Nb,2))
    allocate(qState(Ng,4))
    allocate(qL(Ng1,4))
    allocate(qR(Ng1,4))

    ! ------------------- interior element contribution
    do elem = 1,nelem
        call getQ(q(elem,:,:),p,xy,Ng,qState)
        if (any(qlist == elem)) then ! curved element
            qelem = .true.
            idx = minloc((qlist-elem)**2,dim=1)
        else
            qelem = .false.
        endif
        do ig = 1,Ng
            if (qelem) then ! curved element
                vec = matmul(gphi(ig,:,:),Jinv2(idx,ig,:,:))             
            else
                vec = matmul(gphi(ig,:,:),Jinv(elem,:,:))
            endif
            do ib = 1,Nb
                call eulerFlux(qState(ig,:),F,vec(ib,:),gamma,smax) ! TODO: we actually only need to get F as 4x2 once
                if (qelem) then ! curved element
                    resids(elem,ib,:) = resids(elem,ib,:) - norm2(vec(ib,:))*F(:)*detJ2(idx,ig)*w2(ig)
                else
                    resids(elem,ib,:) = resids(elem,ib,:) - norm2(vec(ib,:))*F(:)*detJ(elem)*w2(ig)
                endif
            enddo
        enddo
    enddo

    ! ------------------- interior edge flux contribution
    do iface = 1, niface
        nrm = In(iface,1:2)
        length = In(iface,3)
        elemL = I2E(iface,1)
        faceL = I2E(iface,2)
        elemR = I2E(iface,3)
        faceR = I2E(iface,4)
        call getQ(q(elemL,:,:),p,xyL(faceL,:,:),Ng1,qL)
        call getQ(q(elemR,:,:),p,xyR(faceR,:,:),Ng1,qR) 
        do ig = 1,Ng1
            call roeFlux(qL(ig,:),qR(ig,:),F,nrm,gamma,smax)
            do ib = 1,Nb
                resids(elemL,ib,:) = resids(elemL,ib,:) + phiL(faceL,ig,ib)*F(:)*length*w1(ig)
                resids(elemR,ib,:) = resids(elemR,ib,:) - phiR(faceR,ig,ib)*F(:)*length*w1(Ng1-ig+1) ! apply weights in reverse
            enddo
            wavespeed(elemL) = wavespeed(elemL) + smax*length * w1(ig)
            wavespeed(elemR) = wavespeed(elemR) + smax*length * w1(Ng1-ig+1)
        enddo
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
        if (any(qlist == elem) .and. btype == 4) then ! curved element
            qelem = .true.
            idx = minloc((qlist-elem)**2,dim=1)
        else
            qelem = .false.
        endif

        do ig = 1,Ng1
            if (qelem) then
                length = norm2(qnrm(idx,ig,:))
                nrm = qnrm(idx,ig,:)/length
            endif
            call getQ(q(elem,:,:),p,xyL(face,ig,:),1,qI) ! we evaluate these one at a time to simplify code
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
                pplus = (gamma-1.d0)*(qI(4) - 0.5d0*qI(1)*(uI**2+vI**2))
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
        
            do ib = 1,Nb
                resids(elem,ib,:) = resids(elem,ib,:) + phiL(face,ig,ib)*F(:)*length*w1(ig)
            enddo ! ib
            wavespeed(elem) = wavespeed(elem) + smax*length*w1(ig)
        enddo ! ig
    enddo ! iface
end subroutine getResidual

subroutine getQ(q,p,xy,n_xy,qState)
    ! -----------------------------------------------------------------------
    ! Purpose: Evaluates the state qState, given by the coefficients
    !          q on a triangular element with Lagrange basis functions
    ! 
    ! Inputs:
    !   q[Nb,4] = basis coefficients for each state rank
    !   p = order of the polynomial. Naturally Nb = (p+1)*(p+2)/2
    !   xy [n_xy,2] = xi and eta locations for each point to be evaluated
    ! 
    ! Outs:
    !   qState[n_xy,4] = the state for each point
    ! 
    ! -----------------------------------------------------------------------
    implicit none
    integer, intent(in) :: p,n_xy
    real(8), intent(in), dimension(n_xy,2) :: xy
    real(8), intent(in), dimension((p+1)*(p+2)/2,4) :: q
    real(8), intent(out), dimension(n_xy,4) :: qState
!f2py intent(in) p,q
!f2py intent(out) qState
    real(8), dimension(n_xy,(p+1)*(p+2)/2) :: phi

    call basis2D(xy, p, phi, n_xy)
    qState = matmul(phi,q)
end subroutine getQ



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
    Minf = sqrt((Tt - 1.d0)*2.d0/(gamma-1.d0))
    u = Minf*c*cos(alpha)
    v = Minf*c*sin(alpha)
    q(1) = rhoinf
    q(2) = rhoinf*u
    q(3) = rhoinf*v
    q(4) = pinf/(gamma-1.d0) + 0.5d0*rhoinf*(u**2.d0+v**2.d0)
end subroutine getIC