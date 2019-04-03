! TODO: make this work for curved elements too
subroutine integrate(q,p,geom,nodes,qlist,E2N1,E2N2,B2E,Bn,rBC,cl,cd,Es,gamma,Rgas,nelem,nnodes,nqelem,nbface)
    ! -----------------------------------------------------------------------
    ! Purpose: Calculates the integrated values: cl, cd, Es.
    ! 
    ! Inputs:
    !   q[nelem,4] = the initial state
    !   grad[nelem,4,2] = the solution gradient
    !   B2E[nbface,3]
    !   Bn[nbface,3] = normal vector for boundary faces, and length
    !   rBC[5] = information on freestream/boundary values
    !   area[nelem] 
    !   centroidvec[nelem,3,2] = vector pointing from midpoint to centroid
    ! 
    ! Outs:
    !   cl, cd, Es
    ! 
    ! -----------------------------------------------------------------------
    implicit none
    integer, intent(in) :: nelem,nbface, p,nnodes,nqelem,geom
    real(8), intent(in), dimension(nelem,(p+1)*(p+2)/2,4) :: q
    real(8), intent(in), dimension(nnodes,2) :: nodes
    integer, intent(in), dimension(nqelem) :: qlist
    integer, intent(in), dimension(nelem,3) :: E2N1
    integer, intent(in), dimension(nqelem,(geom+1)*(geom+2)/2) :: E2N2
    integer, intent(in), dimension(nbface,3) :: B2E
    real(8), intent(in), dimension(nbface,3) :: Bn ! includes length as 3rd component
    real(8), intent(in), dimension(5) :: rBC
    real(8), intent(in) :: gamma, Rgas
    real(8), intent(out) :: cl,cd,Es
!f2py intent(in) q,p,B2E,Bn,rBC,nodes,qlist,E2N1,E2N2
!f2py intent(out) cl,cd,Es,cp
    integer :: iface,btype,elem,idx,face,Nb,ig,ielem,idx2,g,g1,Ng,Ng1
    real(8), dimension(4) :: qState
    real(8), dimension(2) :: nrm,vec,tangent
    real(8), dimension(2,2) :: Jedge
    real(8) :: uI,vI,ub,vb,length,pb,pinf,rhoinf,Tt,pt,Minf,h,rhot,st,s,pdyn
    real(8), dimension(:,:,:), allocatable :: xyL, phiL,J,Jinv,qnrm ! xyL is the edge integration points on T
    real(8), dimension(:,:), allocatable :: xy, qB, qI,detJ2,phi
    real(8), dimension(:), allocatable :: w1,w,x
    real(8), dimension(nelem) :: detJ
    logical :: qelem
    pinf = rBC(1)
    rhoinf = rBC(2)
    Tt = rBC(3)
    pt = rBC(4)
    Minf = sqrt((Tt - 1)*2/(gamma-1))
    h = 0.0625
    pdyn = gamma/2.d0*pinf*Minf**2
    cl = 0
    cd = 0
    idx = 1

    ! -----------------------------------
    ! 2D integration
    ! -----------------------------------
    Nb = (p+1)*(p+2)/2
    g = 2*p+1 + 2*(geom-1) ! 2D integration order for everything
    call Gauss2D_pre(g,Ng)
    allocate(xy(Ng,2))
    allocate(w(Ng))
    allocate(phi(Ng,Nb))
    call Gauss2D(g,Ng,xy,w)
    call basis2D(xy, p, phi, Ng)
    allocate(J(Ng,2,2))
    allocate(Jinv(Ng,2,2))
    allocate(detJ2(nqelem,Ng))

    ! -----------------------------------
    ! 1D integration
    ! -----------------------------------
    g1 = 2*p+1
    Ng1 = (g1+2)/2
    allocate(x(Ng1))
    allocate(w1(Ng1))
    allocate(xyL(3,Ng1,2))
    allocate(phiL(3,Ng1,Nb))
    call Gauss1D(g1,x,w1)
    ! now we map these 1D quadrature points to T, one for each edge in ccw order
    xyL(1,:,1) = 1-x
    xyL(1,:,2) = x
    xyL(2,:,1) = 0
    xyL(2,:,2) = 1-x
    xyL(3,:,1) = x
    xyL(3,:,2) = 0
    call basis2D(xyL(1,:,:), p, phiL(1,:,:), Ng1)
    call basis2D(xyL(2,:,:), p, phiL(2,:,:), Ng1)
    call basis2D(xyL(3,:,:), p, phiL(3,:,:), Ng1)
    ! -----------------------------------
    ! Jacobians
    ! -----------------------------------
    do ielem = 1,nelem
        call getJacobian(nodes(E2N1(ielem,:),:), J(1,:,:), Jinv(1,:,:), detJ(ielem))
        if (any(qlist == ielem)) then ! curved element
            idx = minloc((qlist-ielem)**2,dim=1)
            call getHOJacobian(nodes(E2N2(idx,:),:), geom, xy, J(:,:,:), Jinv(:,:,:), detJ2(idx,:),Ng)
            detJ(ielem) = dot_product(detJ2(idx,:),w(:))*2 ! ! We multiply by 2 so it's the same as detJ = 2*area
        endif
    enddo
    ! -----------------------------------
    ! Normals and Edge Jacobians
    ! -----------------------------------
    allocate(qnrm(nqelem,Ng1,2))
    do iface = 1, nbface
        ielem = B2E(iface,1)
        face = B2E(iface,2)
        if (any(qlist == ielem)) then ! curved element
            idx = minloc((qlist-ielem)**2,dim=1)
            call getdXidSigma(face,vec)
            do ig = 1,Ng1
                call getHOJacobian(nodes(E2N2(idx,:),:), geom, xyL(face,ig,:), Jedge, Jinv(1,:,:), Jinv(1,:,:),1)
                tangent = Jedge(:,1)*vec(1) + Jedge(:,2)*vec(2) ! TODO vectorize
                qnrm(idx,ig,1) = tangent(2)
                qnrm(idx,ig,2) = -tangent(1)
            enddo
        endif
    enddo

    ! -----------------------------------
    ! Boundary integration
    ! -----------------------------------
    allocate(qB(Ng1,4))
    idx2 = 1
    do iface = 1, nbface
        btype = B2E(iface,3)
        if (btype == 4) then
            nrm = Bn(iface,1:2)
            length = Bn(iface,3)
            elem = B2E(iface,1)
            face = B2E(iface,2)
            call getQ(q(elem,:,:),p,phiL(face,:,:),Ng1,qB)
            if (any(qlist == elem)) then ! curved element
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
                qState = qB(ig,:)
                ub = qState(2)/qState(1)
                vb = qState(3)/qState(1)
                pb = (gamma-1.d0)*(qState(4) - 0.5d0*qState(1)*(ub**2+vb**2))
                cl = cl + (pb-pinf)*nrm(2)*length*w1(ig)
                cd = cd + (pb-pinf)*nrm(1)*length*w1(ig)
            enddo
        endif
    enddo
    cl = cl/pdyn/h
    cd = cd/pdyn/h

    ! -----------------------------------
    ! Area integration
    ! -----------------------------------
    allocate(qI(Ng,4))
    rhot = pt/(Rgas*Tt)
    st = pt/rhot**gamma
    Es = 0.d0
    do elem = 1, nelem
        call getQ(q(elem,:,:),p,phi,Ng,qI)
        if (any(qlist == elem)) then ! curved element
            qelem = .true.
            idx = minloc((qlist-elem)**2,dim=1)
        else
            qelem = .false.
        endif
        do ig = 1,Ng
            qState = qI(ig,:)
            uI = qState(2)/qState(1)
            vI = qState(3)/qState(1)
            pb = (gamma-1.d0)*(qState(4) - 0.5d0*qState(1)*(uI**2+vI**2)) ! this is actually just p, but we reuse the variable
            s = pb/qState(1)**gamma
            if (.not. qelem) then
                Es = Es + (s/st - 1)**2 * detJ(elem)*w(ig)
            else
                Es = Es + (s/st - 1)**2 * detJ2(idx,ig)*w(ig)
            endif
        enddo
    enddo
    Es = sqrt(Es / (sum(detJ)/2))
end subroutine

subroutine getM(q,p,mach,xy,n_xy,gamma)
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
    real(8), intent(in) :: gamma
    integer, intent(in) :: p,n_xy
    real(8), intent(in), dimension(n_xy,2) :: xy
    real(8), intent(in), dimension((p+1)*(p+2)/2,4) :: q
    real(8), intent(out),dimension(n_xy) :: mach
!f2py intent(in) p,q,xy
!f2py intent(out) mach
    real(8), dimension(n_xy,(p+1)*(p+2)/2) :: phi
    integer :: ig
    real(8) :: uI,vI,pb,c
    real(8), dimension(n_xy,4) :: qelem
    real(8), dimension(4) :: qState

    call basis2D(xy, p, phi, n_xy)
    call getQ(q,p,phi,n_xy,qelem)
    do ig = 1,n_xy
        qState = qelem(ig,:)
        uI = qState(2)/qState(1)
        vI = qState(3)/qState(1)
        pb = (gamma-1.d0)*(qState(4) - 0.5d0*qState(1)*(uI**2+vI**2))
        c = sqrt(gamma*pb/qState(1))
        mach(ig) = sqrt(uI**2 + vI**2)/c
    enddo
end subroutine getM

subroutine getCp(q,p,cp,rBC,xy,n_xy,gamma)
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
    real(8), intent(in) :: gamma
    integer, intent(in) :: p,n_xy
    real(8), intent(in), dimension(5) :: rBC
    real(8), intent(in), dimension(n_xy,2) :: xy
    real(8), intent(in), dimension((p+1)*(p+2)/2,4) :: q
    real(8), intent(out),dimension(n_xy) :: cp
!f2py intent(in) p,q,xy
!f2py intent(out) mach
    real(8), dimension(n_xy,(p+1)*(p+2)/2) :: phi
    integer :: ig
    real(8) :: uI,vI,pb,pdyn,pinf,Tt,Minf
    real(8), dimension(n_xy,4) :: qelem
    real(8), dimension(4) :: qState
    
    pinf = rBC(1)
    Tt = rBC(3)
    Minf = sqrt((Tt - 1)*2/(gamma-1))
    pdyn = gamma/2.d0*pinf*Minf**2
    call basis2D(xy, p, phi, n_xy)
    call getQ(q,p,phi,n_xy,qelem)
    do ig = 1,n_xy
        qState = qelem(ig,:)
        uI = qState(2)/qState(1)
        vI = qState(3)/qState(1)
        pb = (gamma-1.d0)*(qState(4) - 0.5d0*qState(1)*(uI**2+vI**2))
        cp(ig) = pb-pinf
    enddo
    cp = cp/pdyn
end subroutine getCp