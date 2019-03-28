subroutine integrate(q,p,B2E,Bn,rBC,detJ,cl,cd,Es,cp,mach,gamma,Rgas,nelem,nbface)
    ! -----------------------------------------------------------------------
    ! Purpose: Calculates the integrated values: cl, cd, Es.
    !          Also calculates mach number for each elem, and
    !          cp for the bottom face
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
    !   cp[nelem,2] we don't know how many elements on the bottom, so we initialize
    !               for the whole vector. First entry is elem number, second the cp
    !   mach[nelem]
    ! 
    ! -----------------------------------------------------------------------
    implicit none
    integer, intent(in) :: nelem,nbface, p
    real(8), intent(in), dimension(nelem,(p+1)*(p+2)/2,4) :: q
    integer, intent(in), dimension(nbface,3) :: B2E
    real(8), intent(in), dimension(nbface,3) :: Bn ! includes length as 3rd component
    real(8), intent(in), dimension(5) :: rBC
    real(8), intent(in), dimension(nelem) :: detJ
    real(8), intent(in) :: gamma, Rgas
    real(8), intent(out) :: cl,cd,Es
    real(8), intent(out), dimension(nelem,2) :: cp
    real(8), intent(out), dimension(nelem) :: mach
!f2py intent(in) q, p, B2E,Bn,rBC,detJ
!f2py intent(out) cl,cd,Es,mach,cp
    integer :: iface,btype,elem,idx,face,Nb,Ng,Ng1,g1,ig
    real(8), dimension(4) :: qState
    real(8), dimension(2) :: nrm
    real(8) :: uI,vI,ub,vb,length,pb,pinf,rhoinf,Tt,pt,Minf,h,rhot,st,s,c,pdyn,denom
    real(8), dimension(:,:,:), allocatable :: xy1, phi1 ! xy1 is the edge integration points on T
    real(8), dimension(:,:), allocatable :: xy, qB, qI
    real(8), dimension(:), allocatable :: w1,w2,x

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
    cp(:,:) = -1

    ! -------------------------------------------
    ! Pre-computing Integration Quantities
    ! -------------------------------------------
    Nb = (p+1)*(p+2)/2
    call Gauss2D_pre(2*p+1,Ng)
    allocate(xy(Ng,2))
    allocate(w2(Ng))
    allocate(qI(Ng,4))
    call Gauss2D(2*p+1,Ng,xy,w2)
    
    g1 = 2*p+1
    Ng1 = (g1+2)/2
    allocate(x(Ng1))
    allocate(w1(Ng1))
    allocate(xy1(3,Ng1,2))
    allocate(phi1(3,Ng1,Nb))
    call Gauss1D(g1,x,w1)
    allocate(qB(Ng1,4))
    ! now we map these 1D quadrature points to T, one for each edge in ccw order
    xy1(1,:,1) = 1-x
    xy1(1,:,2) = x
    xy1(2,:,1) = 0
    xy1(2,:,2) = 1-x
    xy1(3,:,1) = x
    xy1(3,:,2) = 0
    call basis2D(xy1(1,:,:), p, phi1(1,:,:), Ng1)
    call basis2D(xy1(2,:,:), p, phi1(2,:,:), Ng1)
    call basis2D(xy1(3,:,:), p, phi1(3,:,:), Ng1)



    do iface = 1, nbface
        btype = B2E(iface,3)
        if (btype == 4) then
            nrm = Bn(iface,1:2)
            length = Bn(iface,3)
            elem = B2E(iface,1)
            face = B2E(iface,2)
            call getQ(q(elem,:,:),p,xy1(face,:,:),Ng1,qB)
            do ig = 1,Ng1
                qState = qB(ig,:)
                ub = qState(2)/qState(1)
                vb = qState(3)/qState(1)
                pb = (gamma-1.d0)*(qState(4) - 0.5d0*qState(1)*(ub**2+vb**2))
                cl = cl + (pb-pinf)*nrm(2)*length*w1(ig)
                cd = cd + (pb-pinf)*nrm(1)*length*w1(ig)
            enddo
            cp(idx,1) = elem
            cp(idx,2) = pb-pinf
            idx = idx + 1
        endif
    enddo
    cl = cl/pdyn/h
    cd = cd/pdyn/h
    cp(:,2) = cp(:,2)/pdyn

    rhot = pt/(Rgas*Tt)
    st = pt/rhot**gamma
    Es = 0
    denom = 0
    do elem = 1, nelem
        call getQ(q(elem,:,:),p,xy,Ng,qI)
        do ig = 1,Ng
            qState = qI(ig,:)
            uI = qState(2)/qState(1)
            vI = qState(3)/qState(1)
            pb = (gamma-1.d0)*(qState(4) - 0.5d0*qState(1)*(uI**2+vI**2)) ! this is actually just p, but we reuse the variable
            s = pb/qState(1)**gamma
            Es = Es + (s/st - 1)**2 * detJ(elem)*w2(ig)
            denom = denom + detJ(elem)*w2(ig)
        enddo
        ! Es(elem) = s/st - 1
        c = sqrt(gamma*pb/qState(1))
        mach(elem) = sqrt(uI**2 + vI**2)/c
    enddo
    Es = sqrt(Es / denom)

    ! print*, cl,cd,Es,mach,cp
    
    ! deallocate(xy1, phi1,xy, qB, qI,w1,w2,x)
end subroutine