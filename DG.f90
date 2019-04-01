!TODO: refactor code to pre-compute 1D and 2D integration points here.
!      Also pre-allocate things like Jacobians etc for different q elements
subroutine DG(q,p,geom,resids,maxres,detJ,nodes,qlist,E2N1,E2N2,I2E,B2E,In,Bn,rBC,gamma,Rgas,CFL,convtol,min_iter,&
    max_iter,nnodes,nelem,niface,nbface,nqelem)
    implicit none
    integer, intent(in) :: p,nnodes,nelem,niface,nbface,min_iter,max_iter,nqelem,geom
    real(8), intent(inout), dimension(nelem,(p+1)*(p+2)/2,4) :: q
    real(8), intent(in), dimension(nnodes,2) :: nodes
    integer, intent(in), dimension(nqelem) :: qlist
    integer, intent(in), dimension(nelem,3) :: E2N1
    integer, intent(in), dimension(nqelem,(geom+1)*(geom+2)/2) :: E2N2
    integer, intent(in), dimension(niface,4) :: I2E
    integer, intent(in), dimension(nbface,3) :: B2E
    real(8), intent(in), dimension(niface,3) :: In ! includes length as 3rd component
    real(8), intent(in), dimension(nbface,3) :: Bn ! includes length as 3rd component
    real(8), intent(out),dimension(nelem,(p+1)*(p+2)/2,4) :: resids
    real(8), intent(in), dimension(5) :: rBC
    real(8), intent(in) :: gamma, Rgas, convtol, CFL
    real(8), intent(out), dimension(max_iter) :: maxres
    real(8), dimension(nelem),intent(out) :: detJ
!f2py intent(in) node,E2N1,E2N2,I2E,B2E,In,Bn,gamma,rBC,qlist
!f2py intent(out) resids, maxres, detJ
!f2py intent(in,out) q
    integer :: ielem, Nb, Ng, g, idx, g1, Ng1,iface,face,ig
    real(8) :: temp2
    real(8), dimension(nelem,2,2) :: J, Jinv
    real(8), dimension(:,:,:,:), allocatable :: J2, Jinv2
    real(8), dimension(:,:,:), allocatable :: Minv, gphi, xyL, phiL,xyR, phiR,qnrm ! xy1 is the edge integration points on T
    real(8), dimension(:,:), allocatable :: xy,phi, detJ2, qlength
    real(8), dimension(:), allocatable :: w,x,w1
    real(8), dimension(nelem) :: wavespeed
    real(8), dimension(2,2) :: temp1,Jedge
    real(8), dimension(2) :: vec,tangent
    ! -----------------------------------
    ! Allocate and precompute basis values at integration points
    ! -----------------------------------
    Nb = (p+1)*(p+2)/2
    allocate(Minv(nelem,Nb,Nb))

    ! -----------------------------------
    ! 2D integration
    ! -----------------------------------

    g = 2*p+1 + 2*(geom-1) ! 2D integration order for everything
    call Gauss2D_pre(g,Ng)
    allocate(xy(Ng,2))
    allocate(w(Ng))
    allocate(phi(Ng,Nb))
    allocate(gphi(Ng,Nb,2))
    call Gauss2D(g,Ng,xy,w)
    call basis2D(xy, p, phi, Ng)
    call gbasis2D(xy, p, gphi, Ng)
    allocate(J2(nqelem,Ng,2,2))
    allocate(Jinv2(nqelem,Ng,2,2))
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
    allocate(xyR(3,Ng1,2))
    allocate(phiR(3,Ng1,Nb))
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
    xyR(:,:,:) = xyL(:,Ng1:1:-1,:) ! Now we flip them for the right element
    call basis2D(xyR(1,:,:), p, phiR(1,:,:), Ng1)
    call basis2D(xyR(2,:,:), p, phiR(2,:,:), Ng1)
    call basis2D(xyR(3,:,:), p, phiR(3,:,:), Ng1)
    ! -----------------------------------
    ! Jacobians
    ! -----------------------------------
    do ielem = 1,nelem
        call getJacobian(nodes(E2N1(ielem,:),:), J(ielem,:,:), Jinv(ielem,:,:), detJ(ielem))
        if (any(qlist == ielem)) then ! curved element
            idx = minloc((qlist-ielem)**2,dim=1)
            call getHOJacobian(nodes(E2N2(idx,:),:), geom, xy, J2(idx,:,:,:), Jinv2(idx,:,:,:), detJ2(idx,:),Ng)
        endif
    enddo
    ! -----------------------------------
    ! Normals and Edge Jacobians
    ! -----------------------------------
    allocate(qnrm(nqelem,Ng1,2))
    allocate(qlength(nqelem,Ng1))
    do iface = 1, nbface
        ielem = B2E(iface,1)
        face = B2E(iface,2)
        if (any(qlist == ielem)) then ! curved element
            idx = minloc((qlist-ielem)**2,dim=1)
            call getdXidSigma(face,vec)
            do ig = 1,Ng1
                call getHOJacobian(nodes(E2N2(idx,:),:), geom, xyL(face,ig,:), Jedge, temp1, temp2,1)
                tangent = Jedge(:,1)*vec(1) + Jedge(:,2)*vec(2) ! TODO vectorize
                qnrm(idx,ig,1) = tangent(2)
                qnrm(idx,ig,2) = -tangent(1)
            enddo
        endif
    enddo

    ! -----------------------------------
    ! Mass matrix
    ! -----------------------------------
    call getMassInv(p,detJ,detJ2,qlist,Minv,nelem,nqelem,Ng)

    ! -----------------------------------
    ! Timestep
    ! -----------------------------------
    
    call timeIntegration(q,p,I2E,B2E,In,Bn,qnrm,Jinv,Jinv2,detJ,detJ2,Minv,xy,w,gphi,w1,rBC,resids,&
    phiL,phiR,xyL,xyR,qlist,maxres,gamma,Rgas,CFL,convtol,min_iter,max_iter,nelem,niface,nbface,nqelem,Ng,Ng1)

    ! call getResidual(q,p,I2E,B2E,In,Bn,qnrm,rBC,resids,Jinv,Jinv2,detJ,detJ2,xy,w,gphi,w1,&
    ! phiL,phiR,xyL,xyR,qlist,wavespeed,gamma,Rgas,nelem,niface,nbface,nqelem,Ng,Ng1)
end subroutine DG