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
    integer :: ielem, Nb, Ng, g, idx
    real(8), dimension(nelem,2,2) :: J, Jinv
    real(8), dimension(:,:,:,:), allocatable :: J2, Jinv2
    real(8), dimension(:,:,:), allocatable :: Minv, gphi
    real(8), dimension(:,:), allocatable :: xy,phi, detJ2
    real(8), dimension(:), allocatable :: w
    
    ! -----------------------------------
    ! Define constants and allocate arrays
    ! -----------------------------------
    Nb = (p+1)*(p+2)/2
    allocate(Minv(nelem,Nb,Nb))

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
    ! Pre-compute quantities
    ! -----------------------------------
    
    do ielem = 1,nelem
        call getJacobian(nodes(E2N1(ielem,:),:), J(ielem,:,:), Jinv(ielem,:,:), detJ(ielem))
        print*, ielem
        if (any(qlist == ielem)) then ! curved element
            idx = minloc((qlist-ielem)**2,dim=1)
            call getHOJacobian(nodes(E2N2(idx,:),:), geom, xy, J2(idx,:,:,:), Jinv2(idx,:,:,:), detJ2(idx,:),Ng)
        endif
    enddo
    call getMassInv(p,detJ,detJ2,qlist,Minv,nelem,nqelem,Ng)
    
    ! -----------------------------------
    ! Timestep
    ! -----------------------------------
    
    call timeIntegration(q,p,I2E,B2E,In,Bn,Jinv,detJ,Minv,rBC,resids,maxres,&
    gamma,Rgas,CFL,convtol,min_iter,max_iter,nelem,niface,nbface)


    ! call getResidual(q,p,I2E,B2E,In,Bn,rBC,resids,Jinv,detJ,wavespeed,gamma,Rgas,nelem,niface,nbface)
end subroutine DG