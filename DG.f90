subroutine DG(q,p,resids,maxres,detJ,nodes,E2N,I2E,B2E,In,Bn,rBC,gamma,Rgas,CFL,convtol,min_iter,&
    max_iter,nnodes,nelem,niface,nbface)
    implicit none
    integer, intent(in) :: p,nnodes,nelem,niface,nbface,min_iter,max_iter
    real(8), intent(inout), dimension(nelem,(p+1)*(p+2)/2,4) :: q
    real(8), intent(in), dimension(nnodes,2) :: nodes
    integer, intent(in), dimension(nelem,3) :: E2N
    integer, intent(in), dimension(niface,4) :: I2E
    integer, intent(in), dimension(nbface,3) :: B2E
    real(8), intent(in), dimension(niface,3) :: In ! includes length as 3rd component
    real(8), intent(in), dimension(nbface,3) :: Bn ! includes length as 3rd component
    real(8), intent(out),dimension(nelem,(p+1)*(p+2)/2,4) :: resids
    real(8), intent(in), dimension(5) :: rBC
    real(8), intent(in) :: gamma, Rgas, convtol, CFL
    real(8), intent(out), dimension(max_iter) :: maxres
    real(8), dimension(nelem),intent(out) :: detJ
!f2py intent(in) node,E2N,I2E,B2E,In,Bn,gamma,rBC
!f2py intent(out) resids, maxres, detJ
!f2py intent(in,out) q

    integer :: ielem, Nb
    real(8), dimension(nelem,2,2) :: J, Jinv
    real(8), dimension(:,:,:), allocatable :: Minv
    real(8), dimension(nelem) :: wavespeed
    
    ! -----------------------------------
    ! Define constants and allocate arrays
    ! -----------------------------------
    Nb = (p+1)*(p+2)/2
    allocate(Minv(nelem,Nb,Nb))

    ! -----------------------------------
    ! Pre-compute quantities
    ! -----------------------------------

    do ielem = 1,nelem
        call getJacobian(nodes(E2N(ielem,:),:), J(ielem,:,:), Jinv(ielem,:,:), detJ(ielem))
    enddo
    call getMassInv(p,detJ,Minv,nelem)
    
    call timeIntegration(q,p,I2E,B2E,In,Bn,Jinv,detJ,Minv,rBC,resids,maxres,&
    gamma,Rgas,CFL,convtol,min_iter,max_iter,nelem,niface,nbface)


    ! call getResidual(q,p,I2E,B2E,In,Bn,rBC,resids,Jinv,detJ,wavespeed,gamma,Rgas,nelem,niface,nbface)
end subroutine DG