subroutine getJacobian(nodes, J, Jinv, detJ)
    implicit none
    real(8), intent(in), dimension(3,2) :: nodes
    real(8), intent(out),dimension(2,2) :: J, Jinv
    real(8), intent(out) :: detJ
!f2py intent(in) nodes
!f2py intent(out) J, Jinv, detJ
    J(1,1) = nodes(2,1) - nodes(1,1)
    J(1,2) = nodes(3,1) - nodes(1,1)
    J(2,1) = nodes(2,2) - nodes(1,2)
    J(2,2) = nodes(3,2) - nodes(1,2)
    detJ = J(1,1) * J(2,2) - J(2,1) * J(1,2)
    call matInv(2,J,Jinv)
end subroutine getJacobian

subroutine getX(nodes,geom,xy,phys,n_xy)
    implicit none
    integer, intent(in) :: geom,n_xy
    real(8), intent(in), dimension((geom+1)*(geom+2)/2,2) :: nodes
    real(8), intent(in), dimension(n_xy,2) :: xy
    real(8), intent(out), dimension(n_xy,2) :: phys
!f2py intent(in) nodes,xy
!f2py intent(out) phys
    real(8), dimension(n_xy,(geom+1)*(geom+2)/2) :: phi

    call basis2D(xy, geom, phi, n_xy)
    phys = matmul(phi,nodes)
end subroutine getX

subroutine getHOJacobian(nodes, geom, xy, J, Jinv, detJ, n_xy)
    implicit none
    integer, intent(in) :: geom, n_xy
    real(8), intent(in), dimension((geom+1)*(geom+2)/2,2) :: nodes
    real(8), intent(in), dimension(n_xy,2) :: xy
    real(8), intent(out),dimension(n_xy,2,2) :: J, Jinv
    real(8), intent(out),dimension(n_xy) :: detJ
!f2py intent(in) nodes
!f2py intent(out) J, Jinv, detJ
    real(8), dimension(n_xy,(geom+1)*(geom+2)/2,2) :: gphi
    integer :: ig

    call gbasis2D(xy, geom, gphi, n_xy)
    do ig = 1,n_xy
        J(ig,:,:) = matmul(transpose(nodes),gphi(ig,:,:))
        detJ(ig) = J(ig,1,1) * J(ig,2,2) - J(ig,2,1) * J(ig,1,2)
        call matInv(2,J(ig,:,:),Jinv(ig,:,:))
    enddo
end subroutine getHOJacobian

! subroutine getEdgeNrmJac()
! end subroutine getEdgeNrmJac

subroutine getdXidSigma(face,vec)
    implicit none
    integer, intent(in) :: face
    real(8), intent(out), dimension(2) :: vec

    if (face == 1) then
        vec = (/-1.d0,1.d0/)
    elseif (face == 2) then
        vec = (/0.d0,-1.d0/)
    elseif (face == 3) then
        vec = (/1.d0,0.d0/)
    endif
end subroutine getdXidSigma