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


subroutine getHOJacobian(nodes, geom, xy, J, Jinv, detJ)
    implicit none
    integer, intent(in) :: geom
    real(8), intent(in), dimension((geom+1)*(geom+2)/2,2) :: nodes
    real(8), intent(in), dimension(2) :: xy
    real(8), intent(out),dimension(2,2) :: J, Jinv
    real(8), intent(out) :: detJ
!f2py intent(in) nodes
!f2py intent(out) J, Jinv, detJ
    real(8), dimension((geom+1)*(geom+2)/2,2) :: gphi

    call gbasis2D(xy, geom, gphi, 1)
    J(:,:) = matmul(transpose(nodes(:,:)),gphi(:,:))
    detJ = J(1,1) * J(2,2) - J(2,1) * J(1,2)
    call matInv(2,J,Jinv)
end subroutine getHOJacobian