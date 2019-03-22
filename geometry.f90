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