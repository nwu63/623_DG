subroutine getRefMassMatrix(p,M)
    implicit none
    integer, intent(in) :: p
    real(8), intent(out), dimension((p+1)*(p+2)/2,(p+1)*(p+2)/2) :: M
!f2py intent(in) p
!f2py intent(out) M
    integer :: Nb, Ng, i,g
    real(8), dimension(:), allocatable :: w
    real(8), dimension(:,:), allocatable :: xy, phi, wmat

    Nb = (p+1)*(p+2)/2

    g = 2*p+1
    call Gauss2D_pre(g,Ng) ! integrand is order 2p + 2*(q-1) but for now we assume q = 1 (linear elements)
    allocate(xy(Ng,2))
    allocate(w(Ng))
    allocate(phi(Ng,Nb))
    allocate(wmat(Ng,Ng))
    call Gauss2D(g,Ng,xy,w)
    wmat(:,:) = 0.d0
    FORALL(i=1:Ng) wmat(i,i) = w(i) ! wmat = diag(w)
    call basis2D(xy, p, phi, Ng)
    M = matmul(matmul(transpose(phi),wmat),phi)
end subroutine getRefMassMatrix

subroutine getMassInv(p,detJ,Minv,nelem)
    implicit none
    integer, intent(in) :: p,nelem
    real(8), intent(in), dimension(nelem) :: detJ
    real(8), intent(out), dimension(nelem,(p+1)*(p+2)/2,(p+1)*(p+2)/2) :: Minv
!f2py intent(in) p,nelem,detJ
!f2py intetn(out) Minv
    integer :: Nb, ielem
    real(8), dimension(:,:), allocatable :: Mref, Mrefinv

    Nb = (p+1)*(p+2)/2
    allocate(Mref(Nb,Nb))
    allocate(Mrefinv(Nb,Nb))
    call getRefMassMatrix(p,Mref)
    call matInv(Nb,Mref,Mrefinv)

    do ielem = 1, nelem
        Minv(ielem, :,:) = 1/detJ(ielem)*Mrefinv
    enddo
end subroutine getMassInv


