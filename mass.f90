subroutine getRefMassMatrix(p,M)
    integer, intent(in) :: p
    real(8), intent(out), dimension((p+1)*(p+2)/2,(p+1)*(p+2)/2) :: M

    integer :: Nb, Ng, i
    real(8), dimension(:), allocatable :: xy, w
    real(8), dimension(:,:), allocatable :: phi, wmat

    Nb = (p+1)*(p+2)/2

    call Gauss2D_pre(2*p,Ng) ! integrand is order 2p + 2*(q-1) but for now we assume q = 1 (linear elements)
    allocate(xy(2*Ng))
    allocate(w(Ng))
    allocate(phi(Ng,Nb))
    allocate(wmat(Ng,Ng))
    call Gauss2D(2*p,Ng,xy,w)
    wmat = 0 
    FORALL(i=1:size(wmat,1)) wmat(i,i) = w(i) ! wmat = diag(w)
    call basis2D(xy, p, phi, Ng)
    M = matmul(matmul(transpose(phi),wmat),phi)
end subroutine getRefMassMatrix