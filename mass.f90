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

subroutine getHORefMassMatrix(p,w,phi,detJ2,M,Ng)
    implicit none
    integer, intent(in) :: p,Ng
    real(8), intent(in), dimension(Ng) :: w
    real(8), intent(in), dimension(Ng,(p+1)*(p+2)/2) :: phi
    real(8), intent(in), dimension(Ng) :: detJ2
    real(8), intent(out), dimension((p+1)*(p+2)/2,(p+1)*(p+2)/2) :: M
!f2py intent(in) p
!f2py intent(out) M
    integer :: Nb, i
    real(8), dimension(:,:), allocatable :: wmat

    Nb = (p+1)*(p+2)/2
    allocate(wmat(Ng,Ng))
    
    wmat(:,:) = 0.d0
    FORALL(i=1:Ng) wmat(i,i) = w(i)*detJ2(i) ! we multiply the detJ weights here in the integration
    M = matmul(matmul(transpose(phi),wmat),phi)
end subroutine getHORefMassMatrix

subroutine getMassInv(p,detJ,detJ2,w,phi,qlist,Minv,nelem,nqelem,Ng)
    implicit none
    integer, intent(in) :: p,nelem,nqelem,Ng
    real(8), intent(in), dimension(nelem) :: detJ
    real(8), intent(in), dimension(nqelem,Ng) :: detJ2
    real(8), intent(in), dimension(nqelem) :: qlist
    real(8), intent(out), dimension(nelem,(p+1)*(p+2)/2,(p+1)*(p+2)/2) :: Minv
    real(8), intent(in), dimension(Ng) :: w
    real(8), intent(in), dimension(Ng,(p+1)*(p+2)/2) :: phi
!f2py intent(in) p,nelem,detJ
!f2py intetn(out) Minv
    integer :: Nb, ielem, idx
    real(8), dimension(:,:), allocatable :: Mref, Mrefinv, M

    Nb = (p+1)*(p+2)/2
    allocate(Mref(Nb,Nb))
    allocate(M(Nb,Nb))
    allocate(Mrefinv(Nb,Nb))
    call getRefMassMatrix(p,Mref)
    call matInv(Nb,Mref,Mrefinv)

    do ielem = 1, nelem
        if (any(qlist == ielem)) then ! curved element
            idx = minloc((qlist-ielem)**2,dim=1)
            call getHORefMassMatrix(p,w,phi,detJ2(idx,:),M,Ng)
            call matInv(Nb,M,Minv(ielem,:,:))
        else
            Minv(ielem, :,:) = 1/detJ(ielem)*Mrefinv
        endif
    enddo
end subroutine getMassInv


