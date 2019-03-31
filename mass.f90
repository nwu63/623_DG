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

subroutine getHORefMassMatrix(p,geom,nodes,M)
    implicit none
    integer, intent(in) :: p,geom
    real(8), intent(in), dimension((geom+1)*(geom+2)/2,2) :: nodes
    real(8), intent(out), dimension((p+1)*(p+2)/2,(p+1)*(p+2)/2) :: M
!f2py intent(in) p
!f2py intent(out) M
    integer :: Nb, Ng, i,g
    real(8), dimension(:), allocatable :: w, detJ
    real(8), dimension(:,:), allocatable :: xy, phi, wmat
    real(8), dimension(2,2) :: J, Jinv

    Nb = (p+1)*(p+2)/2

    g = 2*p+1
    call Gauss2D_pre(g,Ng) ! integrand is order 2p + 2*(q-1) but for now we assume q = 1 (linear elements)
    allocate(xy(Ng,2))
    allocate(w(Ng))
    allocate(phi(Ng,Nb))
    allocate(wmat(Ng,Ng))
    call Gauss2D(g,Ng,xy,w)

    allocate(detJ(Ng))
    do i = 1,Ng
        call getHOJacobian(nodes, geom, xy(i,:), J, Jinv, detJ(i))
        print*, detJ(i)
    enddo


    wmat(:,:) = 0.d0
    FORALL(i=1:Ng) wmat(i,i) = w(i)*detJ(i) ! we multiply the detJ weights here in the integration
    call basis2D(xy, p, phi, Ng)
    M = matmul(matmul(transpose(phi),wmat),phi)
end subroutine getHORefMassMatrix

subroutine getMassInv(nodes,E2N2,p,geom,qlist,detJ,Minv,nelem,nnodes,nqelem)
    implicit none
    integer, intent(in) :: p,nelem,geom, nqelem, nnodes
    real(8), intent(in), dimension(nelem) :: detJ
    real(8), intent(in), dimension(nnodes,2) :: nodes
    integer, intent(in), dimension(nqelem,(geom+1)*(geom+2)/2) :: E2N2
    integer, intent(in), dimension(nqelem) :: qlist
    real(8), intent(out), dimension(nelem,(p+1)*(p+2)/2,(p+1)*(p+2)/2) :: Minv
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
            call getHORefMassMatrix(p,geom,nodes(E2N2(idx,:),:),M)
            call matInv(Nb,M,Minv(ielem,:,:))
        else
            Minv(ielem, :,:) = 1/detJ(ielem)*Mrefinv
        endif
    enddo
end subroutine getMassInv


