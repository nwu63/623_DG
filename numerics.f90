subroutine fullOrder(r,s,p,k)
    implicit none
    integer, intent(in) :: r,s,p
    integer, intent(out) :: k
    integer :: sp
    k = 0
    do sp = 0, s-1
        k = k + (p + 1 - sp)
    enddo
    k = k + r + 1
end subroutine fullOrder

subroutine triLagrange2D(p,C)
    implicit none
    integer, intent(in) :: p
    integer :: i,ix,iy,k,s,r,N
    real(8), dimension((p+1)*(p+2)/2,(p+1)*(p+2)/2) :: A,C
    real(8), dimension(p+1) :: xi, eta
    N = (p+1)*(p+2)/2
    xi = (/(i, i=0,p+1,1)/)
    xi = xi/p
    eta = xi
    A = 0.d0
    C = A
    i = 1
    do iy = 0, p
        do ix = 0, p-iy
            k = 1
            do s = 0, p
                do r = 0, p-s
                    A(i,k) = xi(ix+1)**r * eta(iy+1)**s
                    k = k+1
                enddo
            enddo
            i = i + 1
        enddo
    enddo
    call matInv(N,A,C)
end subroutine triLagrange2D


subroutine lagrange(xn,j, x, phi, n_xn, n_x)
    ! computes values of jth lagrange function based on nodes
    ! in xn, at the x-values given in x

    integer, intent(in) :: j, n_xn, n_x
    real(8), intent(in), dimension(n_xn) :: xn
    real(8), intent(in), dimension(n_x) :: x
    real(8), intent(out), dimension(n_x) :: phi
!f2py intent(in) xn,x
!f2py intent(out) phi
    real(8), dimension(n_xn-1) :: xnj
    integer, dimension(n_xn-1) :: idx
    integer :: i
    real(8) :: den
    real(8), dimension(n_x) :: num

    if (n_xn == 1) then
        phi(:) = 1.d0
        return
    endif

    idx = (/(i,i=1,n_xn-1)/)
    where (idx>=j)
        idx = idx + 1
    endwhere

    xnj = xn(idx)
    den = product(xn(j)-xnj)
    num = product(spread(x,2,n_xn-1) - spread(xnj,1,n_x), 2)
    print*, den,num
    phi = num/den
end subroutine lagrange

subroutine basis1D(xn, x, phi, n_xn, n_x)
    ! evaluates basis functions at x
    integer, intent(in) :: n_xn, n_x
    real(8), intent(in), dimension(n_xn) :: xn
    real(8), intent(in), dimension(n_x) :: x
    real(8), intent(out), dimension(n_x,n_xn) :: phi
!f2py intent(in) xn,x
!f2py intent(out) phi
    integer :: p
    real(8), dimension(n_x) :: B
    
    do p = 0,n_xn-1
        call lagrange(xn,p+1, x, B, n_xn, n_x)
        Phi(:,p+1) = B
    enddo
end subroutine basis1D

subroutine basis2D(xy, p, phi, n_xy)
    ! evaluates 2D basis functions at xy
    integer, intent(in) :: n_xy, p
    real(8), intent(in), dimension(n_xy,2) :: xy
    real(8), intent(out), dimension(n_xy,(p+1)*(p+2)/2) :: phi
!f2py intent(in) xn,x
!f2py intent(out) phi
    integer :: s, r, k, j, i_xy
    real(8), dimension((p+1)*(p+2)/2,(p+1)*(p+2)/2) :: C
    
    call triLagrange2D(p,C)
    j = 1
    do s = 0, p
        do r = 0, p-s
            call fullOrder(r,s,p,k)
            do i_xy = 1,n_xy
                phi(i_xy,j) = C(k,j) * xy(i_xy,1)**r * xy(i_xy,2)**s
            enddo
            j = j + 1
        enddo
    enddo
end subroutine basis2D

subroutine gbasis2D(xy, p, gphi, n_xy)
    ! evaluates gradient of 2D basis functions at xy
    integer, intent(in) :: n_xy, p
    real(8), intent(in), dimension(n_xy,2) :: xy
    real(8), intent(out), dimension(n_xy,(p+1)*(p+2)/2,2) :: gphi
!f2py intent(in) xn,x
!f2py intent(out) gphi
    integer :: s, r, k, j, i_xy
    real(8), dimension((p+1)*(p+2)/2,(p+1)*(p+2)/2) :: C
    
    call triLagrange2D(p,C)
    j = 1
    do s = 0, p
        do r = 0, p-s
            call fullOrder(r,s,p,k)
            do i_xy = 1,n_xy
                gphi(i_xy,j,1) = r * C(k,j) * xy(i_xy,1)**(r-1) * xy(i_xy,2)**s
                gphi(i_xy,j,2) = s * C(k,j) * xy(i_xy,1)**r * xy(i_xy,2)**(s-1)
            enddo
            j = j + 1
        enddo
    enddo
end subroutine gbasis2D


subroutine matInv(n,A,Ainv)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in),dimension(n,n) :: A
    real(8), intent(out),dimension(n,n) :: Ainv
!f2py intent(in) n, A
!f2py intent(out) Ainv
    real(8), dimension(n) :: work
    real(8), dimension(n,n) :: LU
    integer :: info
    integer, dimension(n) :: ipiv
    LU = A
    CALL DGETRF(n,n,LU,n,ipiv,info)
    Ainv = LU
    call DGETRI(n,Ainv,n,ipiv,work,n,info)
end subroutine matInv

subroutine solveLinSys(A,x,b,n)
    implicit none
    real(8), intent(in), dimension(n,n) :: A
    real(8), intent(in), dimension(n) :: b
    real(8), intent(out), dimension(n) :: x
    integer, intent(in) :: n
    integer :: ipiv,info
    call DGESV(n,1,A,n,ipiv,b,n,info)
    x = b
end subroutine solveLinSys