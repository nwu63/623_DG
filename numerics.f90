subroutine fullOrder(r,s,p,k)
    implicit none
    integer, intent(in) :: r,s,p
    integer, intent(out) :: k
    integer :: sp
    k = r + 1
    do sp = 0, s-1
        k = k + (p + 1 - sp)
    enddo
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


subroutine basis2D(xy, p, phi, n_xy)
    implicit none
    ! evaluates 2D basis functions at xy
    integer, intent(in) :: n_xy, p
    real(8), intent(in), dimension(n_xy,2) :: xy
    real(8), intent(out), dimension(n_xy,(p+1)*(p+2)/2) :: phi
!f2py intent(in) xy
!f2py intent(out) phi
    integer :: s, r, k, j, i_xy
    real(8), dimension((p+1)*(p+2)/2,(p+1)*(p+2)/2) :: C
    
    call triLagrange2D(p,C)
    phi(:,:) = 0.d0
    do j = 1,(p+1)*(p+2)/2
        do s = 0, p
            do r = 0, p-s
                call fullOrder(r,s,p,k)
                do i_xy = 1,n_xy
                    phi(i_xy,j) = phi(i_xy,j) + C(k,j) * xy(i_xy,1)**r * xy(i_xy,2)**s
                enddo
            enddo
        enddo
    enddo
end subroutine basis2D

subroutine gbasis2D(xy, p, gphi, n_xy)
    implicit none
    ! evaluates gradient of 2D basis functions at xy
    integer, intent(in) :: n_xy, p
    real(8), intent(in), dimension(n_xy,2) :: xy
    real(8), intent(out), dimension(n_xy,(p+1)*(p+2)/2,2) :: gphi
!f2py intent(in) xy
!f2py intent(out) gphi
    integer :: s, r, k, j, i_xy
    real(8), dimension((p+1)*(p+2)/2,(p+1)*(p+2)/2) :: C
    
    call triLagrange2D(p,C)
    gphi(:,:,:) = 0.d0
    do j = 1, (p+1)*(p+2)/2
        do s = 0, p
            do r = 0, p-s
                call fullOrder(r,s,p,k)
                do i_xy = 1,n_xy
                    gphi(i_xy,j,1) = gphi(i_xy,j,1) + r * C(k,j) * xy(i_xy,1)**(r-1) * xy(i_xy,2)**s
                    gphi(i_xy,j,2) = gphi(i_xy,j,2) + s * C(k,j) * xy(i_xy,1)**r * xy(i_xy,2)**(s-1)
                enddo
            enddo
        enddo
    enddo
end subroutine gbasis2D


subroutine matInv(n,A,Ainv)
    ! -----------------------------------------------------------------------
    ! Purpose: Calculates the matrix inverse of A using LAPACK's LU factorization
    ! 
    ! Inputs:
    !   n = the size of the system
    !   A[n,n] = the matrix to be inverted
    ! 
    ! Outs:
    !   Ainv[n,n] = the inverse of matrix A
    ! -----------------------------------------------------------------------
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
    call DGETRF(n,n,LU,n,ipiv,info)
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