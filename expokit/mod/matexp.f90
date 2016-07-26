module dexpm_wrap

  use iso_c_binding
  implicit none
      
contains
  
  subroutine wrapalldmexpv(n,m,t, v, w, tol, anorm, wsp, lwsp, iwsp, liwsp, itrace, iflag,ia,ja,a,nz,res ) bind(c)
    integer,intent(inout) :: m, lwsp, liwsp, itrace, iflag
    integer,intent(inout) :: iwsp(liwsp)
    double precision,intent(inout) :: t, tol, anorm, v(n), w(n)
    double precision,intent(inout) :: wsp(lwsp)
    integer,intent(in) :: nz,n,ia(nz),ja(nz)
    double precision,intent(in) :: a(nz)
    double precision, intent(out) :: res(n*n)
    
    double precision ZERO, ONE
    parameter( ZERO=0.0d0, ONE=1.0d0 )
    intrinsic ABS
    integer i,j
    
    do i = 1,n
       wsp(i) = ZERO
    enddo
    do i = 1,nz
       wsp(ia(i)) = wsp(ia(i)) + ABS( a(i) )
    enddo
    anorm = wsp(1)
    do i = 2,n
       if ( anorm.lt.wsp(i) ) anorm =  wsp(i)
    enddo
    do i = 1,n
       do j = 1,n
          v(j) = ZERO
       enddo
       v(i) = ONE
       call myDMEXPV( n, m, t,v,w, tol, anorm,wsp,lwsp, iwsp,liwsp, itrace,iflag,ia,ja,a,nz )
       do j = 1,n
          res(((i-1)*n)+j) = w(j)
       enddo
    enddo
  end subroutine wrapalldmexpv
      
  subroutine wrapsingledmexpv(n,t,v,w,wsp,lwsp,iwsp,liwsp,ia,ja,a,nz) bind(c)
    integer(c_int), intent(in), value :: n, nz, lwsp, liwsp
    real(c_double), intent(in), value :: t
    integer(c_int), intent(in) :: iwsp(liwsp)
    real(c_double) :: v(n), wsp(lwsp)
    integer(c_int), intent(in) :: ia(nz), ja(nz)
    real(c_double), intent(in) :: a(nz)
    real(c_double), intent(out) :: w(n)

    real(c_double) ZERO, ONE, anorm, tol
    parameter( ZERO=0.0d0, ONE=1.0d0 )
    intrinsic ABS
    integer(c_int) i, j, iflag
    
    do i = 1,n
       wsp(i) = ZERO
    enddo

    do i = 1,nz
       wsp(ia(i)) = wsp(ia(i)) + ABS( a(i) )
    enddo
    
    anorm = wsp(1)
    do i = 2,n
       if ( anorm.lt.wsp(i) ) anorm =  wsp(i)
    enddo
    tol = ZERO
    iflag = 0
    call myDMEXPV(n, n-1, t, v, w, tol, anorm, wsp, lwsp, iwsp, liwsp, 0, iflag, ia, ja, a, nz)
    end subroutine wrapsingledmexpv

    subroutine wrapdgpadm(ideg,m,t,H,ldh,wsp,lwsp,ipiv,iexph, ns,iflag ) bind(c)
      integer,intent(inout) :: ideg,m,ldh,lwsp,iexph,ns,iflag,ipiv(m)
      double precision,intent(inout) :: t,H(ldh,m),wsp(lwsp)
      integer i,j
      call DGPADM( ideg,m,t,H,ldh,wsp,lwsp,ipiv,iexph,ns,iflag)
    end subroutine wrapdgpadm

end module
