
      program main
          !!!!!!!!!!!!!!!!!!!!!!!
          !!!!   Author Yundi Quan
          !!!!
          !!!!!!!!!!!!!!!!!!!!!!!
          implicit none

          integer                       :: p,q,c1,c2,n1,n2
          integer                       :: c3,c4,il,iu,p_m,period
          integer                       :: lwork,lrwork,liwork,info,m,np,nq
          double precision              :: d,k1,k2,pi,zero,vl,vu,abstol
          double complex,allocatable    :: h(:,:),work(:),z(:,:)
          double precision,allocatable  :: rwork(:),w(:)
          integer,allocatable           :: isuppz(:),iwork(:)
          double precision              :: e_delta,tbar,t_delta,tprime,e3

          open(unit=2000,file='ene.txt')
          open(unit=2001,file='vec.txt')
          open(unit=2002,file='input')
          open(unit=2003,file='parameter.in')
          read(2002,*) n1,n2
          !read(2002,*) np,nq
          read(2003,*) e_delta,tbar,t_delta,tprime,e3

          pi    = 4e0*atan(1e0)
          abstol = 1d-6
          print *,pi

          write(2001,'(4I)') np,nq,n1,n2
          do q = 7,50
            do p = 0,q
              period = 3*q
              do c1=0,n1-1
                k1 = float(c1)/float(n1)*2d0*pi
                do c2=0,n2-1
                  k2 = float(c2)/float(n2)*2d0*pi
                  allocate(isuppz(8*period),source=0)
                  allocate(h(period,period),source=(0d0,0d0))
                  allocate(z(period,period))
                  allocate(work(1),rwork(1),iwork(1),w(period))
                  call find_h(k1,k2,h,p,q,e_delta,tbar,t_delta,tprime,e3)
                  call zheevr('N','A','U',period,h,period,vl,vu,il,iu,abstol,m,w,z,period,&
                      &isuppz,work,-1,rwork,-1,iwork,-1,info)
                  if (info .ne. 0) then
                    print *,'zheevr failed when calculating liwork etc'
                    print *,info
                  end if
                  lwork = int(work(1))
                  lrwork = int(rwork(1))
                  liwork = int(iwork(1))
                  deallocate(work,rwork,iwork,h,w)
                  allocate(work(lwork),rwork(lrwork),iwork(liwork),w(period))
                  allocate(h(period,period),source=(0d0,0d0))
                  call find_h(k1,k2,h,p,q,e_delta,tbar,t_delta,tprime,e3)
                  call zheevr('V','A','U',period,h,period,vl,vu,il,iu,abstol,m,w,z,period,&
                      &isuppz,work,lwork,rwork,lrwork,iwork,liwork,info)
                  do c3 = 1, period
                    write(2000,'(3I,3f25.10)') c3,p,q,float(p)/float(q),w(c3)
                  end do
                  do c3 = 1, 2*q
                    write(2001,'(2I,f25.10)') c1,c2,w(c3)
                    do c4 = 1, period
                        write(2001,'(2f25.10)') real(z(c4,c3)),aimag(z(c4,c3))
                    end do
                    write(2001,*)
                  end do
                  if (info .ne. 0) then
                    print *,'zheevr failed'
                  end if
                  deallocate(h,isuppz,work,rwork,iwork,z,w)
                enddo
              enddo
            enddo
          enddo
      close(2000)
      end

      subroutine find_h(kx,ky,h,p,q,e_delta,tbar,t_delta,tprime,e3)
          implicit none
          integer,intent(in)                :: p,q
          double complex,intent(inout)      :: h(3*q,3*q)
          double precision,intent(in)       :: kx,ky
          double precision,intent(in)       :: e_delta,t_delta,e3,tbar,tprime

          integer           :: i,j,m
          double precision  :: p_f,q_f,pi
          double precision  :: e_s1,e_s2,e_d,t_s1,t_s2
          double precision  :: lambda_1,lambda_2,phi

          p_f = dble(p)
          q_f = dble(q)
          pi  = 4d0*atan(1d0)

          phi = p_f/q_f*2d0*pi

          e_s1 = -e_delta/2d0
          e_s2 = e_delta/2d0
          e_d  = e3
          t_s1 = tbar-t_delta/2d0
          t_s2 = -(tbar+t_delta/2d0)
          lambda_1 = tprime
          lambda_2 = tprime

          do m=1,q-2
            h(3*m+1,3*m-2) = t_s1
            h(3*m+1,3*m)   = lambda_1
            h(3*m+1,3*m+1) = e_s1 + 2d0*t_s1*cos(ky-float(m)*phi)
            h(3*m+1,3*m+3) = -2d0*lambda_1*cos(ky-float(m)*phi)
            h(3*m+1,3*m+4) = t_s1
            h(3*m+1,3*m+6) = lambda_1

            h(3*m+2,3*m-1) = t_s2
            h(3*m+2,3*m)   = lambda_2
            h(3*m+2,3*m+2) = e_s2 + 2d0*t_s2*cos(ky-float(m)*phi)
            h(3*m+2,3*m+3) = -2d0*lambda_2*cos(ky-float(m)*phi)
            h(3*m+2,3*m+5) = t_s2
            h(3*m+2,3*m+6) = lambda_2

            h(3*m+3,3*m-2) = lambda_1
            h(3*m+3,3*m-1) = lambda_2
            h(3*m+3,3*m+1) = -2d0*lambda_2*cos(ky-float(m)*phi)
            h(3*m+3,3*m+2) = -2d0*lambda_2*cos(ky-float(m)*phi)
            h(3*m+3,3*m+3) = e_d
            h(3*m+3,3*m+4) = lambda_1
            h(3*m+3,3*m+5) = lambda_2
          end do

          h(1,1)      = e_s1+2d0*t_s1*cos(ky)
          h(1,3)      = -2d0*lambda_1*cos(ky)
          h(1,4)      = t_s1
          h(1,6)      = lambda_1
          h(1,3*q-2)  = t_s1*exp((0d0,-1d0)*kx)
          h(1,3*q)    = lambda_1*exp((0d0,-1d0)*kx)

          h(2,2)      = e_s2+2d0*t_s2*cos(ky)
          h(2,3)      = -2d0*lambda_2*cos(ky)
          h(2,5)      = t_s2
          h(2,6)      = lambda_2
          h(2,3*q)    = lambda_2*exp((0d0,-1d0)*kx)
          h(2,3*q-1)  = t_s2*exp((0d0,-1d0)*kx)

          h(3,3)      = e_d
          h(3,4)      = lambda_1
          h(3,1)      = -2d0*lambda_1*cos(ky)
          h(3,2)      = -2d0*lambda_2*cos(ky)
          h(3,3*q-2)  = lambda_1*exp((0d0,-1d0)*kx)
          h(3,3*q-1)  = lambda_2*exp((0d0,-1d0)*kx)
          h(3,5)      = lambda_2

          h(3*q-2,3*q-5) = t_s1
          h(3*q-2,3*q-3) = lambda_1
          h(3*q-2,3*q-2) = e_s1 + 2d0*t_s1*cos(ky-float(q-1)*phi)
          h(3*q-2,3*q)   = -2d0*lambda_1*cos(ky-float(q-1)*phi)
          h(3*q-2,1)     = t_s1*exp((0d0,1d0)*kx)
          h(3*q-2,3)     = lambda_1*exp((0d0,1d0)*kx)

          h(3*q-1,3*q-1) = e_s2 + 2d0*t_s2*cos(ky-float(q-1)*phi)
          h(3*q-1,2)     = t_s2*exp((0d0,1d0)*kx)
          h(3*q-1,3*q-4) = t_s2
          h(3*q-1,3)     = lambda_2 *exp((0d0,1d0)*kx)
          h(3*q-1,3*q)   = -2d0*lambda_2*cos(ky-float(q-1)*phi)
          h(3*q-1,3*q-3) = lambda_2

          h(3*q,3*q)     = e_d
          h(3*q,1)       = lambda_1*exp((0d0,1d0)*kx)
          h(3*q,3*q-5)   = lambda_1
          h(3*q,3*q-2)   = -2d0*lambda_1*cos(ky-float(q-1)*phi)
          h(3*q,2)       = lambda_2*exp((0d0,1d0)*kx)
          h(3*q,3*q-4)   = lambda_2
          h(3*q,3*q-1)   = -2d0*lambda_2*cos(ky-float(q-1)*phi)
      end subroutine
