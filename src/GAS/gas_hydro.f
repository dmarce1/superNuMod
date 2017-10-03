      subroutine gas_hydro(tm,dt)


      use gridmod
      use gasmod
      implicit none

      real*8,intent(in) :: tm
      real*8,intent(in) :: dt



      real*8 x(3), rhoinv, tinv, dx(3)
      real*8 ap, am, a0, dxinv, tmp
      integer i, j, k, l, lp(3), lm(3), fp, fm, d, d1, d2

      real*8,allocatable :: vel(:,:)
      real*8,allocatable :: mom_du(:,:), rho_du(:)
c
      allocate(vel(3,gas_ncell))
      allocate(mom_du(3,gas_ncell))
      allocate(rho_du(gas_ncell))

      if( grd_isvelocity ) then
       tinv = 1.0d0 / tm
      else
       tinv = 1.0d0
      endif

      do i=1,grd_nx
       do j=1,grd_ny
        do k=1,grd_nz
         l = grd_icell(i,j,k)
         lp(1) = grd_icell(i+1,j,k)
         lp(2) = grd_icell(i,j+1,k)
         lp(3) = grd_icell(i,j,k+1)
         rhoinv = 1.0d0 / gas_rho(l)
         x(1) = (grd_xarr(l) + grd_xarr(lp(1))) * 0.5d0*tinv
         x(2) = (grd_xarr(l) + grd_yarr(lp(2))) * 0.5d0*tinv
         x(3) = (grd_xarr(l) + grd_zarr(lp(3))) * 0.5d0*tinv
         do d=1,3
          tmp = x(d) * tinv
          vel(l,d) = gas_mom(l,d) * rhoinv - tmp
         enddo
        enddo
       enddo
      enddo

      do i=2,grd_nx-1
       do j=2,grd_ny-1
        do k=2,grd_nz-1
         l = grd_icell(i,j,k)
         lp(1) = grd_icell(i+1,j,k)
         lp(2) = grd_icell(i,j+1,k)
         lp(3) = grd_icell(i,j,k+1)
         lm(1) = grd_icell(i-1,j,k)
         lm(2) = grd_icell(i,j-1,k)
         lm(3) = grd_icell(i,j,k-1)
         dx(1) = (grd_xarr(lp(1)) - grd_xarr(l))*tinv
         dx(2) = (grd_xarr(lp(2)) - grd_xarr(l))*tinv
         dx(3) = (grd_xarr(lp(3)) - grd_xarr(l))*tinv
         dxinv = 1.0d0 / dx(d1)
         do d1=1,3
          mom_du(d1,l) =  - 3.0d0 * gas_mom(l,d1) * tinv
         enddo
         do d1=1,3
          ap = abs(vel(d1,lp(d1)))
          a0 = abs(vel(d1,l))
          am = abs(vel(d1,lm(d1)))
          a0 = max(ap,a0,am)
          tmp = gas_rho(l) * vel(d1,l)
          fp = (gas_rho(lp(d1)) * vel(lp(d1),l) + tmp)*0.5d0
          fm = (gas_rho(lm(d1)) * vel(lm(d1),l) + tmp)*0.5d0
          fp = fp - 0.5d0 * a0 * (gas_rho(lp(d1)) - gas_rho(l))
          fm = fm - 0.5d0 * a0 * (gas_rho(l) - gas_rho(lm(d1)))
          rho_du(l) = rho_du(l) - (fp - fm) * dxinv
          do d2=1,3
           tmp = gas_mom(d2,l) * vel(d1,l)
           fp = (gas_mom(d2,lp(d1)) * vel(d1,lp(d1)) + tmp)*0.5d0
           fm = (gas_mom(d2,lm(d1)) * vel(d1,lm(d1)) + tmp)*0.5d0
           fp = fp - 0.5d0 * a0 * (gas_mom(d2,lp(d1)) - gas_mom(d2,l))
           fm = fm - 0.5d0 * a0 * (gas_mom(d2,l) - gas_mom(d2,lm(d1)))
           mom_du(d2,l) = mom_du(d2,l) - (fp - fm) * dxinv
          enddo
         enddo


        enddo
       enddo
      enddo

      do l=1,gas_ncell
       do d=1,3
        gas_mom(d,l) = gas_mom(d,l) + dt * mom_du(d,l)
        gas_rho(l) = gas_rho(l) + dt * rho_du(l)
       enddo
      enddo


      deallocate(vel)
      deallocate(mom_du)

      end subroutine gas_hydro
