      subroutine gas_hydro(tm,dt)


      use gridmod
      use gasmod
      implicit none

      real*8,intent(in) :: tm
      real*8,intent(in) :: dt



      real*8 x(3), rhoinv, tinv, dx(3)
      real*8 ap, am, a0, dxinv, tmp
      integer i, j, k, l, m, lp(3), lm(3), fp, fm, d, d1, d2

      real*8,allocatable :: vel(:,:)
      real*8,allocatable :: mom_du(:,:), rho_du(:), natom_du(:,:)
c
      allocate(vel(3,gas_ncell))
      allocate(mom_du(3,gas_ncell))
      allocate(rho_du(gas_ncell))
      allocate(natom_du(-2*gas_nchain:gas_nelem,gas_ncell))

      if( grd_isvelocity ) then
       tinv = 1.0d0 / tm
      else
       tinv = 1.0d0
      endif

      do m=-2*gas_nchain,gas_nelem
       if(m.eq.0) then
        continue
       endif
       gas_natom1fr(m,:) = gas_natom1fr(m,:)*gas_natom(:)
      enddo

      do i=1,grd_nx
       x(1) = (grd_xarr(i) + grd_xarr(i+1)) * 0.5d0*tinv
       do j=1,grd_ny
        x(2) = (grd_xarr(j) + grd_yarr(j+1)) * 0.5d0*tinv
        do k=1,grd_nz
         l = grd_icell(i,j,k)
         rhoinv = 1.0d0 / gas_rho(l)
         x(3) = (grd_xarr(k) + grd_zarr(k+1)) * 0.5d0*tinv
         do d=1,3
          vel(d,l) = gas_mom(l,d) * rhoinv
          if( grd_isvelocity ) then
           vel(d,l) = vel(d,l) - x(d) * tinv
          endif
         enddo
        enddo
       enddo
      enddo

      do i=2,grd_nx-1
       dx(1) = (grd_xarr(i+1) - grd_xarr(i))*tinv
       do j=2,grd_ny-1
        dx(2) = (grd_yarr(j+1) - grd_yarr(j))*tinv
        do k=2,grd_nz-1
         l = grd_icell(i,j,k)
         lp(1) = grd_icell(i+1,j,k)
         lp(2) = grd_icell(i,j+1,k)
         lp(3) = grd_icell(i,j,k+1)
         lm(1) = grd_icell(i-1,j,k)
         lm(2) = grd_icell(i,j-1,k)
         lm(3) = grd_icell(i,j,k-1)
         dx(3) = (grd_zarr(k+1) - grd_zarr(k))*tinv
         if( grd_isvelocity ) then
          do d1=1,3
           mom_du(d1,l) =  - 3.0d0 * gas_mom(l,d1) * tinv
          enddo
          rho_du(l) =  - 3.0d0 * gas_rho(l) * tinv
         endif
         do d1=1,3
          dxinv = 1.0d0 / dx(d1)
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
          do m=-2*gas_nchain,gas_nelem
           if(m.eq.0) then
            continue
           endif
           tmp = gas_natom1fr(m,l) * vel(d1,l)
           fp = (gas_natom1fr(m,lp(d1)) * vel(d1,lp(d1)) + tmp)*0.5d0
           fm = (gas_natom1fr(m,lm(d1)) * vel(d1,lm(d1)) + tmp)*0.5d0
           fp = fp-0.5d0*a0*(gas_natom1fr(m,lp(d1)) - gas_natom1fr(m,l))
           fm = fm-0.5d0*a0*(gas_natom1fr(m,l) - gas_natom1fr(m,lm(d1)))
           natom_du(m,l) = natom_du(m,l) - (fp - fm) * dxinv
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

      gas_natom = sum(gas_natom1fr, 1)

      do m=-2*gas_nchain,gas_nelem
       if(m.eq.0) then
        continue
       endif
       gas_natom1fr(m,:) = gas_natom1fr(m,:)/gas_natom(:)
      enddo


      deallocate(vel)
      deallocate(mom_du)
      deallocate(natom_du)

      end subroutine gas_hydro
