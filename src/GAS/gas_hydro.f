      subroutine gas_hydro(tm)


      use gridmod
      use gasmod
      implicit none

      real*8,intent(in) :: tm



      real*8 x, y, z, ux0, uy0, uz0, rhoinv, tinv, dx, dy, dz
      real*8 ap, am, a0, dxinv
      integer i, j, k, l, lpx, lpy, lpz, lmx, lmy, lmz, fp, fm

      real*8,allocatable :: vel_x(:)
      real*8,allocatable :: vel_y(:)
      real*8,allocatable :: vel_z(:)
      real*8,allocatable :: sx_du(:)
      real*8,allocatable :: sy_du(:)
      real*8,allocatable :: sz_du(:)
c
      allocate(vel_x(gas_ncell))
      allocate(vel_y(gas_ncell))
      allocate(vel_z(gas_ncell))
      allocate(sx_du(gas_ncell))
      allocate(sy_du(gas_ncell))
      allocate(sz_du(gas_ncell))

      do i=1,grd_nx
       do j=1,grd_ny
        do k=1,grd_nz
         l = grd_icell(i,j,k)
         lpx = grd_icell(i+1,j,k)
         lpy = grd_icell(i,j+1,k)
         lpz = grd_icell(i,j,k+1)
         rhoinv = 1.0d0 / gas_rho(l)
         tinv = 1.0d0 / tm
         x = (grd_xarr(l) + grd_xarr(lpx)) * 0.5d0
         y = (grd_xarr(l) + grd_yarr(lpy)) * 0.5d0
         z = (grd_xarr(l) + grd_zarr(lpz)) * 0.5d0
         ux0 = x * tinv
         uy0 = y * tinv
         uz0 = z * tinv
         vel_x(l) = gas_sx(l) * rhoinv - ux0
         vel_y(l) = gas_sy(l) * rhoinv - uy0
         vel_z(l) = gas_sz(l) * rhoinv - uz0
        enddo
       enddo
      enddo

      do i=2,grd_nx-1
       do j=2,grd_ny-1
        do k=2,grd_nz-1
         l = grd_icell(i,j,k)
         lpx = grd_icell(i+1,j,k)
         lmx = grd_icell(i-1,j,k)
         dx = (grd_xarr(l) - grd_xarr(lpx))
         dxinv = 1.0d0 / dx
         ap = abs(vel_x(lpx))
         a0 = abs(vel_x(l))
         am = abs(vel_x(lmx))
         a0 = max(ap,a0,am)

         fp = (gas_sx(lpx) * vel_x(lpx) + gas_sx(l) * vel_x(l))*0.5d0
         fm = (gas_sx(lmx) * vel_x(lmx) + gas_sx(l) * vel_x(l))*0.5d0
         fp = fp - 0.5d0 * a0 * (gas_sx(lpx) - gas_sx(l))
         fm = fm - 0.5d0 * a0 * (gas_sx(l) - gas_sx(lmx))
         sx_du(l) = (fp - fm) * dxinv

         fp = (gas_sy(lpx) * vel_x(lpx) + gas_sy(l) * vel_x(l))*0.5d0
         fm = (gas_sy(lmx) * vel_x(lmx) + gas_sy(l) * vel_x(l))*0.5d0
         fp = fp - 0.5d0 * a0 * (gas_sy(lpx) - gas_sy(l))
         fm = fm - 0.5d0 * a0 * (gas_sy(l) - gas_sy(lmx))
         sy_du(l) = (fp - fm) * dxinv

         fp = (gas_sz(lpx) * vel_x(lpx) + gas_sz(l) * vel_x(l))*0.5d0
         fm = (gas_sz(lmx) * vel_x(lmx) + gas_sz(l) * vel_x(l))*0.5d0
         fp = fp - 0.5d0 * a0 * (gas_sz(lpx) - gas_sz(l))
         fm = fm - 0.5d0 * a0 * (gas_sz(l) - gas_sz(lmx))
         sz_du(l) = (fp - fm) * dxinv

        enddo
       enddo
      enddo



      deallocate(vel_x)
      deallocate(vel_y)
      deallocate(vel_z)
      deallocate(sx_du)
      deallocate(sy_du)
      deallocate(sz_du)

      end subroutine gas_hydro
