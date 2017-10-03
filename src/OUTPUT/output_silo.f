
      subroutine putvar( db, var, name, namelen )

      use gasmod
      use gridmod
      implicit none

      include '/usr/include/silo.inc'

      integer, intent(in) :: db, namelen
      character(len=namelen) :: name
      real*8, intent(in)  :: var(gas_ncell)
      real*8, allocatable :: v3d(:,:,:)


      integer :: rc, i, j, k, l
      integer :: ndims = 3
      integer, dimension(3) :: dims
      integer :: ctt = DB_NODECENT
      integer :: datat = DB_DOUBLE



      dims(1) = grd_nx
      dims(2) = grd_ny
      dims(3) = grd_nz

      allocate(v3d(grd_nx,grd_ny,grd_nz))

      do i=1,grd_nx
       do j=1,grd_ny
        do k=1,grd_nz
         l = grd_icell(i,j,k)
         v3d(i,j,k) = var(l)
        enddo
       enddo
      enddo

      rc = DBPutQV1(db,name,'mesh',v3d,dims,ndims,0,0,datat,ctt,0)

      deallocate(v3d)
      end subroutine putvar


      subroutine output_silo(oname, onamelen, tm)

      use gridmod
      use gasmod
      use elemdatamod
      implicit none

      integer, intent(in) :: onamelen
      real*8, intent(in) :: tm
      character(len=onamelen), intent(in) :: oname

      include '/usr/include/silo.inc'
      integer :: db, rc, i
      real*8 :: tinv

      character(len=4) :: mname = 'mesh'
      character(len=1), dimension(3) :: cnames
      integer, dimension(3) :: dims
      integer :: datat = DB_DOUBLE
      integer :: crdt = DB_COLLINEAR
      integer :: ndims = 3
      real*8 :: coords(3,grd_nx)

      if( grd_isvelocity ) then
       tinv = 1.0d0 / tm
      else
       tinv = 1.0d0
      endif

      dims(1) = grd_nx
      dims(2) = grd_ny
      dims(3) = grd_nz
      do i=1,grd_nx
       coords(1,i) = (grd_xarr(i) + grd_xarr(i+1)) * 0.5d0 * tinv
      enddo
      do i=1,grd_ny
       coords(2,i) = (grd_yarr(i) + grd_yarr(i+1)) * 0.5d0 * tinv
      enddo
      do i=1,grd_nz
       coords(3,i) = (grd_zarr(i) + grd_zarr(i+1)) * 0.5d0 * tinv
      enddo
      cnames(1) = 'x'
      cnames(2) = 'y'
      cnames(3) = 'z'

      db = DBCreate(oname, DB_CLOBBER, DB_LOCAL, 'Euler Mesh', DB_PDB)

      rc = DBPutQM(db,mname,cnames,coords,dims,ndims,datat,crdt,0)
      do i=1,elem_neldata
       call putvar(db,gas_natom1fr(i,:),elem_data(i)%sym,2)
      enddo
      call putvar(db,gas_natom1fr(-1,:),'ni56',4)
      call putvar(db,gas_natom1fr(-2,:),'co56',4)
      call putvar(db,gas_natom1fr(-3,:),'fe52',4)
      call putvar(db,gas_natom1fr(-4,:),'mn52',4)
      call putvar(db,gas_natom1fr(-5,:),'cr48',4)
      call putvar(db,gas_natom1fr(-6,:),'v48',3)
      call putvar(db,gas_temp,'temp',4)
      call putvar(db,gas_eraddens,'eraddens',8)
      call putvar(db,gas_ur,'ur',2)
      call putvar(db,gas_rho,'rho',3)
      call putvar(db,gas_bcoef,'bcoef',5)
      call putvar(db,gas_decaygamma,'decaygamma',10)
      call putvar(db,gas_decaybeta,'decaybeta',9)
      call putvar(db,gas_ye,'Ye',2)
      call putvar(db,gas_natom,'natom',5)
      call putvar(db,gas_nelec,'nelec',5)
      call putvar(db,gas_matsrc,'matsrc',6)
      call putvar(db,gas_edep,'matsrc',4)

      call putvar(db,gas_mom(1,:),'sx',2)
      call putvar(db,gas_mom(2,:),'sy',2)
      call putvar(db,gas_mom(3,:),'sz',2)

      rc = DBClose(db)

      end subroutine output_silo
