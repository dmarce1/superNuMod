
      subroutine putvar( db, var, name0, namelen0, opt, ndims )

      use gasmod
      use gridmod
      implicit none

      include '/home/dmarce1/include/silo.inc'

      integer, intent(in) :: db, namelen0, opt, ndims
      character(len=namelen0) :: name0
      real*8, intent(in)  :: var(gas_ncell)
      real*8, allocatable :: v3d(:,:,:)
      character(len=len_trim(adjustl(name0))) :: name
      integer :: namelen


      integer :: rc, i, j, k, l
      integer, dimension(3) :: dims
      integer :: ctt = DB_ZONECENT
      integer :: datat = DB_DOUBLE
      integer :: dum

      name = trim(adjustl(name0))
      namelen = len(name)

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
c          DBPUTQV1_FC (int *dbid, FCD_DB name,
c           int *lname, FCD_DB meshname, int *lmeshname,
c           void const *var, int *dims, int *ndims, void const *mixvar, int *mixlen,
c           int *datatype, int *centering, int *optlist_id, int *status)c

      rc=DBPutQV1(db,name,
     & namelen,'mesh',4,
     & v3d,dims,ndims,0,0,datat,ctt,opt,dum)

      deallocate(v3d)
      end subroutine putvar







      subroutine output_silo(tm,ocnt)

      use gridmod
      use gasmod
      use elemdatamod
      implicit none

      include '/home/dmarce1/include/silo.inc'

      real*8, intent(in) :: tm

      character(len=1024) :: oname
      integer, intent(in) :: ocnt

      integer :: db, rc, i, opt
      real*8 :: fac, tmp

      integer :: dims(3)
      integer :: datat = DB_DOUBLE
      integer :: crdt = DB_COLLINEAR
      integer :: ndims
      real*8 :: crdx(grd_nx+1)
      real*8 :: crdy(grd_ny+1)
      real*8 :: crdz(grd_nz+1)
      logical :: issphere

      integer :: dum

      if( grd_igeom .eq. 11 ) then
       ndims = 1
      else
       ndims = 3
      endif

c DBMKOPTLIST_FC (int *maxopts, int *optlist_id)
      rc = DBMkOptlist(1,opt)

      if( grd_isvelocity ) then
       fac = tm
      else
       fac = 1.0d0
      endif

      dims(1) = grd_nx+1
      dims(2) = grd_ny+1
      dims(3) = grd_nz+1

      if( grd_igeom .eq. 11 .or. grd_igeom .eq. 1 ) then
       issphere = .true.
      else
       issphere = .false.
      endif

      do i=1,grd_nx+1
       crdx(i) = grd_xarr(i) * fac
      enddo
      do i=1,grd_ny+1
       if( issphere ) then
        crdy(i) = 3.14159265359d0 - acos(grd_yarr(i))
       else
        crdy(i) = grd_yarr(i) * fac
       endif
      enddo
      do i=1,grd_nz+1
       if( issphere ) then
        crdz(i) = 2.0d0 * 3.14159265359d0 - grd_zarr(i)
       else
        crdz(i) = grd_zarr(i) * fac
       endif
      enddo



100   format(A2,I0,A5)
      write (oname, 100) 'X.', ocnt, '.silo'
      print *, trim(oname), tm
      rc=DBCreate(trim(oname),len_trim(oname),
     & DB_CLOBBER, DB_LOCAL,'',0,DB_PDB,db)

c DBPUTQM_FC (int *dbid, FCD_DB name, int *lname, FCD_DB xname, int *lxname,
c          FCD_DB yname, int *lyname, FCD_DB zname, int *lzname, void const *x,
c          void const *y, void const *z, int *dims, int *ndims, int *datatype,
c          int *coordtype, int *optlist_id, int *status)

      rc = DBPutQM(db,'mesh',4,'x',1,
     & 'y',1,'z',1,crdx,
     & crdy,crdz,dims,ndims,datat,crdt,opt,dum)
      do i=1,elem_neldata
       call putvar(db,gas_natom1fr(i,:),elem_data(i)%sym,2,opt,ndims)
      enddo
      call putvar(db,gas_natom1fr(-1,:),'ni56',4,opt,ndims)
      call putvar(db,gas_natom1fr(-2,:),'co56',4,opt,ndims)
      call putvar(db,gas_natom1fr(-3,:),'fe52',4,opt,ndims)
      call putvar(db,gas_natom1fr(-4,:),'mn52',4,opt,ndims)
      call putvar(db,gas_natom1fr(-5,:),'cr48',4,opt,ndims)
      call putvar(db,gas_natom1fr(-6,:),'v48',3,opt,ndims)
      call putvar(db,gas_temp,'temp',4,opt,ndims)
      call putvar(db,gas_eraddens,'eraddens',8,opt,ndims)
      call putvar(db,gas_ur,'ur',2,opt,ndims)
      call putvar(db,gas_rho,'rho',3,opt,ndims)
      call putvar(db,gas_bcoef,'bcoef',5,opt,ndims)
      call putvar(db,gas_decaygamma,'decaygamma',10,opt,ndims)
      call putvar(db,gas_decaybeta,'decaybeta',9,opt,ndims)
      call putvar(db,gas_ye,'Ye',2,opt,ndims)
      call putvar(db,gas_natom,'natom',5,opt,ndims)
      call putvar(db,gas_nelec,'nelec',5,opt,ndims)
      call putvar(db,gas_matsrc,'matsrc',6,opt,ndims)
      call putvar(db,gas_edep,'matsrc',4,opt,ndims)

      call putvar(db,gas_mom(1,:),'sx',2,opt,ndims)
      call putvar(db,gas_mom(2,:),'sy',2,opt,ndims)
      call putvar(db,gas_mom(3,:),'sz',2,opt,ndims)

      rc = DBClose(db)
      rc = DBFreeOptlist(opt)

      end subroutine output_silo
