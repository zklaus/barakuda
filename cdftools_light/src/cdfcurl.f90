PROGRAM cdfcurl
   !!---------------------------------------------------------------------------
   !!         ***  PROGRAM  cdfcurl  ***
   !!
   !!  **  Purpose: Compute the curl on F-points for given gridU gridV files and variables
   !!
   !! history :
   !!   Original :  J.M. Molines (May 2005)
   !!---------------------------------------------------------------------
   !!  $Rev: 256 $
   !!  $Date: 2009-07-21 17:49:27 +0200 (Tue, 21 Jul 2009) $
   !!  $Id: cdfcurl.f90 256 2009-07-21 15:49:27Z molines $
   !!--------------------------------------------------------------
   !! * Modules used
   USE cdfio

   !! * Local variables
   IMPLICIT NONE
   INTEGER :: ji,jj,jk, jt, ilev, istatus
   INTEGER :: npiglo, npjglo, npk, nt
   INTEGER :: narg, iargc, ncout, ierr
   INTEGER, DIMENSION(2) ::  ipk, id_varout         !lolo

   REAL(kind=8), DIMENSION(:,:), ALLOCATABLE  :: e2v, e1u, e1f, e2f
   REAL(kind=4), DIMENSION(:,:), ALLOCATABLE  :: un_t, vn_t, rotn, taum, tmask, fmask, zun, zvn
   REAL(KIND=4) ,DIMENSION(:), ALLOCATABLE    ::  tim

   CHARACTER(LEN=256) :: cfilu, cfilv, cvaru, cvarv, cdum, ctim
   CHARACTER(LEN=256) :: cf_mm='mesh_mask.nc', cfileout='curl.nc'
   TYPE (variable), DIMENSION(2) :: typvar         !: structure for attibutes

   !!
   narg = iargc()
   IF ( narg /= 5 ) THEN
      PRINT *,' USAGE : cdfcurl fileU fileV varU varV lev'
      PRINT *,'        lev is the level where the curl will be computed'
      PRINT *,'        Produce a cdf file curl.nc with socurl variable'
      PRINT *,'        Need mesh_hgr.nc'
      PRINT *,'                         '
      STOP
   ENDIF

   CALL getarg(1, cfilu)
   CALL getarg(2, cfilv)
   CALL getarg(3, cvaru)
   CALL getarg(4, cvarv)
   CALL getarg(5, cdum)
   READ(cdum,*) ilev

   ! define new variables for output ( must update att.txt)
   typvar(1)%name='socurl'
   typvar(1)%units='10^-6 s-1'
   typvar(1)%missing_value=-9999.
   typvar(1)%valid_min= -1000.
   typvar(1)%valid_max= 1000.
   typvar(1)%long_name='Relative_Vorticity (curl)'
   typvar(1)%short_name='socurl'
   typvar(1)%online_operation='N/A'
   typvar(1)%axis='TYX'

   typvar(2)%name='sotaum'
   typvar(2)%units='N/m^2'
   typvar(2)%missing_value=-9999.
   typvar(2)%valid_min= -1000.
   typvar(2)%valid_max= 1000.
   typvar(2)%long_name='Wind stress module (on T-point)'
   typvar(2)%short_name='sotaum'
   typvar(2)%online_operation='N/A'
   typvar(2)%axis='TYX'


   ipk(:) = 1  !  2D

   npiglo = getdim(cfilu,'x')
   npjglo = getdim(cfilu,'y')
   npk    = getdim(cfilu,'depth')

   ctim = 'none'
   nt    = getdim (cfilu,'time',cdtrue=ctim,kstatus=istatus) !LB

   PRINT *, 'npiglo =',npiglo
   PRINT *, 'npjglo =',npjglo
   PRINT *, 'npk    =',npk
   PRINT *, 'nt     =',nt  !PM
   PRINT *, 'ilev   =',ilev

   !test if lev exists
   IF ((npk==0) .AND. (ilev .GT. 0) ) THEN
      PRINT *, 'Problem : npk = 0 and lev > 0 STOP'
      STOP
   END IF

   ! if forcing field (PM)
   IF (ilev==0 .AND. npk==0) THEN
      !lforcing=.TRUE.
      npk = 1 ; ilev=1
      PRINT *, 'npk =0, assume 1'
   END IF

   IF (nt==0) THEN
      PRINT *, 'nt=0, assume 1'
      nt=1
   END IF
   !end (PM)

   ! check files and determines if the curl will be 2D of 3D

   ! create output fileset
   ncout =create(cfileout, cfilu, npiglo,npjglo,0)
   ierr= createvar(ncout ,typvar, 2, ipk, id_varout )
   ierr= putheadervar(ncout, cfilu, npiglo, npjglo, 0)


   ! Allocate the memory
   ALLOCATE ( e1u(npiglo,npjglo) , e1f(npiglo,npjglo) )
   ALLOCATE ( e2v(npiglo,npjglo) , e2f(npiglo,npjglo) )
   ALLOCATE ( un_t(npiglo,npjglo)  , vn_t(npiglo,npjglo)  )
   ALLOCATE ( zun(npiglo,npjglo) , zvn(npiglo,npjglo) )
   ALLOCATE ( rotn(npiglo,npjglo) , taum(npiglo,npjglo) )
   ALLOCATE ( tmask(npiglo,npjglo) , fmask(npiglo,npjglo) )
   ALLOCATE ( tim(nt) )

   e1u=  getvar(cf_mm, 'e1u', 1,npiglo,npjglo)
   e1f=  getvar(cf_mm, 'e1f', 1,npiglo,npjglo)
   e2v=  getvar(cf_mm, 'e2v', 1,npiglo,npjglo)
   e2f=  getvar(cf_mm, 'e2f', 1,npiglo,npjglo)
   fmask =  getvar(cf_mm, 'fmask', 1,npiglo,npjglo) ; !lolo
   tmask =  getvar(cf_mm, 'tmask', 1,npiglo,npjglo) ; !lolo


   tim=getvar1d(cfilu,trim(ctim),nt)
   ierr=putvar1d(ncout,tim,nt,'T')

   DO jt=1,nt

      PRINT *, ' * [cdfcurl] jt = ', jt

      ! if files are forcing fields
      jk = ilev
      zun(:,:) =  getvar(cfilu, cvaru, jk ,npiglo,npjglo, ktime=jt)
      zvn(:,:) =  getvar(cfilv, cvarv, jk ,npiglo,npjglo, ktime=jt)

      !! VU and VV on T-point:
      un_t = 0.
      DO ji=2, npiglo
         un_t(ji,:) = 0.5*(zun(ji-1,:) + zun(ji,:))*tmask(ji,:)
      END DO
      vn_t = 0.
      DO jj=2, npjglo
         vn_t(:,jj) = 0.5*(zvn(:,jj-1) + zvn(:,jj))*tmask(:,jj)
      END DO
      
      rotn(:,:) = 0.
      DO jj = 1, npjglo -1
         DO ji = 1, npiglo -1   ! vector opt.
            rotn(ji,jj) = 1.E6 * (  e2v(ji+1,jj  ) * zvn(ji+1,jj  ) - e2v(ji,jj) * zvn(ji,jj)    &
               &                   - e1u(ji  ,jj+1) * zun(ji  ,jj+1) + e1u(ji,jj) * zun(ji,jj)  ) &
               &           * fmask(ji,jj) / ( e1f(ji,jj) * e2f(ji,jj) )
         END DO
      END DO
      !
      rotn(npiglo,:) = rotn(2,:) ; ! periodicity lolo
      WHERE ( fmask == 0. ) rotn = -9999.

      ! write rotn on file at level k and at time jt
      ierr = putvar(ncout, id_varout(1), rotn, 1 ,npiglo, npjglo, jt)
      
      !! Wind stress module at T-point
      taum(:,:) = 0.
      taum(:,:) = SQRT( un_t(:,:)*un_t(:,:) + vn_t(:,:)*vn_t(:,:) )
      taum(1,:) = taum(npiglo-1,:) ; ! periodicity lolo
      WHERE ( tmask == 0. ) taum = -9999.
      ierr = putvar(ncout, id_varout(2), taum, 1 ,npiglo, npjglo, jt)

   END DO
   ierr = closeout(ncout)
   
END PROGRAM cdfcurl

