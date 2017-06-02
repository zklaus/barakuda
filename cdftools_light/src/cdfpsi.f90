PROGRAM cdfpsi
   !!-------------------------------------------------------------------
   !!               ***  PROGRAM cdfpsi  ***
   !!
   !!  **  Purpose  :  Compute Barotropic Stream Function
   !!                  PARTIAL STEPS
   !!
   !!  **  Method   :  Compute the 2D fields ztrpu, ztrpv
   !!                  as the integral on the vertical of u, v on their
   !!                  respective points.
   !!                  Then integrate from West to East   : ==> psiv
   !!                  (should be almost the same (if no error ))
   !!   Default (appropriate for global model): output psiu;
   !!                    normalizes the values setting psi (jpi,jpj) = 0
   !!   If option "V" is given as last argument, output psiv,
   !!                    normalizes values setting psi(jpi,1) = 0.
   !!                    This is appropriate for North Atlantic
   !!
   !! history ;
   !!  Original :  J.M. Molines (May 2005 )
   !!-------------------------------------------------------------------
   !!  $Rev: 256 $
   !!  $Date: 2009-07-21 17:49:27 +0200 (Tue, 21 Jul 2009) $
   !!  $Id: cdfpsi.f90 256 2009-07-21 15:49:27Z molines $
   !!--------------------------------------------------------------
   !! * Modules used
   USE cdfio
   USE io_ezcdf

   IMPLICIT NONE

   !INTEGER, PARAMETER :: nb_basins = 3
   !CHARACTER(len=3), DIMENSION(3), PARAMETER :: cb_lab = (/ 'atl','pac','ind' /)
   !CHARACTER(len=8), DIMENSION(3), PARAMETER :: cb_lgn = (/ 'Atlantic','Pacific ','Indian  ' /)
   INTEGER, PARAMETER :: nb_basins = 2
   CHARACTER(len=3), DIMENSION(nb_basins), PARAMETER :: cb_lab = (/ 'atl','pac' /)
   CHARACTER(len=8), DIMENSION(nb_basins), PARAMETER :: cb_lgn = (/ 'Atlantic','Pacific ' /)

   INTEGER   :: ji,jj,jk,jt                         !: dummy loop index
   INTEGER   :: ierr                                !: working integer
   INTEGER   :: narg, iargc                         !: command line
   INTEGER   :: npiglo,npjglo, npk, nt              !: size of the domain
   INTEGER   :: ncout
   INTEGER, DIMENSION(nb_basins) ::  ipk, id_varout         !
   TYPE(variable), DIMENSION(nb_basins)  :: typvar         !: structure for attributes

   REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE ::  zmask, e1v, e2u, glamf, gphif
   REAL(KIND=4), DIMENSION (:,:,:),   ALLOCATABLE ::  zmask_basin
   REAL(KIND=4) ,DIMENSION(:),      ALLOCATABLE ::  tim

   REAL(4), DIMENSION(:,:,:), ALLOCATABLE :: E3V_3D, V_3D


   REAL(KIND=8),   DIMENSION (:,:), ALLOCATABLE ::  ztrpv, psiv

   CHARACTER(LEN=256) :: cfileu ,cfilev, cfileoutnc='psi.nc', ctim
   CHARACTER(LEN=256) :: cf_mm='mesh_mask.nc'
   CHARACTER(LEN=256) :: cf_bm='new_maskglo.nc'
   CHARACTER(LEN=1)  :: coption


   CHARACTER(LEN=64) :: cv_u, cv_v

   INTEGER    :: istatus, jb

   INTEGER :: idf_0=0, idv_0=0,        &
      &    idf_u=0, idv_u=0, idf_v=0, idv_v=0

   ! constants

   !!  Read command line and output usage message if not compliant.
   narg= iargc()
   IF ( narg == 0 ) THEN
      PRINT *,' Usage : cdfpsi  Ufile Vfile nameU nameV (optional argument)'
      PRINT *,' Computes the barotropic stream function as the integral of the transport'
      PRINT *,' PARTIAL CELLS VERSION'
      PRINT *,' Files mesh_mask.nc and new_maskglo.nc must be in te current directory'
      PRINT *,' Output on psi.nc, variables sobarstf on f-points'
      PRINT *,' Default works well for a global ORCA grid. use V 3rdargument for North Atlantic'
      STOP
   ENDIF

   CALL getarg (1, cfileu  )
   CALL getarg (2, cfilev  )
   CALL getarg (3, cv_u)
   CALL getarg (4, cv_v)

   npiglo= getdim (cfileu,'x')
   npjglo= getdim (cfileu,'y')
   npk   = getdim (cfileu,'depth')

   ctim = 'none'
   nt    = getdim (cfileu,'time',cdtrue=ctim,kstatus=istatus) !LB

   !lolo
   DO jb = 1, nb_basins
      ! define new variables for output ( must update att.txt)
      typvar(jb)%name='psi_'//TRIM(cb_lab(jb))
      typvar(jb)%units='10^6 m^3/s'
      typvar(jb)%missing_value=-9999.
      typvar(jb)%valid_min= -300.e6
      typvar(jb)%valid_max= 300.e6
      typvar(jb)%long_name='Barotropic_Stream_Function_'//TRIM(cb_lgn(jb))
      typvar(jb)%short_name='psi_'//TRIM(cb_lab(jb))
      typvar(jb)%online_operation='N/A'
      typvar(jb)%axis='TYX'
      ipk(jb) = 1  !  2D ( X, Y , T )
   END DO


   PRINT *, 'npiglo=', npiglo
   PRINT *, 'npjglo=', npjglo
   PRINT *, 'npk   =', npk
   PRINT *, 'nt   =', nt !LB


   ! Allocate arrays
   ALLOCATE ( zmask(npiglo,npjglo), zmask_basin(npiglo,npjglo,nb_basins) )
   ALLOCATE ( e1v(npiglo,npjglo), e2u(npiglo,npjglo), E3V_3D(npiglo,npjglo,npk))
   ALLOCATE ( V_3D(npiglo,npjglo,npk) )
   ALLOCATE ( ztrpv(npiglo,npjglo), psiv(npiglo,npjglo))
   ALLOCATE ( glamf(npiglo,npjglo), gphif(npiglo,npjglo))
   ALLOCATE ( tim(nt) ) !LB

   glamf(:,:) = getvar(cf_mm, 'glamf',1,npiglo,npjglo)
   gphif(:,:) = getvar(cf_mm, 'gphif',1,npiglo,npjglo)

   ! create output fileset
   ncout =create(cfileoutnc, cfileu, npiglo,npjglo,1)
   ierr= createvar(ncout ,typvar,nb_basins, ipk,id_varout )
   ierr= putheadervar(ncout, cfileu,npiglo, npjglo,1,glamf, gphif)
   tim=getvar1d(cfileu,trim(ctim),nt) !LB
   ierr=putvar1d(ncout,tim,nt,'T') !LB


   e1v(:,:)   = getvar(cf_mm, 'e1v', 1,npiglo,npjglo)
   e2u(:,:)   = getvar(cf_mm, 'e2u', 1,npiglo,npjglo)
   zmask(:,:) = getvar(cf_mm, 'fmask', 1,npiglo,npjglo)



   ! get e3u, e3v  at all levels
   CALL GETVAR_3D(idf_0, idv_0, cf_mm, 'e3v_0', 0, 0, E3V_3D) ; idf_0 = 0. ; idv_0 = 0.

   ! get rid of the free-slip/no-slip condition
   WHERE ( zmask >= 2 ) zmask = 1



   DO jt=1,nt

      PRINT *, ' * [cdfpsi] jt = ', jt

      CALL GETVAR_3D(idf_v, idv_v, cfilev, cv_v, nt, jt, V_3D)

      !! Loop over basins:
      DO jb = 1, nb_basins
         
         PRINT *, ''
         
         IF ( jt == 1 ) THEN
            PRINT *, ' Getting mask for basin ', TRIM(cb_lgn(jb)),' => ', TRIM(cb_lab(jb))
            zmask_basin(:,:,jb) = getvar(cf_bm, 'tmask'//TRIM(cb_lab(jb)), 1,npiglo,npjglo) * zmask(:,:)
         END IF
         
         ztrpv(:,:)= 0.d0
         !psiv(npiglo,:)= 0.d0
         psiv(:,:) = 0.d0

         DO jk = 1,npk
            ztrpv(:,:) = ztrpv(:,:) + V_3D(:,:,jk)*e1v(:,:)*E3V_3D(:,:,jk)  ! meridional transport of each grid cell
         END DO  ! loop to next level


         ! integrate zonally from east to west            
         DO ji = npiglo-1, 1, -1
            psiv(ji,:) = psiv(ji+1,:) - ztrpv(ji,:)  ! psi at f point
         END DO
         psiv(:,:) = 1.E-6 * psiv(:,:)*zmask_basin(:,:,jb)
         WHERE ( zmask_basin(:,:,jb) == 0. ) psiv(:,:) = -9999.
         WHERE ( psiv(:,:) < -10000. ) psiv(:,:) = 0.
            
         ierr = putvar(ncout, id_varout(jb), SNGL(psiv), 1, npiglo, npjglo, jt) !LB
            
      END DO  ! DO jb = 1, nb_basins

   END DO  ! DO jt=1,nt
   

   istatus = closeout (ncout)

END PROGRAM cdfpsi
