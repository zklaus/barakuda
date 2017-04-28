PROGRAM cdficeflux

   !!---------------------------------------------------------------------
   !!               ***  PROGRAM cdficeflux  ***
   !!
   !!  ** Method  : The begining and end point of the section are given in
   !!               term of F-points index. A broken line joining successive
   !!               F-points is defined between the begining and end point
   !!               of the section. Therefore each segment between F-points
   !!               is either a zonal or meridional segment corresponding to
   !!               V or U velocity component. Doing so, the volume conservation
   !!               is ensured as velocities are not interpolated, and stay
   !!               on the native model grid.
   !!                 The section name and the begin/end point of a section are
   !!               read from standard input, till 'EOF' is given as section
   !!               name. This make possible to give a bunch of sections in
   !!               an ASCII files and use the < redirection.
   !!            SIGN CONVENTION : The transport is positive when the flow cross
   !!               the section to the right, negative otherwise. This depends
   !!               on the sense the section is described.  With this convention
   !!               The algebric sum of transports accross sections forming a
   !!               closed area is 0.
   !!            OPTIONS :
   !!               -time   : specify the time frame to be used
   !!            REQUIREMENT :
   !!               mesh-mask file are required in the current directory.
   !!
   !!
   !! History : 2.1  : 01/2005  : J.M. Molines : Original code
   !!           2.1  : 07/2009  : R. Dussin : add cdf output
   !!           2.1  : 01/2010  : M.A. Balmaseda : Change integration signs
   !!                             so that the transport across a segment is
   !!                             independent of the chosen trajectory.
   !!           3.0  : 04/2011  : J.M. Molines : Doctor norm + Lic.
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   routines      : description
   !!  interm_pt  : choose intermediate points on a broken line.
   !!----------------------------------------------------------------------

   USE netcdf

   USE cdfio

   IMPLICIT NONE

   REAL(8) :: &
      &        ref_sali0 = 34.8 ! reference salinity (in PSU)
   CHARACTER(len=128) :: cref2='34.8'


   INTEGER(KIND=4)                             :: jseg   ! dummy loop index
   INTEGER(KIND=4)                             :: ji, jj, jk     ! dummy loop index
   INTEGER(KIND=4), DIMENSION(:),  ALLOCATABLE :: ipk, id_varout ! Netcdf output
   INTEGER(KIND=4)                             :: ipos           ! working integer (position of ' ' in strings)
   INTEGER(KIND=4)                             :: nvarout=12     ! number of values to write in cdf output
   INTEGER(KIND=4)                             :: ivtrp          ! var index of volume transport (barotrope)
   INTEGER(KIND=4)                             :: iptrp          ! var index of volume transport (barotrope)
   INTEGER(KIND=4)                             :: imtrp          ! var index of volume transport (barotrope)
   INTEGER(KIND=4)                             :: ivtrpcl        ! var index of volume transport (p. class)
   INTEGER(KIND=4)                             :: iptrpcl        ! var index of volume transport (p. class)
   INTEGER(KIND=4)                             :: imtrpcl        ! var index of volume transport (p. class)
   INTEGER(KIND=4)                             :: ilonmin        ! var index of starting section longitude
   INTEGER(KIND=4)                             :: ilonmax        ! var index of ending section longitude
   INTEGER(KIND=4)                             :: ilatmin        ! var index of starting section latitude
   INTEGER(KIND=4)                             :: ilatmax        ! var index of ending section latitude
   INTEGER(KIND=4)                             :: itop           ! var index of top depth class
   INTEGER(KIND=4)                             :: ibot           ! var index of bottom depth class
   INTEGER(KIND=4)                             :: ikx=1, iky=1   ! dims of netcdf output file
   INTEGER(KIND=4)                             :: numin   = 10   ! logical unit for input section file (overall) !LB
   INTEGER(KIND=4)                             :: numout  = 11   ! logical unit for output file (overall)
   INTEGER(KIND=4)                             :: narg, iargc    ! command line
   INTEGER(KIND=4)                             :: ijarg, nxtarg  !  "       "
   INTEGER(KIND=4)                             :: npiglo, npjglo ! size of the domain
   INTEGER(KIND=4)                             :: npk, npt       ! size of the domain
   INTEGER(KIND=4)                             :: iimin, iimax   ! i-limit of the section
   INTEGER(KIND=4)                             :: ijmin, ijmax   ! j-limit of the section
   INTEGER(KIND=4)                             :: ivar, jt    ! working integer
   INTEGER(KIND=4)                             :: ii, ij, ik     ! working integer
   INTEGER(KIND=4), PARAMETER                  :: jpseg=10000    ! used for broken line algorithm
   INTEGER(KIND=4)                             :: ii0, ij0       !  "        "             "
   INTEGER(KIND=4)                             :: ii1, ij1       !  "        "             "
   INTEGER(KIND=4)                             :: iitmp, ijtmp   !  "        "             "
   INTEGER(KIND=4)                             :: np, nn         ! segment counters,
   INTEGER(KIND=4)                             :: iist, ijst     ! local point offset for velocity
   INTEGER(KIND=4)                             :: norm_u, norm_v ! normalization factor (sign of normal transport)
   INTEGER(KIND=4)                             :: idirx, idiry   ! sense of description of the section

   REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: e1t, e2t       ! horizontal metric
   REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: glamf          ! longitudes of F points
   REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: gphif          ! latitudes of F points
   !LB:
   REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: zu_i, zv_i, zc_i, zh_i ! velocities, concentration and thickness of sea-ice
   !LB.
   REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: rdum           ! dummy (1x1) array for ncdf output
   REAL(KIND=4), DIMENSION(:),     ALLOCATABLE :: tim            ! time counter
   REAL(KIND=4), DIMENSION(2)                  :: gla, gphi      ! lon/lat of the begining/end of section (f point)
   REAL(KIND=4), DIMENSION(jpseg)              :: rxx, ryy       ! working variables
   REAL(KIND=4)                                :: dvoltrpsum, ryear     ! LB
   REAL(KIND=4)                                :: rxi0, ryj0     ! working variables
   REAL(KIND=4)                                :: rxi1, ryj1     ! working variables
   REAL(KIND=4)                                :: ai, bi         ! equation of line (y=ai.x +bi)
   REAL(KIND=4)                                :: aj, bj         ! equation of line (x=aj.y +bj
   REAL(KIND=4)                                :: rd, rd1, rd2   ! distance between point, between vertical layers
   REAL(KIND=4)                                :: udum, vdum     ! dummy velocity components for tests

   ! at every model point
   REAL(KIND=8), DIMENSION(:,:),   ALLOCATABLE :: dwku, dwkv, xtmp  ! volume transport at each cell boundary

   ! for a given section

   REAL(KIND=8), DIMENSION(jpseg)              :: dvoltrp        ! volume transport across each segment of a section
   REAL(KIND=8), DIMENSION(jpseg)              :: dvoltrpp       ! volume transport across each segment of a section
   REAL(KIND=8), DIMENSION(jpseg)              :: dvoltrpm       ! volume transport across each segment of a section
   REAL(KIND=8)                                :: dvolalleg      ! over all leg sum of volume transport
   REAL(KIND=8)                                :: dvolallegp     ! over all leg sum of volume transport +
   REAL(KIND=8)                                :: dvolallegm     ! over all leg sum of volume transport -
   REAL(KIND=8)                                :: dheatalleg     ! over all leg sum of heat transport
   REAL(KIND=8)                                :: dsaltalleg     ! over all leg sum of salt transport

   COMPLEX, DIMENSION(jpseg)                   :: yypt           ! array of points coordinates in a section
   COMPLEX                                     :: yypti          ! working point

   TYPE(variable), DIMENSION(:),   ALLOCATABLE :: stypvar        ! structure of output

   CHARACTER(LEN=256) , PARAMETER              :: cn_fmm = 'mesh_mask.nc'

   CHARACTER(LEN=256)                          :: cf_sections='transport_ice.dat'  ! input file containing sections !LB
   CHARACTER(LEN=256)                          :: cf_icefil, cdum      ! seaice file   (in) !LB
   CHARACTER(LEN=4)                            :: cyear
   CHARACTER(LEN=256)                          :: cv_iceF, cv_iceU, cv_iceV, cv_iceH
   CHARACTER(LEN=256)                          :: cf_out         ! output file name (ASCII)
   CHARACTER(LEN=256)                          :: csection            ! section names
   CHARACTER(LEN=256)                          :: cvarname            ! variable names (root)
   CHARACTER(LEN=256)                          :: clongname           ! variable longname (root)
   CHARACTER(LEN=512)                          :: cglobal             ! global attribute
   CHARACTER(LEN=256)                          :: cldum               ! dummy char variable
   CHARACTER(LEN=256)                          :: cline               ! dummy char variable
   CHARACTER(LEN=256), DIMENSION(3)            :: cldumt              ! dummy char variable

   CHARACTER(LEN=20) :: ce1, ce2, cvmask


   LOGICAL                 :: ltest   = .FALSE.   ! flag for test case
   LOGICAL                 :: lchk    = .FALSE.   ! flag for missing files

   CHARACTER(LEN=64)  :: cv_dum
   CHARACTER(LEN=512) :: cd_out
   LOGICAL                 :: lfncout = .FALSE.
   INTEGER :: ierr, jt_pos, idf_out, idd_t, idv_time, id_volu

   REAL(4), DIMENSION(:), ALLOCATABLE :: X_trsp   ! lolo

   !!----------------------------------------------------------------------
   !CALL ReadCdfNames()

   narg= iargc()
   ! Print usage if no argument
   IF ( narg < 7 ) THEN
      PRINT *,' usage : cdficeflux <ICEMOD-file> <ice_conc> <ice_veloc_u> <ice_veloc_v> <ice_volume> <year> <DIROUT>'
      PRINT *,'      '
      PRINT *,'    PURPOSE :'
      PRINT *,'      Compute the transports of solid freshwater accross a section.'
      PRINT *,'      The name of the section and the imin, imax, jmin, jmax for the section '
      PRINT *,'      is read from the standard input. To finish the program use the key name'
      PRINT *,'      ''EOF'' for the section name.'
      PRINT *,'      '
      PRINT *,'     ARGUMENTS :'
      PRINT *,'      [ICEMOD-file] : netcdf file with the zonal velocity component.'
      PRINT *,'      [ice_conc]    : name of ice concentration (0<= ... <=1)'
      PRINT *,'      [ice_veloc_u] : name of ice velocity along i (m/s)'
      PRINT *,'      [ice_veloc_v] : name of ice velocity along j (m/s)'
      PRINT *,'      [ice_volume]  : name of ice volume/thickness (m) '
      PRINT *,'      [year]        : Year'
      PRINT *,'      [DIROUT]      : directory where to write netcdf file'
      PRINT *,'      '
      PRINT *,'      '
      PRINT *,'     REQUIRED FILES :'
      PRINT *,'      File ',TRIM(cn_fmm),' must be in the current directory.'
      PRINT *,'      '
      PRINT *,'     OUTPUT : '
      PRINT *,'      - Netcdf files for each section. name of the file is buildt'
      PRINT *,'          from section name.'
      PRINT *,'      '
      PRINT *,'     SEE ALSO :'
      PRINT *,'       cdfsigtrp'
      PRINT *,'      '
      STOP
   ENDIF

   jt  = 1
   ijarg  = 1

   CALL getarg (1, cf_icefil)
   CALL getarg (2, cv_iceF)
   CALL getarg (3, cv_iceU)
   CALL getarg (4, cv_iceV)
   CALL getarg (5, cv_iceH)
   CALL getarg (6, cyear)    ; READ(cyear,*) ryear
   CALL getarg (7, cd_out)

   PRINT *, '   *** Ice concentration => ', TRIM(cv_iceF)
   PRINT *, '   *** Ice U-velocity    => ', TRIM(cv_iceU)
   PRINT *, '   *** Ice V-velocity    => ', TRIM(cv_iceV)
   PRINT *, '   *** Ice thickness     => ', TRIM(cv_iceH)
   PRINT *, ''

   ce1='e1t'
   ce2='e2t'
   cvmask='tmask'


   ! adjust the number of output variables according to options
   nvarout = 1


   npiglo = getdim (cf_icefil,'x')
   npjglo = getdim (cf_icefil,'y')
   npt    = getdim (cf_icefil,'time')

   PRINT *, 'npiglo =', npiglo
   PRINT *, 'npjglo =', npjglo
   PRINT *, 'npt    =', npt

   ! define new variables for output
   ALLOCATE ( stypvar(nvarout), ipk(nvarout), id_varout(nvarout) )
   ALLOCATE ( rdum(1,1) )

   rdum(:,:)=0.e0

   ! Allocate arrays
   ALLOCATE ( zu_i(npiglo,npjglo),  zv_i(npiglo,npjglo), zc_i(npiglo,npjglo),  zh_i(npiglo,npjglo) )
   ALLOCATE ( dwku(npiglo,npjglo), dwkv(npiglo,npjglo), xtmp(npiglo,npjglo) )


   ALLOCATE ( e1t(npiglo,npjglo) )
   ALLOCATE ( e2t(npiglo,npjglo) )
   !
   ALLOCATE ( gphif(npiglo,npjglo) )
   ALLOCATE ( glamf(npiglo,npjglo) )
   ALLOCATE ( tim(npt)                       )
   !
   ALLOCATE( X_trsp(npt) )


   ! read metrics and grid position
   e1t(:,:)   = getvar(cn_fmm, ce1, 1, npiglo, npjglo)
   e2t(:,:)   = getvar(cn_fmm, ce2, 1, npiglo, npjglo)

   glamf(:,:) = getvar(cn_fmm, 'glamf', 1,npiglo, npjglo)
   gphif(:,:) = getvar(cn_fmm, 'gphif', 1,npiglo, npjglo)


   DO WHILE ( 1 == 1 )
      OPEN(numin, FILE=cf_sections, status='old')
      READ(numin,'(a)') cline
      ii = 0
      cldumt(:) = 'none'
      ipos = index(cline,' ')
      DO WHILE ( ipos > 1 )
         ii = ii + 1
         cldumt(ii) = cline(1:ipos - 1 )
         cline = TRIM ( cline(ipos+1:) )
         ipos  = index( cline,' ' )
         IF ( ii >= 3 ) EXIT
      END DO
      csection = TRIM(cldumt(1) )
      cvarname = TRIM(cldumt(2) )
      clongname = TRIM(cldumt(3) )

      IF (TRIM(csection) == 'EOF' ) THEN
         !PRINT *, 'LB! exiting'
         CLOSE(numin)
         GOTO 1111 ; ! LB
         !EXIT  ! infinite DO loop
      ENDIF
      PRINT *, ''
      PRINT *, '   cdficeflux: doing section '//trim(csection)


      DO jt = 1, npt   !LB

         PRINT *, ' * [cdficeflux] jt = ', jt

         ! compute the transports at each grid cell

         ! Get velocities, concentration and thickness for sea-ice:
         IF ( ltest ) THEN
            zu_i (:,:) = udum ; zv_i (:,:) = vdum
         ELSE
            zu_i (:,:) = getvar(cf_icefil, cv_iceU, 1, npiglo, npjglo, ktime=jt) ; ! LB
            zv_i (:,:) = getvar(cf_icefil, cv_iceV, 1, npiglo, npjglo, ktime=jt) ; ! LB
         ENDIF
         zc_i (:,:) = getvar(cf_icefil, cv_iceH, 1, npiglo, npjglo, ktime=jt) ; ! LB
         zh_i (:,:) = getvar(cf_icefil, cv_iceF, 1, npiglo, npjglo, ktime=jt) ; ! LB
         !LB2JM: add all the ice variables into modcdfnames?


         !LB:
         !! Volume of ice transport at each grid point:
         !!  => concentration * thickness * u * dy
         !!  => concentration * thickness * v * dx
         xtmp(:,:) = zh_i(:,:)*zc_i(:,:)
         dwku(:,:) = xtmp(:,:)*zu_i(:,:)*e2t(:,:)*1.d0
         dwkv(:,:) = xtmp(:,:)*zv_i(:,:)*e1t(:,:)*1.d0


         IF ( jt == 1 ) THEN   !LB

            READ(numin,*) iimin, iimax, ijmin, ijmax
            !! Find the broken line between P1 (iimin,ijmin) and P2 (iimax, ijmax)
            ! ... Initialization
            ii0  = iimin ; ij0  = ijmin ; ii1  = iimax ;  ij1 = ijmax
            rxi0 = ii0   ; ryj0 = ij0   ; rxi1 = ii1   ; ryj1 = ij1
            
            !By defining the direction of the integration as
            idirx = SIGN(1,ii1-ii0) !positive to the east or if ii1=ii0
            idiry = SIGN(1,ij1-ij0) !positive to the north or if ij1=ij0
            
            !Then dS=(e2t*idiry,-e1t*idirx)
            !This will produce the following sign convention:
            !    West-to-est line (dx>0, dy=0)=> -My*dx (-ve for a northward flow)
            !    South-to-north   (dy>0, dx=0)=>  Mx*dy (+ve for an eastward flow)
            norm_u =  idiry
            norm_v = -idirx

            ! .. Compute equation:  ryj = aj rxi + bj [valid in the (i,j) plane]
            IF ( (rxi1 -rxi0) /=  0 ) THEN
               aj = (ryj1 - ryj0 ) / (rxi1 -rxi0)
               bj = ryj0 - aj * rxi0
            ELSE
               aj = 10000.  ! flag value
               bj = 0.
            END IF

            ! .. Compute equation:  rxi = ai ryj + bi [valid in the (i,j) plane]
            IF ( (ryj1 -ryj0) /=  0 ) THEN
               ai = (rxi1 - rxi0 ) / ( ryj1 -ryj0 )
               bi = rxi0 - ai * ryj0
            ELSE
               ai = 10000. ! flag value
               bi = 0.
            END IF

            ! ..  Compute the integer pathway: a succession of F points
            np=0
            ! .. Chose the strait line with the smallest slope
            IF (ABS(aj) <=  1 ) THEN
               ! ... Here, the best line is y(x)
               ! ... If ii1 < ii0 swap points [ always describe section from left to right ]
               IF (ii1 <  ii0 ) THEN
                  iitmp = ii0   ; ijtmp = ij0
                  ii0   = ii1   ; ij0   = ij1
                  ii1   = iitmp ; ij1   = ijtmp
               END IF

               ! iist,ijst is the grid offset to pass from F point to either U/V point
               IF ( ij1 >= ij0 ) THEN     ! line heading NE
                  iist = 1 ; ijst = 1
               ELSE                       ! line heading SE
                  iist = 1 ; ijst = 0
               END IF

               ! ... compute the nearest ji point on the line crossing at ji
               DO ji=ii0, ii1
                  np=np+1
                  IF (np > jpseg) STOP 'np > jpseg !'
                  ij=NINT(aj*ji + bj )
                  yypt(np) = CMPLX(ji,ij)
               END DO
            ELSE
               ! ... Here, the best line is x(y)
               ! ... If ij1 < ij0 swap points [ always describe section from bottom to top ]
               IF (ij1 <  ij0 ) THEN
                  iitmp = ii0   ; ijtmp = ij0
                  ii0   = ii1   ; ij0   = ij1
                  ii1   = iitmp ; ij1   = ijtmp
               END IF

               ! iist,ijst is the grid offset to pass from F point to either U/V point
               IF ( ii1 >=  ii0 ) THEN
                  iist = 1 ; ijst = 1
               ELSE
                  iist = 0 ; ijst = 1
               END IF

               ! ... compute the nearest ji point on the line crossing at jj
               DO jj=ij0,ij1
                  np=np+1
                  IF (np > jpseg) STOP 'np > jpseg !'
                  ii=NINT(ai*jj + bi)
                  yypt(np) = CMPLX(ii,jj)
               END DO
            END IF

            !!
            !! Look for intermediate points to be added.
            !  ..  The final positions are saved in rxx,ryy
            rxx(1) = REAL(yypt(1))
            ryy(1) = IMAG(yypt(1))
            nn     = 1

            DO jk=2,np
               ! .. distance between 2 neighbour points
               rd=ABS(yypt(jk)-yypt(jk-1))
               ! .. intermediate points required if rd > 1
               IF ( rd > 1 ) THEN
                  CALL interm_pt(yypt, jk, ai, bi, aj, bj, yypti)
                  nn=nn+1
                  IF (nn > jpseg) STOP 'nn>jpseg !'
                  rxx(nn) = REAL(yypti)
                  ryy(nn) = IMAG(yypti)
               END IF
               nn=nn+1
               IF (nn > jpseg) STOP 'nn>jpseg !'
               rxx(nn) = REAL(yypt(jk))
               ryy(nn) = IMAG(yypt(jk))
            END DO
            ! record longitude and latitude of initial en endind point of the section
            gla (1) = glamf( INT(rxx(1)),  INT(ryy(1))  )
            gphi(1) = gphif( INT(rxx(1)),  INT(ryy(1))  )
            gla (2) = glamf( INT(rxx(nn)), INT(ryy(nn)) )
            gphi(2) = gphif( INT(rxx(nn)), INT(ryy(nn)) )

            ! Now extract the transport through a section
            ! ... Check whether we need a u velocity or a v velocity
            !   Think that the points are f-points and delimit either a U segment
            !   or a V segment (iist and ijst are set in order to look for the correct
            !   velocity point on the C-grid

         END IF ! jt == 1


         ! Now really calculate the flux for time jt!

         dvoltrpsum   = 0.d0

         ! segment jseg is a line between (rxx(jseg),ryy(jseg))  and (rxx(jseg+1),ryy(jseg+1))
         DO jseg = 1, nn-1

            ii0=rxx(jseg)
            ij0=ryy(jseg)

            IF ( rxx(jseg) ==  rxx(jseg+1) ) THEN    ! meridional segment, use U velocity
               dvoltrp(jseg) = dwku(ii0,ij0+ijst)*norm_u

            ELSE IF ( ryy(jseg) == ryy(jseg+1) ) THEN ! zonal segment, use V velocity
               dvoltrp(jseg) = dwkv(ii0+iist,ij0)*norm_v

            ELSE
               PRINT *,' ERROR :',  rxx(jseg),ryy(jseg),rxx(jseg+1),ryy(jseg+1) ! likely to never happen !
            END IF

            dvoltrpsum = dvoltrpsum + dvoltrp(jseg)

         END DO   ! next segment


         !! Correcting value, due to expension of ice and salinity of ice:
         dvoltrpsum = 0.92*dvoltrpsum   ; ! 1m^3 of ice makes 0.92 m^3 of water

         dvoltrpsum = dvoltrpsum * ( ref_sali0 - 4. )/ref_sali0  ; ! salinity of sea-ice taken as 4 PSU
         !!                                              ! referencesalinity for Arctic: 34.8 PSU

         X_trsp(jt) = dvoltrpsum/1.e6  ! (Sv)


      END DO ! loop on time steps ( jt = 1 => npt


      PRINT *, ''

      !! Time to create or append X_trsp into the netcdf file for current section
      !! -----------------------------------------------------

      id_volu = 0

      WRITE(cf_out, '(a,"/transport_ice_sect_",a,".nc")') trim(cd_out), trim(csection)

      INQUIRE( FILE=cf_out, EXIST=lfncout )

      IF ( .NOT. lfncout ) THEN

         !! Creating file
         PRINT *, ' Creating file '//TRIM(cf_out)//' !!!'
         ierr = NF90_CREATE(cf_out, NF90_CLOBBER, idf_out)
         ierr = NF90_DEF_DIM(idf_out, 'time', NF90_UNLIMITED, idd_t)
         ierr = NF90_DEF_VAR(idf_out, 'time', NF90_DOUBLE,    idd_t, idv_time)


         ierr = NF90_DEF_VAR(idf_out, 'trsp_volu', NF90_FLOAT, (/idd_t/), id_volu)
         ierr = NF90_PUT_ATT(idf_out, id_volu, 'long_name', 'Transport of (solid) freshwater due to sea-ice drift')
         ierr = NF90_PUT_ATT(idf_out, id_volu, 'units', 'Sv')
         ierr = NF90_PUT_ATT(idf_out, id_volu, 'Sref', cref2)
         ierr = NF90_PUT_ATT(idf_out, NF90_GLOBAL, 'About', &
            &   'Created by BaraKuda (cdficeflux.x) => https://github.com/brodeau/barakuda')

         ierr = NF90_ENDDEF(idf_out)
         jt_pos = 0


      ELSE

         !! Opening already existing file
         ierr = NF90_OPEN  (cf_out, NF90_WRITE,   idf_out)

         !! Need IDs of variables to append... NF90_INQ_VARID
         ierr = NF90_INQ_VARID(idf_out, 'trsp_volu', id_volu)

         ! Get ID of unlimited dimension
         ierr = NF90_INQUIRE(idf_out, unlimitedDimId = idv_time)

         ! Need to know jt_pos, record number of the last time record writen in the file
         ierr = NF90_INQUIRE_DIMENSION(idf_out, idv_time, name=cv_dum, len=jt_pos)

      END IF

      WRITE(*,'("Going to write record ",i4.4," to ",i4.4," into ",a)') jt_pos+1, jt_pos+npt, trim(cf_out)

      DO jt = 1, npt

         ! Writing record jt for time vector and 1d fields:
         ierr = NF90_PUT_VAR( idf_out, idv_time, (/ryear+1./12.*(REAL(jt)-1.+0.5)/), start=(/jt_pos+jt/), count=(/1/) )

         !! Default variable is the only one present (index = 1) :
         ierr = NF90_PUT_VAR(idf_out, id_volu, (/ X_trsp(jt) /), start=(/jt_pos+jt/), count=(/1/))

      END DO

      ierr = NF90_CLOSE(idf_out)



   END DO ! infinite loop : gets out when input is EOF




1111 CONTINUE








CONTAINS



   SUBROUTINE interm_pt (ydpt, kk, pai, pbi, paj, pbj, ydpti)
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE nterm_pt  ***
      !!
      !! ** Purpose : Find the best intermediate points on a pathway.
      !!
      !! ** Method  : ydpt : complex vector of the positions of the nearest points
      !!               kk  : current working index
      !!          pai, pbi : slope and original ordinate of x(y)
      !!          paj, pbj : slope and original ordinate of y(x)
      !!             ydpti : Complex holding the position of intermediate point
      !!
      !! ** Reference : 19/07/1999 : J.M. Molines in Clipper
      !!----------------------------------------------------------------------
      COMPLEX, DIMENSION(:), INTENT(in ) :: ydpt
      COMPLEX,               INTENT(out) :: ydpti
      REAL(KIND=4),          INTENT(in ) :: pai, pbi, paj, pbj
      INTEGER(KIND=4),       INTENT(in ) :: kk
      ! ... local
      COMPLEX                            :: ylptmp1, ylptmp2
      REAL(KIND=4)                       :: za0, zb0
      REAL(KIND=4)                       :: za1, zb1
      REAL(KIND=4)                       :: zd1, zd2
      REAL(KIND=4)                       :: zxm, zym
      REAL(KIND=4)                       :: zxp, zyp
      !!----------------------------------------------------------------------
      ! ... Determines whether we use y(x) or x(y):
      IF (ABS(paj) <=  1) THEN
         ! .....  use y(x)
         ! ... possible intermediate points:
         ylptmp1=ydpt(kk-1)+(1.,0.)                 ! M1
         ylptmp2=ydpt(kk-1)+CMPLX(0.,SIGN(1.,paj))  ! M2
         !
         ! ... M1 is the candidate point:
         zxm=REAL(ylptmp1)
         zym=IMAG(ylptmp1)
         za0=paj
         zb0=pbj
         !
         za1=-1./za0
         zb1=zym - za1*zxm
         ! ... P1 is the projection of M1 on the strait line
         zxp=-(zb1-zb0)/(za1-za0)
         zyp=za0*zxp + zb0
         ! ... zd1 is the distance M1P1
         zd1=(zxm-zxp)*(zxm-zxp) + (zym-zyp)*(zym-zyp)
         !
         ! ... M2 is the candidate point:
         zxm=REAL(ylptmp2)
         zym=IMAG(ylptmp2)
         za1=-1./za0
         zb1=zym - za1*zxm
         ! ... P2 is the projection of M2 on the strait line
         zxp=-(zb1-zb0)/(za1-za0)
         zyp=za0*zxp + zb0
         ! ... zd2 is the distance M2P2
         zd2=(zxm-zxp)*(zxm-zxp) + (zym-zyp)*(zym-zyp)
         ! ... chose the smallest (zd1,zd2)
         IF (zd2 <=  zd1) THEN
            ydpti=ylptmp2   ! use M2
         ELSE
            ydpti=ylptmp1   ! use M1
         END IF
         !
      ELSE
         ! ...  use x(y)
         ! ... possible intermediate points:
         ylptmp1=ydpt(kk-1)+CMPLX(SIGN(1.,pai),0.)  ! M1
         ylptmp2=ydpt(kk-1)+(0.,1.)                 ! M2
         !
         ! ... M1 is the candidate point:
         zxm=REAL(ylptmp1)
         zym=IMAG(ylptmp1)
         za0=pai
         zb0=pbi
         !
         za1=-1./za0
         zb1=zxm - za1*zym
         zyp=-(zb1-zb0)/(za1-za0)
         zxp=za0*zyp + zb0
         zd1=(zxm-zxp)*(zxm-zxp) + (zym-zyp)*(zym-zyp)
         !
         zxm=REAL(ylptmp2)
         zym=IMAG(ylptmp2)
         za1=-1./za0
         zb1=zxm - za1*zym
         zyp=-(zb1-zb0)/(za1-za0)
         zxp=za0*zyp + zb0
         zd2=(zxm-zxp)*(zxm-zxp) + (zym-zyp)*(zym-zyp)
         IF (zd2 <=  zd1) THEN
            ydpti=ylptmp2
         ELSE
            ydpti=ylptmp1
         END IF
      END IF
   END SUBROUTINE interm_pt





END PROGRAM cdficeflux
