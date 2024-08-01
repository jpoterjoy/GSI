module pf
!$$$  module documentation block
!
! module: pf                           Update model state variables with local PF.
!                                      
!
! prgmmr: poterjoy         org: noaa aoml               date: 2017-04-26
!        
! NOTE: This code is modified from the enkf module written by Jeff Whitaker. 
!
! Abstract: This module contains all the necessary code for running the local PF
!   of Poterjoy (2016) and Poterjoy et al. (2019) in GSI. It is based on 
!   the pre-existing enkf module written by Jeff Whitaker. For details on 
!   parallelization strategy and code structure, read description in enkf.f90.
!   For now, weights are calculated for each state variable, but they are only
!   needed for each model grid point. To reduce the computational cost, future 
!   implimentations of the local PF should calculate weights only for grid points.
!
!
! Public Subroutines:
!   pf_update: performs the pf update (calls update_biascorr to perform
!   the bias coefficient update.  The EnKF/bias coefficient update is 
!   iterated numiter times. This 
!
! Public Variables: None
!
! Modules Used: kinds, constants, params, covlocal, mpisetup, loadbal, statevec,
!               kdtree2_module, enkf_obsmod, radinfo, radbias, gridinfo
!
! program history log:
!   2017-04-26:  JPOTERJOY: started modification for local PF 
!   2018-06-01:  JPOTERJOY: finalized code for testing with HWRF
!   2021-01-20:  JPOTERJOY: updated to reflect recent GSI upgrades 
!
! attributes:
!   language: f95
!
!$$$

! JPOTERJOY NOTES: - The current version is hard-coded to use iterative resampling
!                    with a max of three iterations. 
!
!                  - To maintain diversity in particles, the PF pre-calculates "beta"
!                    coefficients for each location in model and obs space. Dividing 
!                    the log of pre-normalized weights by each beta guarantees the 
!                    resulting weights have an effective ensemble size greater than 
!                    or equal to frac_neff*nanals, where frac_neff is specified in the
!                    namelist.
!
!                  - Radiance bias correction coefficients can be calculated
!                    by uncommenting lines indicated below (search for "bias
!                    correction" in iterative resampling loop. This part of the
!                    code has not yet been tested.
!
!                  - The code is not written for efficiency. More work is needed to 
!                    optimize various portions of the filter update steps.


use mpimod, only: mpi_comm_world
use mpisetup, only: mpi_real4,mpi_sum,mpi_comm_io,mpi_in_place,numproc,nproc,&
                mpi_integer,mpi_wtime,mpi_status,mpi_real8,mpi_max,mpi_realkind,&
                mpi_2real,mpi_minloc,mpi_real
use covlocal, only:  taper
use kinds, only: r_double,i_kind,r_single,r_single,r_kind
use kdtree2_module, only: kdtree2_r_nearest, kdtree2_result
use loadbal, only: numobsperproc, numptsperproc, indxproc_obs, iprocob, &
                   indxproc, lnp_chunk, kdtree_obs, kdtree_grid, &
                   ensmean_obchunk, indxob_chunk, oblnp_chunk, nobs_max, &
                   obtime_chunk, grdloc_chunk, obloc_chunk, &
                   npts_max, anal_obchunk_prior, ensmean_chunk, anal_chunk, &
                   ensmean_chunk_prior
use controlvec, only: cvars3d,  ncdim, index_pres
use enkf_obsmod, only: oberrvar, ob, ensmean_ob, obloc, oblnp, &
                  nobstot, nobs_conv, nobs_oz, nobs_sat,&
                  obfit_prior, obfit_post, obsprd_prior, obsprd_post, obtime,&
                  obtype, oberrvarmean, numobspersat, deltapredx, biaspreds,&
                  oberrvar_orig, probgrosserr, prpgerr,&
                  corrlengthsq,lnsigl,obtimel,obloclat,obloclon,stattype
use constants, only: pi, one, zero, two
use params, only: sprd_tol, datapath, nanals,&
                  sortinc,numiter,nlevs,nvars,&
                  zhuberleft,zhuberright,varqc,lupd_satbiasc,huber,univaroz,&
                  nbackgrounds,nhr_anal,fhr_assim,&
                  iseed_perturbed_obs, frac_neff, pf_alpha, pf_kddm, pf_bwnorm, min_res
use radinfo, only: npred,nusis,nuchan,jpch_rad,predx
use radbias, only: apply_biascorr, update_biascorr
use gridinfo, only: nlevs_pres
use sorting, only: quicksort, isort
use mpeu_util, only: getindex

implicit none

private
public :: pf_update

contains

subroutine pf_update()
!------------------------------------------------------------------------
!
!  Main assimilation routine for local particle filter: J. Poterjoy Apr. 2017
! 

use random_normal, only : rnorm, set_random_seed

! local variables.
integer(i_kind) nob,nob1,nob2,npob,nf,nf2,ii,nobx,nskip,i,nrej, &
                niter,maxiter,npt,nuse,ncount,nb,outniter,maxoutiter
real(r_single) dist,lnsig,obt,corrlengthinv,lnsiglinv,obtimelinv,four
real(r_single) corrsqr
real(r_double) :: t1,t2,t3,t4,t5,t6,t7,tbegin,tend
real(r_single) :: term1,term2,term3,pf_alpha2, max_res,m1,m2,m3
real(r_single) r_nanals,r_nanalsm1,neff
real(r_single),allocatable, dimension(:,:) :: anal_obchunk
real(r_single),dimension(nobstot):: oberrvaruse
real(r_single) r
integer(i_kind), allocatable, dimension(:) :: indxassim,iskip
real(r_single), allocatable, dimension(:) :: taper_disob,taper_disgrd
integer(i_kind) ierr
! kd-tree search results
type(kdtree2_result),dimension(:),allocatable :: sresults1,sresults2 
integer(i_kind) nn,nnn,nobm,nsame,nn1,nn2,oz_ind
real(r_single),dimension(nlevs_pres):: taperv
logical lastiter, kdgrid, kdobs, do_update

! New variables for local PF 
character(8)  :: date
character(10) :: time
character(255) :: msgstring
real(r_single) :: ens_var, ens_mean, em, ev, norm, sq_nanals
real(r_single) :: w(nanals), hw(nanals), wn(nanals), ws, taper_coeff, dr, drsq, taper_thresh
real(r_single) :: d(nanals)
real(r_single) :: prior_ens(nanals), post_ens(nanals)
real(r_single) :: t8, wns1, wns2
integer(i_kind) :: indx(nanals), n
real(r_single),parameter :: beta_max = 1e10
real(r_single),parameter :: cutoff = 0.98
real(r_single),allocatable, dimension(:,:,:,:) :: prior_chunk
real(r_single),allocatable, dimension(:) :: ensmean_obchunk_prior, buffertmp
real(r_single),allocatable, dimension(:,:,:,:) :: wo
real(r_single),allocatable, dimension(:,:) :: hwo, hwop
real(r_single),allocatable, dimension(:,:,:) :: beta_x, res_x
real(r_single),allocatable, dimension(:) :: beta_y, res_y
real(r_single),allocatable, dimension(:,:,:) :: ensmean_chunk_prior_temp

allocate(ensmean_chunk_prior_temp(npts_max,ncdim,nbackgrounds))

! allocate temporary arrays.
allocate(anal_obchunk(nanals,nobs_max))
allocate(sresults1(numptsperproc(nproc+1)),taper_disgrd(numptsperproc(nproc+1)))
allocate(sresults2(numobsperproc(nproc+1)),taper_disob(numobsperproc(nproc+1)))
allocate(buffertmp(nobstot))
! index array that controls assimilation order
allocate(indxassim(nobstot),iskip(nobstot))

! New arrays for local PF
allocate(prior_chunk(nanals,npts_max,ncdim,nbackgrounds))
allocate(ensmean_obchunk_prior(nobs_max))

allocate(wo(nanals,npts_max,nlevs_pres,nbackgrounds))
allocate(beta_x(npts_max,nlevs_pres,nbackgrounds))
allocate(hwo(nanals,nobs_max),beta_y(nobs_max))
allocate(hwop(nanals,nobstot))
allocate(res_y(nobs_max))
allocate(res_x(npts_max,nlevs_pres,nbackgrounds))

! define a few frequently used parameters
r_nanals=one/float(nanals)
r_nanalsm1=one/float(nanals-1)
taper_thresh = epsilon(taper_coeff)
sq_nanals=sqrt((float(nanals-1)))
neff = frac_neff*float(nanals)

! initialize some arrays with first-guess values.
obfit_post(1:nobstot) = obfit_prior(1:nobstot)
obsprd_post(1:nobstot) = obsprd_prior(1:nobstot)
anal_obchunk = anal_obchunk_prior

! Check to see if kdtree structures are associated
kdgrid=associated(kdtree_grid)
kdobs=associated(kdtree_obs)

! Initalize obs index vector
do nob=1,nobstot
  indxassim(nob) = nob
end do

! The satellite bias correction loop is borrowed for the PF beta-inflation
! algorithm. Here, the PF is forced to go through at least 2 iterations 
! for sat radiance loop; beta factors for particle weights are estimated 
! on first iteration, thus requiring the minimum number of 2.

! Initialize residuals
res_x = one - min_res
res_y = one - min_res

if (nproc .eq. 0) then
  write(msgstring, '(A,F4.2)') &
       '   PF weights use min_res of ',min_res
  write(6,'(A)') trim(msgstring)
end if

! Initialize vector for indicating whether a low-impact ob should be skipped.
iskip = 0

! The outer iterations are for "iterative resampling." This approach
! break observations into parts, like iterative ensemble kalman smoothers.
maxoutiter=3
!maxoutiter=1

! Tempering loop
temperingloop: do outniter=1,maxoutiter

  if (outniter == 1) then
    maxiter = max(numiter,2)
  else
    maxiter = 2
  endif

  oberrvaruse(1:nobstot) = oberrvar(1:nobstot)

  ! Update first-guess arrays between outer iterations
  prior_chunk = anal_chunk
  anal_obchunk_prior = anal_obchunk
  ensmean_chunk_prior_temp = ensmean_chunk
  ensmean_obchunk_prior = ensmean_obchunk
  do nob1=1,numobsperproc(nproc+1)
    nob2 = indxproc_obs(nproc+1,nob1)
    ensmean_ob(nob2) = ensmean_obchunk(nob1)
  end do

  if (nproc .eq. 0) then
    write(msgstring, '(A,I3,A,I3)') &
                      '   Performing outer PF iteration number ',outniter,' of ',maxoutiter
    write(6,'(A)') trim(msgstring)
  end if

  ! Iterative sat radiance bias calculation; obs are assimilated on last iteration.
  ! First iteration caclulates beta coefficients for preventing weight collapse.
  do niter=1,maxiter

    if (nproc .eq. 0) then
      write(msgstring, '(A,I3,A,I3)') &
                        '   Performing PF inner iteration number ',niter,' of ',maxiter
      write(6,'(A)') trim(msgstring)
    end if

    lastiter = niter == maxiter
    ! Apply bias correction with latest estimate of bias coeffs.
    ! (already done for first iteration)
  
    ! JPOTERJOY: uncomment to apply bias correction
    !if (nobs_sat > 0 .and. numiter > 2) then
    !  if (outniter > 2 .or. niter > 2) call apply_biascorr()
    !end if
  
    ! Reset variables at start of each iteration.
    nrej = 0
    nsame = 0

    ! Obs space mean is updated with latest bias coefficient perturbations.
    ! nob1 is the index of the obs to be processed on this rank
    ! nob2 maps nob1 to 1:nobstot array (nob)
    do nob1=1,numobsperproc(nproc+1)
       nob2 = indxproc_obs(nproc+1,nob1)
       ensmean_obchunk(nob1) = ensmean_ob(nob2)
    enddo
    ensmean_obchunk_prior = ensmean_obchunk
  
    ! Calculate obs-space weights and determines whether to skip obs
    if (niter == 1) then

       ! Initialize obs-space weights
       hwop = r_nanals

       infloop: do nobx=1,nobstot

         ! Current ob and the task owning its prior
         nob = indxassim(nobx)
         npob = iprocob(nob)

         if (nproc == npob) then

           ! Calculate scalar weights used for vector weight formulation
           nob1 = indxob_chunk(nob)

           d = ( ob(nob) - ( anal_obchunk_prior(:,nob1) + ensmean_obchunk_prior(nob1)) )
           d = d**2 / (two * oberrvaruse(nob))
           d = d - minval(d)

           ! Cap min likelihood at 1E-40 to prevent errors with regularization step
           hw = exp( - d/(em) ) + 1E-40
           hw = hw / sum(hw)

!if (nobx < 1000) then
!  write(*,*) ' '
!  write(*,*) 'nobx is',nobx
!  write(*,*) 'Min/Max hw:',minval(hw),maxval(hw)
!  write(*,*) 'innov is:',ob(nob) - ( anal_obchunk_prior(:,nob1) + ensmean_obchunk_prior(nob1))
!  write(*,*) 'Infl:',ws*em
!end if

         end if

         call mpi_bcast(hw,nanals,mpi_real4,npob,mpi_comm_world,ierr)

         ! Everyone gets copies of original obs-space weights
         hwop(1:nanals,nob) = hw

         ! Determine whether to skip ob
         if (one > 0.98_r_single*nanals*sum(hw**2) ) then
           iskip(nob) = 1
         end if

       end do infloop
    end if

    nobm = 1
    ncount = 0
    t2 = zero
    t3 = zero
    t4 = zero
    t5 = zero
    t6 = zero
    t7 = zero
    t8 = zero
    nf    = 0
    nf2   = 0
    tbegin = mpi_wtime()
  
    ! Reset weights each iteration
    hwo = zero
    wo = zero

    ! Reset regularization values every outer iteration
    if (niter == 1) then
      beta_x = one
      beta_y = one
    end if
    
    ! Loop over obs that have not been rejected
    obsloop: do nobx=1,nobstot

! JPOTERJOY
!if (nobx > 10000) then
!  cycle obsloop
!end if


      if (lastiter .and. nproc == 0 .and. mod(nobx, 1000) == 0 ) then
        write(msgstring, '(A,I7)') '   Assimilating ob number ',nobx
        write(6,'(A)') trim(msgstring)
      end if
      if ( niter == 1 .and. nproc == 0 .and. mod(nobx, 1000) == 0 ) then
        write(msgstring, '(A,I7)') '   Processing ob number ',nobx
        write(6,'(A)') trim(msgstring)
      end if

! TEMP
!if (nobx > 3000) cycle obsloop


      t1 = mpi_wtime()

      ! Current ob and the task owning its prior
      nob = indxassim(nobx)
      npob = iprocob(nob)

      ! Skip flagged obs
      if ( iskip(nob) == 1 ) then
       cycle obsloop
      end if

      !-------------------------------------------------------------------
      ! Values needed for time and space localization
      corrsqr = corrlengthsq(nob)
      corrlengthinv=one/corrlengthsq(nob)
      lnsiglinv=one/lnsigl(nob)
      obtimelinv=one/obtimel(nob)
      four=4.0_r_single

      t2 = t2 + mpi_wtime() - t1
      t1 = mpi_wtime()

      !  Only need to recalculate nearest points when lat/lon is different
      if ( nobx == 1 .or. &
          abs(obloclat(nob)-obloclat(nobm)) .gt. taper_thresh .or. &
          abs(obloclon(nob)-obloclon(nobm)) .gt. taper_thresh .or. &
          abs(corrlengthsq(nob)-corrlengthsq(nobm)) .gt. taper_thresh) then

        ! Store ob index for next iteration
        nobm=nob

        ! Find observations and model points near current ob being assimilated
        ! (using a kd-tree for speed).

        ! Obs space part
        nf = 0
        if (kdobs) then
          call kdtree2_r_nearest(tp=kdtree_obs,qv=obloc(:,nob),r2=corrsqr,&
               nfound=nf,nalloc=numobsperproc(nproc+1),results=sresults2)
        else

          ! use brute force search if number of obs on this proc <= 3
          do nob1=1,numobsperproc(nproc+1)
            r = sum( (obloc(:,nob)-obloc_chunk(:,nob1))**2, 1 )
            if (r < corrsqr) then
              nf = nf + 1
              sresults2(nf)%idx = nob1
              sresults2(nf)%dis = r
            end if     
          end do
        end if

        do nob1=1,nf
          ! ozone obs only affect ozone (if univaroz is .true.).
          nob2 = sresults2(nob1)%idx
          if (univaroz .and. obtype(nob)(1:3) .eq. ' oz' .and. obtype(indxproc_obs(nproc+1,nob2))(1:3) .ne. ' oz') then
             taper_disob(nob1) = 0.0_r_single
          else
             dist = sqrt(sresults2(nob1)%dis*corrlengthinv)
             taper_disob(nob1) = taper(dist)
          endif
        end do

        ! Model space part
        nf2=0
        if(niter == 1 .or. lastiter) then

          ! Search analysis grid points for those within corrlength of 
          ! ob being assimilated (using a kd-tree for speed).
          if (kdgrid) then
            call kdtree2_r_nearest(tp=kdtree_grid,qv=obloc(:,nob),r2=corrsqr,&
                  nfound=nf2,nalloc=numptsperproc(nproc+1),results=sresults1)
          else
            ! use brute force search if number of grid points on this proc <= 3
            do npt=1,numptsperproc(nproc+1)
              r = sum( (obloc(:,nob)-grdloc_chunk(:,npt))**2, 1 )
              if (r < corrsqr) then
                nf2 = nf2 + 1
                sresults1(nf2)%idx = npt
                sresults1(nf2)%dis = r
              end if
            end do

          end if

          ! Get localization coefficient corresponding to neighboring grid points
          do nob1=1,nf2
            dist = sqrt(sresults1(nob1)%dis*corrlengthinv)
            taper_disgrd(nob1) = taper(dist)
          end do

        end if

      else

        ! Keep track of how many obs are at same lat/lon location
        nsame=nsame+1

      end if

      t3 = t3 + mpi_wtime() - t1
      t1 = mpi_wtime()

      !-------------------------------------------------------------------

      ! Grab current obs-space weights
      hw = hwop(1:nanals,nob)

      t4 = t4 + mpi_wtime() - t1
      t1 = mpi_wtime()

      ! Perform sampling after first iteration
      if (niter > 1) then

        ! Owner finds indices of particles for sampling
        if (nproc == npob) then

          prior_ens = anal_obchunk_prior(:,nob1) + ensmean_obchunk_prior(nob1)

          ! Sample from weighted prior particles
          w = hwo(1:nanals,nob1) - log(nanals*hw)
          w = w - minval(w)
          w = exp( - w / (beta_y(nob1)) )
          w = w / sum(w)

          ! Determine whether to skip sampling step
          if (one > 0.98_r_single*nanals*sum(w**2) ) then

            indx = 0

          else

            call pf_sample(prior_ens, w, nanals, indx)

          end if

        end if

        call mpi_bcast(indx,nanals,mpi_integer,npob,mpi_comm_world,ierr)

      end if

      t5 = t5 + mpi_wtime() - t1
      t1 = mpi_wtime()

      ! Skip update on first iteration and if Neff exceeds threshold
      if ( niter == 1 .or. indx(1) == 0 ) then
        do_update = .false.
      else
        do_update = .true.
      end if

      ! Keep track of number of obs
      ncount = ncount + 1

      ! only need to update state variables on last iteration.
      oz_ind = getindex(cvars3d, 'oz')
      if (univaroz .and. obtype(nob)(1:3) .eq. ' oz' .and. oz_ind > 0) then ! ozone obs only affect ozone
          nn1 = (oz_ind-1)*nlevs+1
          nn2 = oz_ind*nlevs
      else
          nn1 = 1
          nn2 = ncdim
      end if

      ! Model space update. Called for inflation on first iteration and state update on last iteration.

      ! Redundant part of weight calculation
      w = float(nanals)*hw - one

      if (nf2 > 0) then

!$omp parallel do schedule(dynamic,1) private(ii,i,nb,obt,n,nn,nnn,lnsig,taperv,taper_coeff,prior_ens,post_ens,ens_mean,ens_var,wn,dr,drsq,norm,wns1,wns2,pf_alpha2,em,ev,term1,term2,term3,m1,m2,m3)
        do ii=1,nf2 ! loop over nearby horiz grid points
          do nb=1,nbackgrounds ! loop over background time levels
            !taper_disgrd is horizontal localization and taper(obt*obtimelinv) is time localization
            obt = abs(obtime(nob)-(nhr_anal(nb)-fhr_assim))
            taper_coeff=taper_disgrd(ii)*taper(obt*obtimelinv)
            i = sresults1(ii)%idx

            taperv = zero
            do nn=1,nlevs_pres
              lnsig = abs(lnp_chunk(i,nn)-oblnp(nob))
              !taper(lnsig*lnsiglinv) is vertical localization
              if(lnsig < lnsigl(nob)) then
                taperv(nn)=taper_coeff*taper(lnsig*lnsiglinv)
              end if
            end do

            wo_update: do nn=1,nlevs_pres

              if ( beta_x(i,nn,nb) == beta_max ) cycle wo_update

              taper_coeff = taperv(nn)

              if (taper_coeff .gt. taper_thresh) then

                if (taper_coeff == one) then
                  wn = log(nanals*hw)
                else
                  wn = w*taper_coeff
                  do n = 1,nanals
                    if (abs(wn(n)) > 0.3_r_single) then
                      wn(n) = log( wn(n) + one )
                    end if
                  end do
               end if

                wo(:,i,nn,nb) = wo(:,i,nn,nb) - wn
                wo(:,i,nn,nb) = wo(:,i,nn,nb) - minval(wo(:,i,nn,nb))

              end if

            end do wo_update


            if (do_update) then

              state_update: do nn=nn1,nn2

                nnn=index_pres(nn)
                taper_coeff = taperv(nnn)

                if ( beta_x(i,nnn,nb) == beta_max ) cycle state_update

                if (taper_coeff .gt. taper_thresh) then

! TEMP
!wn = exp(-wo(:,i,nnn,nb)/beta_x(i,nnn,nb))
!wn = wn/sum(wn)
!anal_chunk(1:nanals,i,nn,nb) = wn
!ensmean_chunk(i,nn,nb) =  sum(wn) * r_nanals
!anal_chunk(1:nanals,i,nn,nb) = zero
!ensmean_chunk(i,nn,nb) =  one / sum (wn**2)
!cycle state_update

                 ! Update only when sample variance and localization coefficient are non-zero
                  post_ens = anal_chunk(1:nanals,i,nn,nb) + ensmean_chunk(i,nn,nb)

                  if ( maxval( post_ens ) .ne. minval( post_ens ) ) then

                    ! Update weights and normalization coefficient for update equations
                    wn = exp(-wo(:,i,nnn,nb)/beta_x(i,nnn,nb))
                    wn = wn/sum(wn)

                    ! Calculate posterior mean
                    prior_ens = prior_chunk(:,i,nn,nb) + ensmean_chunk_prior_temp(i,nn,nb)
                    ens_mean = sum(wn*prior_ens)

!                    if (abs(sum(post_ens)*r_nanals - ens_mean) .lt. taper_thresh) cycle state_update
                    if (sum(post_ens)*r_nanals .eq. ens_mean) cycle state_update

   
                    ! Calculate posterior variance and update coefficients
                    dr = (one - taper_coeff) / taper_coeff

                    drsq=dr*dr
                    ens_var = zero
                    norm  = zero
                    term1 = zero
                    term2 = zero
                    term3 = zero

                    do n = 1,nanals
                      ens_var = ens_var + wn(n) * ( prior_ens(n) - ens_mean )**2
                      term1 = term1 + ( prior_ens(indx(n)) - ens_mean )**2
                      term2 = term2 + ( post_ens(n) - ens_mean )**2
                      term3 = term3 + ( post_ens(n) - ens_mean )*( prior_ens(indx(n)) - ens_mean )
                      norm  = norm + wn(n)**2
                    end do
                    ens_var = ens_var / (one - norm)
                    wns1 = term1 + drsq*term2 + two*term3*dr
                    wns2 = drsq/wns1
    
                    ! Alpha reduces update to maintain particle diversity near observation.
                    ! The coeffiecient pf_alpha is specified, while pf_alpha2 is derived to 
                    ! maintain the first two moments given by ens_mean and ens_var
                    wns1 = pf_alpha*sq_nanals*sqrt(ens_var/wns1)
                    wns2 = sq_nanals*sqrt(ens_var*wns2)

                    m1 = sum(prior_ens(indx) - ens_mean)*r_nanals
                    m2 = sum(post_ens - ens_mean)*r_nanals

                    term1 = term1 - nanals*m1**2
                    term2 = term2 - nanals*m2**2
                    term3 = term3 - nanals*m1*m2

                    m1 = term2
                    m2 = two*( wns1*term3 + wns2*term2 )
                    m3 = term1*wns1**2 + term2*wns2**2 + two*term3*wns1*wns2 - (nanals-one)*ens_var

                    pf_alpha2 = ( - m2 + sqrt( m2**2 - four*m1*m3 ) ) / (two*m1)
                    wns2 = wns2 + pf_alpha2

                    ! Update particles
                    post_ens = wns1*(prior_ens(indx) - ens_mean) + wns2*(post_ens - ens_mean)

                    ! Ensure perturbations are zero
                    em = sum(post_ens)*r_nanals
                    anal_chunk(:,i,nn,nb) = post_ens - em
                    ensmean_chunk(i,nn,nb) = ens_mean
  
                  end if ! check if update is needed
                end if ! check if state is within localization region
              end do state_update ! end loop over vertical levels
            end if ! first iteration 
          end do ! end loop over background time levels
        end do ! end loop over nearby horiz grid points
!$omp end parallel do
      end if ! nf2 > 0

      t6 = t6 + mpi_wtime() - t1
      t1 = mpi_wtime()

      ! Obs-space update performed on each iteration
      if (nf > 0) then

!$omp parallel do schedule(dynamic,1) private(nob1,nob2,lnsig,obt,n,taper_coeff,prior_ens,post_ens,ens_mean,ens_var,wn,dr,drsq,norm,wns1,wns2,pf_alpha2,em,ev,term1,term2,term3,m1,m2,m3)
        obs_update: do nob1=1,nf
          nob2 = sresults2(nob1)%idx
          lnsig = abs(oblnp(nob)-oblnp_chunk(nob2))

          if ( beta_y(nob2) == beta_max ) cycle obs_update

          if (lnsig .lt. lnsigl(nob) .and. taper_disob(nob1) .gt. zero) then

            obt = abs(obtime(nob)-obtime_chunk(nob2))
            if (obt .lt. obtimel(nob)) then

              ! Update only when sample variance and localization coefficient are non-zero
              post_ens = anal_obchunk(:,nob2) + ensmean_obchunk(nob2)
              taper_coeff = taper_disob(nob1)*taper(lnsig*lnsiglinv)*taper(obt*obtimelinv)

              if ( maxval( post_ens ) .ne. minval( post_ens ) .and. taper_coeff .gt. taper_thresh ) then

                ! Negative log likelihoods
                if (taper_coeff == one) then
                  wn = log(nanals*hw)
                else
                  wn = w*taper_coeff
                  do n = 1,nanals
                    if (abs(wn(n)) > 0.3_r_single) then
                      wn(n) = log( wn(n) + one )
                    end if
                  end do
                end if

                hwo(:,nob2) = hwo(:,nob2) - wn
                hwo(:,nob2) = hwo(:,nob2) - minval(hwo(:,nob2))

                if (do_update) then

                  ! Update weights and normalization coefficient for update equations
                  wn = exp(-hwo(:,nob2)/beta_y(nob2))
                  wn = wn/sum(wn)

                  ! Calculate posterior mean
                  prior_ens = anal_obchunk_prior(:,nob2) + ensmean_obchunk_prior(nob2)
                  ens_mean = sum(wn*prior_ens)

                  if ( sum(post_ens)*r_nanals .eq. ens_mean ) cycle obs_update

                  ! Calculate posterior variance and update coefficients
                  dr = (one - taper_coeff) / taper_coeff
                  drsq=dr*dr
                  ens_var = zero
                  term1 = zero
                  term2 = zero
                  term3 = zero
                  norm  = zero
                  do n = 1,nanals
                    ens_var = ens_var + wn(n) * ( prior_ens(n) - ens_mean )**2
                    term1 = term1 + ( prior_ens(indx(n)) - ens_mean )**2
                    term2 = term2 + ( post_ens(n) - ens_mean )**2
                    term3 = term3 + ( post_ens(n) - ens_mean )*( prior_ens(indx(n)) - ens_mean )
                    norm = norm + wn(n)**2
                  end do
                  ens_var = ens_var / (one - norm)
                  wns1 = term1 + drsq*term2 + two*term3*dr
                  wns2 = drsq/wns1
    
                  ! Alpha reduces update to maintain particle diversity near observation.
                  ! The coeffiecient pf_alpha is specified, while pf_alpha2 is derived to 
                  ! maintain the first two moments given by ens_mean and ens_var
                  wns1 = pf_alpha*sq_nanals*sqrt(ens_var/wns1)
                  wns2 = sq_nanals*sqrt(ens_var*wns2)

                  m1 = sum(prior_ens(indx) - ens_mean)*r_nanals
                  m2 = sum(post_ens - ens_mean)*r_nanals

                  term1 = term1 - nanals*m1**2
                  term2 = term2 - nanals*m2**2
                  term3 = term3 - nanals*m1*m2

                  m1 = term2
                  m2 = two*( wns1*term3 + wns2*term2 )
                  m3 = term1*wns1**2 + term2*wns2**2 + two*term3*wns1*wns2 - (nanals-one)*ens_var

                  pf_alpha2 = ( - m2 + sqrt( m2**2 - four*m1*m3 ) ) / (two*m1)
                  wns2 = wns2 + pf_alpha2

                  ! Update particles 
                  post_ens = wns1*(prior_ens(indx) - ens_mean) + wns2*(post_ens - ens_mean)

                  ! Ensure perturbations are zero
                  em = sum(post_ens)*r_nanals
                  anal_obchunk(:,nob2) = post_ens - em
                  ensmean_obchunk(nob2) = ens_mean

                end if ! first iteration
              end if ! check if update is needed
            end if ! check if within time cutoff
          end if ! check if in localization region
        end do obs_update ! ob loop
!$omp end parallel do
      end if ! no close obs.

      t7 = t7 + mpi_wtime() - t1
      t1 = mpi_wtime()
  
! JPOTERJOY: sanity check
!if (nproc .eq. npob .and. lastiter) then
!  nob1 = indxob_chunk(nob)
!  print *, 'ob:',ob(nob)
!  print *, 'ob err:',oberrvaruse(nob)
!  print *, 'ob infl:',obs_err_infl(nob1)
!  print *, 'hwo:',hwo(:,nob1)
!  print *, 'hw:',hw
!  print *, 'indx:',indx
!  print *, 'hx prior:',anal_obchunk_prior(:,nob1) + ensmean_obchunk_prior(nob1)
!  print *, 'hx post:',anal_obchunk(:,nob1) + ensmean_obchunk(nob1)
!end if
 
    end do obsloop ! loop over obs to assimilate

 
    if (niter == 1) then

      ! Calculate regularization terms
      if (nproc == 0) then
        write(6,*) 'Calculating obs space regularization factors'
      end if

!$omp parallel do schedule(dynamic,1) private(nob1,ws)

      do nob1=1,numobsperproc(nproc+1)

        if ( res_y(nob1) > zero ) then

          call pf_regularization(hwo(:,nob1), nanals, neff, ws, beta_max)
          beta_y(nob1) = ws

          ! Track residual of beta factors
          if (res_y(nob1) <= one/beta_y(nob1)) then
             beta_y(nob1) = one/res_y(nob1)
             res_y(nob1) = zero
          else
             res_y(nob1) = res_y(nob1) - one/beta_y(nob1)
          end if

        else

          beta_y(nob1) = beta_max

        end if

      end do

!$omp end parallel do
  
      ! Calculate regularization terms
      if (nproc == 0) then
        write(6,*) 'Calculating state space regularization factors'
      end if

      max_res = zero

!$omp parallel do schedule(dynamic,1) private(nb,nn,i,ws) reduction(max:max_res)

     do i=1,npts_max
       do nn=1,nlevs_pres
         do nb=1,nbackgrounds

            if ( res_x(i,nn,nb) > zero ) then

              call pf_regularization(wo(:,i,nn,nb), nanals, neff, ws, beta_max)
              beta_x(i,nn,nb) = ws

              ! Track residual of beta factors
              if (res_x(i,nn,nb) <= one/beta_x(i,nn,nb)) then
                beta_x(i,nn,nb) = one/res_x(i,nn,nb)
                res_x(i,nn,nb) = zero
              else
                res_x(i,nn,nb) = res_x(i,nn,nb) - one/beta_x(i,nn,nb)
              end if

              max_res = max(res_x(i,nn,nb),max_res)

            else

              beta_x(i,nn,nb) = beta_max

            end if

          end do
        end do
      end do

!$omp end parallel do

      ! Find global max residual
      call mpi_allreduce(max_res,ws,1,mpi_real4,mpi_sum,mpi_comm_world,ierr)
      max_res = ws

    end if  ! first iteration
    
    ! JPOTERJOY: sanity check
    !if (nproc .eq. 0) then
    !  do n = 1,nobs_max
    !    hw = one / sum( hwo(:,n)**2 )
    !    print *, n,'min Neff:',minval(hw),'sum hwo:',sum(hw)
    !  end do
    !end if
  
    ! Additional corrections to state variables using KDDM
    if (pf_kddm > 0 .and. lastiter) then
  
      if (nproc .eq. 0) then
        call date_and_time( date, time )
        write(msgstring,*) 'KDDM begin time: ',time(1:2),':',time(3:4),':',time(5:6)
        write(6,'(A)') trim(msgstring)
      endif

      ! Loop through everything in obs array
!$omp parallel do schedule(dynamic,1) private(nob1,prior_ens,post_ens,wn)
      obs_kddm: do nob1=1,numobsperproc(nproc+1)

        if ( beta_y(nob1) == beta_max ) cycle obs_kddm

        ! Update only when posterior ensemble is outside span of prior ensemble
        prior_ens = anal_obchunk_prior(:,nob1) + ensmean_obchunk_prior(nob1)
        post_ens = anal_obchunk(:,nob1) + ensmean_obchunk(nob1)

        if ( maxval(post_ens) > maxval(prior_ens) .or. minval(post_ens) < minval(prior_ens) ) then
                
          wn = exp(-hwo(:,nob1)/beta_y(nob1))
          wn = wn/sum(wn)

          if ( maxval(wn) .ne. minval(wn) ) then

            call pf_kddm_update(post_ens,prior_ens,wn,nanals,pf_bwnorm)

            if ( maxval(post_ens) .ne. minval(post_ens) ) then
              ensmean_obchunk(nob1) = sum(post_ens)*r_nanals
              anal_obchunk(:,nob1) = post_ens - ensmean_obchunk(nob1)
            endif

          endif

        endif

      end do obs_kddm
!$omp end parallel do
  
      ! Loop through everything in state array
!$omp parallel do schedule(dynamic,1) private(i,nn,nnn,nb,prior_ens,post_ens,wn)
      do i=1,npts_max

        do nb=1,nbackgrounds        

          state_kddm: do nn=1,ncdim

            ! Get index in weighting matrix
            nnn=index_pres(nn)

            if ( beta_x(i,nnn,nb) == beta_max ) cycle state_kddm

            ! Update only when posterior ensemble is outside span of prior ensemble
            prior_ens = prior_chunk(:,i,nn,nb) + ensmean_chunk_prior_temp(i,nn,nb)
            post_ens = anal_chunk(:,i,nn,nb) + ensmean_chunk(i,nn,nb)

            if ( maxval(post_ens) > maxval(prior_ens) .or. minval(post_ens) < minval(prior_ens) ) then

              wn = exp(-wo(:,i,nnn,nb)/beta_x(i,nnn,nb))
              wn = wn/sum(wn)

              if ( maxval(wn) .ne. minval(wn) ) then

                call pf_kddm_update(post_ens,prior_ens,wn,nanals,pf_bwnorm)

                if ( maxval(post_ens) .ne. minval(post_ens) ) then
                  ensmean_chunk(i,nn,nb) = sum(post_ens)*r_nanals
                  anal_chunk(:,i,nn,nb) = post_ens - ensmean_chunk(i,nn,nb)
                endif

              endif

            endif

          enddo state_kddm

        enddo

      enddo
!$omp end parallel do

      if (nproc .eq. 0) then
        call date_and_time( date, time )
        write(msgstring,*) 'KDDM end time: ',time(1:2),':',time(3:4),':',time(5:6)
        write(6,'(A)') trim(msgstring)
      endif

    endif ! KDDM step
  
    t8 = t8 + mpi_wtime() - t1
  
    ! Additional bookkeeping following serial assimilation
    tend = mpi_wtime()
    if (nproc .eq. 0) then

      write(6,8003) niter,'timing on proc',nproc,' = ',tend-tbegin,t2,t3,t4,t5,t6,t7,t8

      nuse = 0
      do nob1=1,ncount
        nob = indxassim(nob1)
        if (iskip(nob) .ne. 1) then
          nuse = nuse + 1
        endif
      enddo
      nskip = nobstot-nuse
  
      if (nskip > 0) print *,nskip,' out of',nobstot,'obs skipped,',nuse,' obs used'
      if (nsame > 0) print *,nsame,' out of',nobstot-nskip,' same lat/long'
      if (nrej >  0) print *,nrej,' obs rejected by varqc'

    endif
    8003  format(i2,1x,a14,1x,i5,1x,a3,9(f7.2,1x))

    t1 = mpi_wtime()

    ! distribute the O-A stats to all processors.
    buffertmp=zero
    do nob1=1,numobsperproc(nproc+1)
      nob2=indxproc_obs(nproc+1,nob1)
      buffertmp(nob2) = ensmean_obchunk(nob1)
    end do

    call mpi_allreduce(buffertmp,obfit_post,nobstot,mpi_real4,mpi_sum,mpi_comm_world,ierr)

    ! JPOTERJOY: sanity check
    !if (nproc == 0 .and. lastiter) then
    !  do nob1=1,nobstot
    !    write(*,*)'nob: ',nob1,'hx:',obfit_post(nob1),'y:',ob(nob1)
    !  end do
    !end if
  
    obfit_post = ob - obfit_post
    if (nproc == 0) print *,'time to broadcast obfit_post = ',mpi_wtime()-t1,' secs, niter =',niter
  
    ! satellite bias correction update.

    ! NOTE: Coefficients are only estimated on first outer iteration.
  
    ! JPOTERJOY: uncomment for satellite bias correction
    !  if (nobs_sat > 0 .and. lupd_satbiasc .and. niter > 1 .and. outniter == 1) call update_biascorr(niter)  
  
  enddo ! niter loop

  if (max_res == zero) then
    if (nproc == 0) then
      write(msgstring, '(A,I3)') &
                        '   Total number of PF iterations: ',outniter
      write(6,'(A)') trim(msgstring)
    end if
    exit temperingloop
  end if

enddo temperingloop ! tempering loop


! TEMP
!!$omp parallel do schedule(dynamic,1) private(prior_ens,i,nb,nn)
!do i=1,npts_max
!  do nb=1,nbackgrounds
!    do nn=1,ncdim
!      prior_ens = wo(1:nanals,i,nn,nb)
!      anal_chunk(1:nanals,i,nn,nb) = prior_ens
!      ensmean_chunk(i,nn,nb) = sum(prior_ens)*r_nanals
!    end do
!  end do
!end do
!!$omp end parallel do

! distribute the HPaHT stats to all processors.
t1 = mpi_wtime()
buffertmp=zero
do nob1=1,numobsperproc(nproc+1)
  nob2=indxproc_obs(nproc+1,nob1)
  buffertmp(nob2) = sum(anal_obchunk(:,nob1)**2)*r_nanalsm1
end do

call mpi_allreduce(buffertmp,obsprd_post,nobstot,mpi_real4,mpi_sum,mpi_comm_world,ierr)

if (nproc == 0) print *,'time to broadcast obsprd_post = ',mpi_wtime()-t1

predx = predx + deltapredx ! add increment to bias coeffs.

! free local temporary arrays.
deallocate(ensmean_chunk_prior_temp)
deallocate(taper_disob,taper_disgrd)
deallocate(prior_chunk)
deallocate(beta_x)
deallocate(beta_y)
deallocate(res_x)
deallocate(res_y)
deallocate(wo)
deallocate(hwo)
deallocate(hwop)

! these allocated in loadbal, no longer needed
if (min_res == 0.0) then
  deallocate(anal_obchunk)
  deallocate(anal_obchunk_prior)
end if
deallocate(sresults1,sresults2)
deallocate(indxassim,buffertmp)

end subroutine pf_update

! JPOTERJOY: additional subroutines for local PF

subroutine pf_regularization(lw, ens_size, Neff, beta, beta_max)
!------------------------------------------------------------------------
!
!  Calculate regularization factors for particle filter: J. Poterjoy Jun. 2019
! 

integer,  intent(in)          :: ens_size
real(r_single), intent(in)    :: Neff, beta_max, lw(ens_size)
real(r_single), intent(out)   :: beta

real(r_single) :: Neff_init, Neff_final, ke, km, ks, ws
real(r_single) :: tol, fke, fkm, fks, w(ens_size)
integer  :: i

! Initial weights and Neff
w = exp(-lw)
ws = sum(w)
w = w/ws
Neff_init = one / sum( w**2 )

! Inflate if effective ensemble size is smaller than threshold
if ( ( Neff_init < Neff ) .or. ( ws == zero ) ) then

   ! Initial start and end bounds
   ks = 1
   ke = max(1.0_r_single,10*maxval(lw))

   ! Apply bisection method to solve for k
   tol = 1E-4_r_single
   do i = 1,1000
 
      ! Mid point
      km = (ke + ks) / two
  
      ! Evaluate function at end points
      w = exp(-lw/ks)
      if (sum(w) == zero) then
         fks = Neff - one
      else
         w = w / sum(w)
         fks = Neff - one / sum(w**2)
      end if

      w = exp(-lw/ke)
      w = w / sum(w)
      fke = Neff - one / sum(w**2)

      ! Evaluate function at mid points
      w = exp(-lw/km)
      if (sum(w) == zero) then
         fkm = Neff - one
      else
         w = w / sum(w)
         fkm = Neff - one / sum(w**2)
      end if

      ! Exit critera
      if ( abs(ke-ks) < tol ) exit
 
      ! New end points 
      if ( fkm * fks > zero ) then
        ks = km
      else
        ke = km
      end if

   end do

   beta = km

   ! Underflow errors can still lead to wrong result
   w = exp( -lw/beta)
   w = w / sum(w)
   Neff_final = one / sum(w**2)
   ! Target Neff is not always obtainable when multiple members have zero weights.
   ! Set beta to max value when min value of 2 is not reached.
   if (Neff_final < Neff - 0.1_r_single) then
      beta = beta_max
      write(*,*) 'Warning: setting beta to beta_max'
      write(*,*) 'Neff:',Neff_final
      write(*,*) 'min w:',minval(w)
      write(*,*) 'neff iter:',minval(w)
  end if

  ! Sanity check
  !write(*,*) ' Starting Neff: ',Neff_init,' Target Neff: ',Neff,'New Neff: ',Neff_final,'beta: ',beta,'iterations: ',i

else

   beta = one

end if


end subroutine pf_regularization


subroutine pf_regularization_minw(lw, minwt, ens_size, beta)
!------------------------------------------------------------------------
!
!  Calculate regularization factors for particle filter: J. Poterjoy Jun. 2019
! 

integer,  intent(in)    :: ens_size
real(r_single), intent(in)    :: lw(ens_size), minwt
real(r_single), intent(out)   :: beta

real(r_single) :: ke, km, ks, ws
real(r_single) :: tol, fke, fkm, fks, w(ens_size)
integer  :: i

! Initial weights and minw
w = exp(-lw)
ws = sum(w)
w = w/ws

! Inflate if min normalized weight is smaller than threshold
if ( ( minval(w) < minwt ) .or. ( ws == zero ) ) then

   ! Initial start and end bounds
   ks = 1
   ke = maxval(lw)
  
   ! Apply bisection method to solve for k
   tol = 0.001_r_single

   do i = 1,1000
 
      ! Mid point
      km = (ke + ks) / two
  
      ! Evaluate function at end points
      w = exp(-lw/ks)
      if (sum(w) == zero) then
         fks = minwt
      else
         w = w / sum(w)
         fks = minwt - minval(w)
      end if

      w = exp(-lw/ke)
      w = w / sum(w)
      fke = minwt - minval(w)

      ! Evaluate function at mid points
      w = exp(-lw/km)
      w = w / sum(w)
      fkm = minwt - minval(w)

      ! Exit critera
      if ( abs(ke-ks)/two < tol ) exit
 
      ! New end points 
      if ( fkm * fks > zero ) then
        ks = km
      else
        ke = km
      end if

   end do

   beta = km

   w = exp(-lw/beta)
   w = w / sum(w)

   !Sanity check
   !write(*,*) ' Target min w: ',minwt,'New min w: ',minval(w),'beta: ',beta,'iterations: ',i

else

   beta = one

end if


end subroutine pf_regularization_minw


subroutine pf_sample(ens, w, ens_size, indx2)
!------------------------------------------------------------------------
!
!  Perform sampling step of particle filter: J. Poterjoy Apr. 2017
! 

integer(i_kind),  intent(in)    :: ens_size
real(r_single),   intent(in)    :: w(ens_size), ens(ens_size)
integer(i_kind),  intent(out)   :: indx2(ens_size)

real(r_single) :: cw(0:ens_size), base, frac
integer(i_kind)  :: i, j, indx1(ens_size), m, ind(ens_size)

! Find sorting indices and sort weights
call quicksort(ens_size, real(ens), ind)

! Perform deterministic resampling
cw(0) = zero
do i = 1, ens_size
   cw(i) = cw(i - 1) + w(ind(i))
end do

! Divide interval into ens_size parts and choose new particles
! based on the interval they accumulate in
base = one / ens_size / two

j = 1
do i = 1, ens_size

   frac = base + (i - one) / ens_size

   ! Search in the cumulative range to see where frac falls
   m = 0
   do while (m == 0)
     if(cw(j - 1) < frac .and. frac <= cw(j)) then
       indx1(i) = j
       m = 1
     else
       j = j + 1
     end if
   end do

end do

! Unsort indices
indx1 = ind(indx1)

! If a particle is removed, it is replaced by a duplicated
! particle. This is accomplished by looping through indx1
! and flagging replicated indices with a zero, and
! indicating their location in indx2

! Locate the removed indices in indx1
do i = 1, ens_size

   ! Locate first occurance of index i in indx1
   m = minloc(indx1, 1, mask=indx1.eq.i)

   if (m == 0) then
      ! If i is not in indx1, flag the index with a zero in indx2
      indx2(i) = 0
   else
      ! If i is in indx1, indicate value in indx2
      indx2(i) = i
      ! Flag value in indx1 with a zero to show it was removed
      indx1(m) = 0
   endif

end do

! Replace the removed indices with duplicated ones
do i = 1, ens_size

  if (indx2(i) == 0) then
    do m = 1,ens_size
      if (indx1(m) /= 0) exit
    end do
    indx2(i) = indx1(m)
    indx1(m) = 0
  endif

end do

end subroutine pf_sample


subroutine pf_kddm_update(ens1, ens2, w, ens_size, pf_bwnorm)
!------------------------------------------------------------------------
!
!  Apply kernel density distribution mapping method proposed by Seth McGinnis
!  to map a sample of particles into posterior particles:  J. Poterjoy Apr. 2017
! 

integer(i_kind),  intent(in)  :: ens_size
real(r_single),   intent(in)  :: ens2(ens_size)
real(r_single),intent(inout)  :: ens1(ens_size)
real(r_single),    intent(in) :: w(ens_size), pf_bwnorm
integer(i_kind),  parameter   :: npoints = 300
integer                       :: i, j, m, ind(ens_size)
real(r_single)                :: xd(npoints), cda(npoints), qf(ens_size), x(ens_size)
real(r_single)                :: w2, w1, d, q

! NOTE: Specification of npoints needs to consider bandwidth and range of domain
!       for cdfs estimated during probability mapping


! Use kernels to approximate quantiles and posterior cdf
call pf_get_q_cda(ens1,ens2,ens_size,npoints,w,xd,qf,cda,pf_bwnorm)

!! Invert posterior cdf to find values at prior quantiles
!do i = 1,ens_size
!
!   q = qf(i)
!
!   if ( q <= cda(1) ) then 
!      x(i) = xd(1)
!   else if ( q >= cda(npoints) ) then 
!      x(i) = xd(npoints)
!   else
!
!      m = minloc(cda, 1, mask=cda.gt.q)
!
!      if ( q == cda(m) .or. m == 1) then
!         x(i) = xd(m)
!      else
!         d = cda(m) - cda(m-1)
!         if ( d > 1E-10 ) then
!            w1 = ( cda(m) - q ) / d
!            w2 = ( q - cda(m-1) ) / d
!            x(i) = w1 * xd(m-1) + w2 * xd(m)
!         else
!            x(i) = xd(m-1)
!         end if
!      end if
!
!   end if
!
!end do


!! Sorting indices for quantiles
call quicksort(ens_size, qf, ind)

! Invert posterior cdf to find values at prior quantiles
m = 1
do i = 1,ens_size

  ! Get index of sorted quantiles
  j = ind(i)
  q = qf(j)

  ! Check if qf is on edge of domain 
  if ( q <= cda(1) ) then
    x(j) = xd(1)
  else if ( q >= cda(npoints) ) then
    x(j) = xd(npoints)
  else 

    ! Advance index for matching quantile in target cdf
    do while (cda(m) < q )
      m = m + 1
    end do

    ! Interpolate to get x
    d = cda(m) - cda(m-1)
    if ( d > 1E-10_r_single ) then
      w1 = ( cda(m) - q ) / d
      w2 = ( q - cda(m-1) ) / d
      x(j) = w1 * xd(m-1) + w2 * xd(m)
    else
      x(j) = xd(m)
    end if

  end if

end do

! Return remapped samples in x
ens1 = x

end subroutine pf_kddm_update


subroutine pf_get_q_cda(ens1,ens2,ens_size,npoints,w,x,qf,cda,pf_bwnorm)
!------------------------------------------------------------------------
!
! Gaussian kernel density estimation:  J. Poterjoy Apr. 2017
! 
! This subroutine returns prior quantiles and posterior cdf
! estimated using Gaussian kernels.

integer(i_kind),intent(in)    :: ens_size, npoints
real(r_single), intent(in)    :: ens1(ens_size), ens2(ens_size)
real(r_single), intent(in)    :: w(ens_size), pf_bwnorm
real(r_single), intent(out)   :: x(npoints)
real(r_single), intent(out)   :: qf(ens_size), cda(npoints)
integer(i_kind)               :: i
real(r_single)                :: bw, xmin, xmax, v2, range

! Minimum bandwidth is set to be sample standard deviation
v2 = sum(ens2)/real(ens_size)
v2 = sum( ( ens2 - v2 )**2 ) / real(ens_size - 1)
bw = sqrt(v2)*pf_bwnorm
!bw = ( maxval(ens1) - minval(ens1) ) / 6.0_r_single

! Domain for calculating posterior cdf
xmin = min(minval(ens2),minval(ens1))
xmax = max(maxval(ens2),maxval(ens1))
range = xmax-xmin
do i=1,npoints
  x(i) = xmin + (i-one)*range/(npoints-one)
end do

! Estimate quantiles and cdfs by taking sum over Gaussian cdfs
qf  = zero
cda = zero
do i = 1,ens_size

   ! Prior quantiles
   qf = qf + ( one + erf( (ens1 - ens1(i) )/sqrt(two)/bw ) )/two/ens_size

   ! Posterior cdf
   cda = cda + w(i) * ( one + erf( (x - ens2(i) )/sqrt(two)/bw ) )/two

end do

end subroutine pf_get_q_cda


end module pf
