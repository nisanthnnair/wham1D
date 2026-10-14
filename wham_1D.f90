!***********************************************!
! 1D WHAM (Weighted Histogram Analysis Method)
! Stitches 1D biased histograms from umbrella sampling.
! Umbrella bias: (1/2) k (cv - cv0)^2 along a single CV.
!
! Optional:
!   itseries = 1 : build normalized histograms from time series
!                  (time, cv); lines starting with # ignored
!   nboot    > 0 : Poisson/multinomial bootstrap error on F(s)
!***********************************************!
program fes_calc
implicit none
integer i, iter, ncv, umbr_n, nbin1, iwrap, max_iter, ios
integer itseries, nboot, iboot, iref, quiet
real*8, allocatable :: grid0(:,:), prob(:), biased_prob(:,:), grid(:,:)
real*8, allocatable :: hist_count(:,:), dF(:), Fref(:), Fboot(:,:)
real*8, allocatable :: biased_boot(:,:), prob_boot(:), a_boot(:)
real*8 :: kt, toler, pnorm
logical :: parent
integer :: i_umbr, i_s1, isp, ibin, ios_hist, ios2, npts
character (len=256) cvfile, f1, line
real*8 cnvg, cvval, pval, tval, meanF, varF
real*8, allocatable :: umbr_mean(:), umbr_k(:), a(:)
integer, allocatable :: nmax(:)
integer :: rank, gleng1_min, gleng1_max
integer :: ibuf(1)
real*8 :: rbuf(1)
real*8, parameter :: kb=1.9872041E-3 !kcal K-1 mol-1
real*8, parameter :: au_to_kcal = 627.51
real*8, parameter :: kj_to_kcal = 0.239006
integer, parameter :: max_iter_default = 2000
real*8, parameter :: norm_tol = 0.1d0

CALL WHAM_MPI_Init
CALL Set_Parent(parent)

max_iter=max_iter_default
itseries=0
nboot=0

if(parent)then
  open (unit=1, file='input',status='old')
  ! ncv, umbr_n, T, iwrap [, max_iter [, itseries [, nboot]]]
  read(1,'(a)') line
  read(line,*,iostat=ios) ncv, umbr_n, kt, iwrap, max_iter, itseries, nboot
  if(ios.ne.0)then
     read(line,*,iostat=ios) ncv, umbr_n, kt, iwrap, max_iter, itseries
     if(ios.ne.0)then
        read(line,*,iostat=ios) ncv, umbr_n, kt, iwrap, max_iter
        if(ios.ne.0)then
           read(line,*) ncv, umbr_n, kt, iwrap
           max_iter=max_iter_default
        end if
        itseries=0
     end if
     nboot=0
  end if
  if (ncv.gt.1) STOP 'ncv > 1 is not supported in 1D WHAM'
  if ((iwrap.ne.0).and.(iwrap.ne.1)) STOP 'iwrap must be 0 (non-periodic) or 1 (2pi periodic)'
  if (max_iter.lt.1) STOP 'max_iter must be >= 1'
  if ((itseries.ne.0).and.(itseries.ne.1)) STOP 'itseries must be 0 (histogram) or 1 (time series)'
  if (nboot.lt.0) STOP 'nboot must be >= 0'
end if
  kt=kb*kt

ibuf(1)=ncv;      CALL IBcast(ibuf,1); ncv=ibuf(1)
ibuf(1)=umbr_n;   CALL IBcast(ibuf,1); umbr_n=ibuf(1)
rbuf(1)=kt;       CALL RBcast(rbuf,1); kt=rbuf(1)
ibuf(1)=iwrap;    CALL IBcast(ibuf,1); iwrap=ibuf(1)
ibuf(1)=max_iter; CALL IBcast(ibuf,1); max_iter=ibuf(1)
ibuf(1)=itseries; CALL IBcast(ibuf,1); itseries=ibuf(1)
ibuf(1)=nboot;    CALL IBcast(ibuf,1); nboot=ibuf(1)

if (ncv.ne.1) STOP 'This code implements 1D WHAM only (ncv = 1)'

allocate(grid(3,ncv))
allocate(grid0(3,ncv))
allocate(a(umbr_n))
allocate(umbr_mean(umbr_n))
allocate(umbr_k(umbr_n))
allocate(nmax(umbr_n))
a=1.d0
nmax=0

if (parent) then
 do i=1,ncv
    read(1,*)grid0(1:3,i)
    write(*,'(i10,3f16.6)')i,grid0(1:3,i)
 end do
end if

if (parent) then
nbin1=nint((grid0(2,1)-grid0(1,1))/grid0(3,1))+1
end if

CALL RBcast(grid0(1,1),3*ncv)
ibuf(1)=nbin1; CALL IBcast(ibuf,1); nbin1=ibuf(1)

allocate(biased_prob(nbin1,umbr_n))
allocate(hist_count(nbin1,umbr_n))
allocate(prob(nbin1))
allocate(dF(nbin1))
biased_prob=0.d0
hist_count=0.d0
prob=0.d0
dF=0.d0

if (parent) then
   call random_seed()
   open( unit =2, file= 'whaminput', status = 'old')
   read(2,*)toler
   do i_umbr=1, umbr_n
      write(*,*) 'umbrella simulation #', i_umbr
      read(2,'(a)') line
      line=adjustl(line)
      isp=index(line,' ')
      if(isp.gt.0)then
         cvfile=line(1:isp-1)
      else
         cvfile=trim(line)
      end if
      write(*,*)'CVFILE',trim(cvfile)

      read(2,'(a)') line
      if(itseries.eq.1)then
         read(line,*) umbr_mean(i_umbr),umbr_k(i_umbr)
         nmax(i_umbr)=0
      else
         read(line,*) umbr_mean(i_umbr),umbr_k(i_umbr),nmax(i_umbr)
      end if
      umbr_k(i_umbr)=umbr_k(i_umbr)*kj_to_kcal
      if(iwrap.eq.1)then
         if(umbr_mean(i_umbr).gt. 3.14d0)umbr_mean(i_umbr)=umbr_mean(i_umbr)- 6.28d0
         if(umbr_mean(i_umbr).lt.-3.14d0)umbr_mean(i_umbr)=umbr_mean(i_umbr)+ 6.28d0
      end if

      f1 = cvfile
      open( unit=3, file=f1, status='old' )
      biased_prob(:,i_umbr)=0.0d0
      hist_count(:,i_umbr)=0.0d0

      if(itseries.eq.1)then
         ! Build normalized histogram from time series: time, cv
         npts=0
         do
            read(3,'(a)',iostat=ios_hist) line
            if(ios_hist.ne.0) exit
            line=adjustl(line)
            if(len_trim(line).eq.0) cycle
            if(line(1:1).eq.'#') cycle
            read(line,*,iostat=ios2) tval, cvval
            if(ios2.ne.0) cycle
            npts=npts+1
            ibin=nint((cvval-grid0(1,1))/grid0(3,1))+1
            if(ibin.ge.1 .and. ibin.le.nbin1) then
               hist_count(ibin,i_umbr)=hist_count(ibin,i_umbr)+1.0d0
            end if
         end do
         nmax(i_umbr)=npts
         write(*,*) ' time-series points used (nmax) =', nmax(i_umbr)
         if(nmax(i_umbr).lt.1)then
            write(*,*) 'ERROR: no data points in ', trim(cvfile)
            STOP 'empty time series'
         end if
         do i_s1=1,nbin1
            biased_prob(i_s1,i_umbr)=hist_count(i_s1,i_umbr) &
                 /(dfloat(nmax(i_umbr))*grid0(3,1))
         end do
      else
         ! Precomputed histogram: cv, probability
         do
            read(3,*,iostat=ios_hist) cvval, pval
            if(ios_hist.ne.0) exit
            ibin=nint((cvval-grid0(1,1))/grid0(3,1))+1
            if(ibin.ge.1 .and. ibin.le.nbin1) then
               biased_prob(ibin,i_umbr)=pval
            end if
         end do
         ! Effective counts for bootstrap: nmax * P * dcv
         do i_s1=1,nbin1
            hist_count(i_s1,i_umbr)=dfloat(nmax(i_umbr)) &
                 *biased_prob(i_s1,i_umbr)*grid0(3,1)
         end do
      end if
      close(3)

      pnorm=0.0d0
      do i_s1=1,nbin1
         pnorm=pnorm+biased_prob(i_s1,i_umbr)
      end do
      pnorm=pnorm*grid0(3,1)
      write(*,'(A,I5,A,E16.8)') ' umbrella', i_umbr, &
           '  int P(s) ds (before normalization) =', pnorm
      if(dabs(pnorm-1.0d0).gt.norm_tol)then
         write(*,*) 'WARNING: int P(s) ds is not 1 within tolerance ', norm_tol
         write(*,*) ' file = ', trim(cvfile)
         write(*,*) ' performing brute-force normalization: P(s) <- P(s) / int P(s) ds'
         if(pnorm.le.0.0d0)then
            write(*,*) 'ERROR: int P(s) ds <= 0; cannot normalize'
            STOP 'non-positive probability integral'
         end if
         do i_s1=1,nbin1
            biased_prob(i_s1,i_umbr)=biased_prob(i_s1,i_umbr)/pnorm
         end do
         ! keep hist_count consistent with normalized P
         do i_s1=1,nbin1
            hist_count(i_s1,i_umbr)=dfloat(nmax(i_umbr)) &
                 *biased_prob(i_s1,i_umbr)*grid0(3,1)
         end do
         write(*,'(A,I5,A,E16.8)') ' umbrella', i_umbr, &
              '  int P(s) ds (after normalization) =', 1.0d0
      end if
   enddo
   close(2)
   close(1)
end if

rbuf(1)=toler; call RBcast(rbuf,1); toler=rbuf(1)
call RBcast(umbr_mean,umbr_n)
call RBcast(umbr_k,umbr_n)
call IBcast(nmax,umbr_n)
call RBcast(biased_prob(1,1),nbin1*umbr_n)
call RBcast(hist_count(1,1),nbin1*umbr_n)

 if (parent) then
    write(*,*) 'wham begins'
    write(*,*) 'max_iter =', max_iter
    write(*,*) 'itseries =', itseries, ' (0=histogram, 1=time series)'
    write(*,*) 'nboot    =', nboot, ' (0=no FE error bars)'
    if(iwrap.eq.1)then
       write(*,*) 'iwrap = 1 : applying 2pi periodic wrapping'
    else
       write(*,*) 'iwrap = 0 : no periodic wrapping'
    end if
 end if
 CALL DistributeGrids(ncv,grid0,grid,rank,gleng1_min,gleng1_max)

 write(*,*)'new_grid', gleng1_min, gleng1_max, rank, grid(1,1)

 quiet=0
 call run_wham(a,umbr_k,umbr_mean,nmax,grid0,biased_prob,umbr_n,nbin1,ncv,kt, &
      prob,cnvg,gleng1_min,gleng1_max,iwrap,toler,max_iter,parent,quiet)

 ! reference free energy on full grid
 allocate(Fref(nbin1))
 call gather_prob_to_F(nbin1,grid0,kt,prob,Fref,parent)

 !---------- bootstrap error estimate ----------
 ! align reference FE to its minimum
 ! (save Fmin first: in-place Fref(i)=Fref(i)-Fref(iref) zeros the
 !  reference bin too early and leaves later bins unshifted)
 iref=1
 do i_s1=2,nbin1
    if(Fref(i_s1).lt.Fref(iref)) iref=i_s1
 enddo
 meanF=Fref(iref)
 do i_s1=1,nbin1
    Fref(i_s1)=Fref(i_s1)-meanF
 enddo

 if(nboot.gt.0)then
    if(parent) write(*,*) 'bootstrap FE error: nboot =', nboot
    allocate(Fboot(nbin1,nboot))
    allocate(biased_boot(nbin1,umbr_n))
    allocate(prob_boot(nbin1))
    allocate(a_boot(umbr_n))
    Fboot=0.d0
    quiet=1
    do iboot=1,nboot
       if(parent)then
          write(*,*) ' bootstrap sample ', iboot, ' / ', nboot
          call bootstrap_histograms(nbin1,umbr_n,nmax,grid0(3,1), &
               hist_count,biased_boot)
       end if
       call RBcast(biased_boot(1,1),nbin1*umbr_n)
       a_boot=1.d0
       call run_wham(a_boot,umbr_k,umbr_mean,nmax,grid0,biased_boot,umbr_n, &
            nbin1,ncv,kt,prob_boot,cnvg,gleng1_min,gleng1_max,iwrap,toler, &
            max_iter,parent,quiet)
       call gather_prob_to_F(nbin1,grid0,kt,prob_boot,Fboot(1,iboot),parent)
       ! align to same reference bin as Fref (save value before overwriting)
       meanF=Fboot(iref,iboot)
       do i_s1=1,nbin1
          Fboot(i_s1,iboot)=Fboot(i_s1,iboot)-meanF
       enddo
    enddo
    ! Clip unphysical FE (empty bins) before taking std. err.
    ! Use a ceiling a bit above the largest finite reference barrier.
    meanF=0.d0
    do i_s1=1,nbin1
       if(Fref(i_s1).lt.1.0d8 .and. Fref(i_s1).gt.meanF) meanF=Fref(i_s1)
    enddo
    varF=meanF+10.0d0
    if(varF.lt.20.0d0) varF=20.0d0
    do iboot=1,nboot
       do i_s1=1,nbin1
          if(Fboot(i_s1,iboot).gt.varF) Fboot(i_s1,iboot)=varF
       enddo
    enddo
    do i_s1=1,nbin1
       if(Fref(i_s1).gt.varF) Fref(i_s1)=varF
    enddo
    ! standard deviation over bootstrap samples
    do i_s1=1,nbin1
       meanF=0.d0
       do iboot=1,nboot
          meanF=meanF+Fboot(i_s1,iboot)
       enddo
       meanF=meanF/dfloat(nboot)
       varF=0.d0
       do iboot=1,nboot
          varF=varF+(Fboot(i_s1,iboot)-meanF)**2
       enddo
       if(nboot.gt.1)then
          dF(i_s1)=dsqrt(varF/dfloat(nboot-1))
       else
          dF(i_s1)=0.d0
       end if
    enddo
    deallocate(Fboot,biased_boot,prob_boot,a_boot)
 end if

 call print_pmf(ncv,nbin1,grid0,kt,prob,Fref,dF,nboot,parent)
call WHAM_MPI_Finalize
end program

!***********************************************!

subroutine run_wham(a,umbr_k,umbr_mean,nmax,grid0,biased_prob,umbr_n,nbin1, &
     ncv,kt,prob,cnvg,gleng1_min,gleng1_max,iwrap,toler,max_iter,parent,quiet)
implicit none
integer :: umbr_n,nbin1,ncv,gleng1_min,gleng1_max,iwrap,max_iter,quiet,iter
integer :: nmax(umbr_n)
real*8 :: a(umbr_n),umbr_k(umbr_n),umbr_mean(umbr_n)
real*8 :: grid0(3,ncv),biased_prob(nbin1,umbr_n),prob(nbin1)
real*8 :: kt,cnvg,toler
logical :: parent
 a=1.d0
 iter=0
 scf_loop : do
    iter=iter+1
    call wham_scf(a,umbr_k,umbr_mean,nmax,grid0,biased_prob,umbr_n,nbin1,ncv, &
         kt,prob,cnvg,gleng1_min,gleng1_max,iwrap)
    if (parent .and. quiet.eq.0) &
         write(*,*) 'iteration #', iter, 'convergence =',cnvg
    if(cnvg.lt.toler)then
       if (parent .and. quiet.eq.0) write(*,*)'** convergence achieved **'
       exit scf_loop
    end if
    if(iter.ge.max_iter)then
       if (parent .and. quiet.eq.0) then
          write(*,*)'** WARNING: convergence not reached **'
          write(*,*)'** maximum iteration reached: max_iter =', max_iter, &
               '  final convergence =', cnvg, '  tolerance =', toler
       end if
       exit scf_loop
    end if
 end do scf_loop
end subroutine

!***********************************************!

subroutine gather_prob_to_F(nbin1,grid0,kt,prob,F,parent)
implicit none
integer :: nbin1,i_s1,ncv_dummy
real*8 :: grid0(3,*), kt, prob(nbin1), F(nbin1)
logical :: parent
real*8, allocatable :: prob1(:)
 ncv_dummy=1
 allocate(prob1(nbin1))
 prob1=0.d0
 CALL GlobSumR(prob,prob1,nbin1)
 do i_s1=1,nbin1
    if(prob1(i_s1).gt.0.0d0)then
       F(i_s1)=-kt*dlog(prob1(i_s1))
    else
       F(i_s1)=1.0d16
    end if
 enddo
 call RBcast(F,nbin1)
 deallocate(prob1)
end subroutine

!***********************************************!

subroutine bootstrap_histograms(nbin1,umbr_n,nmax,dcv,hist_count,biased_boot)
! Poisson bootstrap of histogram counts, then renormalize to density.
implicit none
integer :: nbin1,umbr_n,nmax(umbr_n),i_umbr,i_s1
real*8 :: dcv, hist_count(nbin1,umbr_n), biased_boot(nbin1,umbr_n)
real*8 :: cnew, pnorm, lambda
 do i_umbr=1,umbr_n
    pnorm=0.d0
    do i_s1=1,nbin1
       lambda=hist_count(i_s1,i_umbr)
       call rand_poisson(lambda,cnew)
       biased_boot(i_s1,i_umbr)=cnew
       pnorm=pnorm+cnew
    enddo
    if(pnorm.le.0.d0)then
       biased_boot(:,i_umbr)=hist_count(:,i_umbr)
       pnorm=0.d0
       do i_s1=1,nbin1
          pnorm=pnorm+biased_boot(i_s1,i_umbr)
       enddo
    end if
    ! density: counts / (N_eff * dcv); use resampled total as N_eff
    do i_s1=1,nbin1
       biased_boot(i_s1,i_umbr)=biased_boot(i_s1,i_umbr)/(pnorm*dcv)
    enddo
 enddo
end subroutine

!***********************************************!

subroutine rand_poisson(lambda, kout)
! Draw Poisson(lambda). For lambda>30 use Gaussian approximation.
implicit none
real*8 :: lambda, kout, u, L, p, g, x1, x2
integer :: k
 if(lambda.le.0.0d0)then
    kout=0.0d0
    return
 end if
 if(lambda.lt.30.0d0)then
    ! Knuth algorithm
    L=dexp(-lambda)
    k=0
    p=1.0d0
    do
       k=k+1
       call random_number(u)
       p=p*u
       if(p.lt.L) exit
    enddo
    kout=dfloat(k-1)
 else
    ! Box-Muller Gaussian N(lambda, sqrt(lambda))
    call random_number(x1)
    call random_number(x2)
    if(x1.lt.1.0d-16) x1=1.0d-16
    g=dsqrt(-2.0d0*dlog(x1))*dcos(6.283185307179586d0*x2)
    kout=lambda+dsqrt(lambda)*g
    if(kout.lt.0.0d0) kout=0.0d0
 end if
end subroutine

!***********************************************!

 subroutine print_pmf(ncv,nbin1,grid0,kt,prob,Fref,dF,nboot,parent)
!prints 1D pmf to free_energy.dat
 implicit none
 integer:: i_s1, nbin1, ncv, nboot
 real*8 :: prob(nbin1), grid0(3,ncv), Fref(nbin1), dF(nbin1)
 real*8 :: s1, kt
 logical:: parent
 real*8, allocatable :: prob1(:)
 character (len=50):: f2

 allocate(prob1(nbin1))
 prob1=0.0
 CALL GlobSumR(prob,prob1,nbin1)
 if (parent)then
 f2= 'free_energy.dat'
 open( unit =7 , file = f2, status =  'unknown' )
 write(7,'(A)') '# cv   F(cv)   dF(cv)   P(cv)   [kcal/mol]'
 do i_s1=1,nbin1
  s1=DFLOAT(i_s1-1)*grid0(3,1)+grid0(1,1)
  write(7,'(4E16.8)') s1, Fref(i_s1), dF(i_s1), prob1(i_s1)
 enddo
 write(*,*) 'free energy written in ',trim(f2)
 if(nboot.gt.0)then
    write(*,*) ' dF column = bootstrap std. err. from ', nboot, ' samples'
 else
    write(*,*) ' dF column = 0 (set nboot>0 in input to estimate errors)'
 end if
 close(7)
 write(*,*)
 write(*,*) '------------------------------------------------------------'
 write(*,*) ' Example: plot free energy with errors in gnuplot'
 write(*,*) '------------------------------------------------------------'
 write(*,*) ' Interactive:'
 write(*,*) '   gnuplot'
 write(*,*) '   set xlabel "CV"'
 write(*,*) '   set ylabel "F (kcal/mol)"'
 write(*,*) '   plot "free_energy.dat" u 1:2:3 with yerrorlines title "F +/- dF"'
 write(*,*)
 write(*,*) ' Or write plot_fe.gp containing:'
 write(*,*) '   set terminal pngcairo size 800,500'
 write(*,*) '   set output "free_energy.png"'
 write(*,*) '   set xlabel "CV"'
 write(*,*) '   set ylabel "F (kcal/mol)"'
 write(*,*) '   plot "free_energy.dat" u 1:2:3 with yerrorlines title "F +/- dF"'
 write(*,*) ' then run:  gnuplot plot_fe.gp'
 write(*,*) '------------------------------------------------------------'
 end if
 endsubroutine
!**************************************************************************************************!


subroutine wham_scf(a,umbr_k,umbr_mean,nmax,grid0,biased_prob,umbr_n,nbin1,ncv,kt,prob,cnvg,gleng1_min,gleng1_max,iwrap)
!performs 1D wham scf.
! P(s) is renormalized each iteration (int P ds = 1).
implicit none
integer:: i_s1, i_umbr, nbin1, umbr_n, ncv, iwrap
integer:: nmax(umbr_n)
real*8 :: umbr_k(umbr_n), umbr_mean(umbr_n), dum
real*8 :: num, den, dummy_v, dummy_s1, del_s1, kt, cnvg
real*8 :: prob(nbin1), a(umbr_n)
real*8 :: grid0(3,ncv), biased_prob(nbin1,umbr_n)
integer:: gleng1_min, gleng1_max
real*8,allocatable :: dummy_a1(:), dummy_a(:)
real*8 :: pnorm_loc(1), pnorm_glob(1)

allocate(dummy_a(umbr_n))
allocate(dummy_a1(umbr_n))

dummy_a = 0.0d0
prob = 0.0d0

do i_s1 =gleng1_min,gleng1_max
   num = 0.0d0
   den = 0.0d0
   dummy_s1 = grid0(1,1)+dfloat(i_s1-1)*grid0(3,1)

   do i_umbr=1,umbr_n
      del_s1=dummy_s1-umbr_mean(i_umbr)
      if(iwrap.eq.1)then
         if ( del_s1 .gt. 3.14d0 ) del_s1 = del_s1 - 6.28d0
         if ( del_s1 .lt.-3.14d0 ) del_s1 = del_s1 + 6.28d0
      end if
      dummy_v=dexp(-(0.50d0*umbr_k(i_umbr)*del_s1*del_s1/kt))
      num=num+dfloat(nmax(i_umbr))*biased_prob(i_s1,i_umbr)
      den=den+dfloat(nmax(i_umbr))*a(i_umbr)*dummy_v
   enddo

   if(den.gt.0.0d0)then
      prob(i_s1)=num/den
   else
      prob(i_s1)=0.0d0
   end if
   if(prob(i_s1).ne.prob(i_s1)) prob(i_s1)=0.0d0
   if(prob(i_s1)+1.eq.prob(i_s1)) prob(i_s1)=0.0d0
   if(prob(i_s1).lt.0.0d0) prob(i_s1)=0.0d0
enddo

pnorm_loc(1)=0.0d0
do i_s1=gleng1_min,gleng1_max
   pnorm_loc(1)=pnorm_loc(1)+prob(i_s1)*grid0(3,1)
enddo
CALL GlobSumR(pnorm_loc,pnorm_glob,1)
if(pnorm_glob(1).gt.0.0d0)then
   do i_s1=gleng1_min,gleng1_max
      prob(i_s1)=prob(i_s1)/pnorm_glob(1)
   enddo
end if

dum=grid0(3,1)
do i_s1 =gleng1_min,gleng1_max
   dummy_s1 = grid0(1,1)+dfloat(i_s1-1)*grid0(3,1)
   do i_umbr=1,umbr_n
      del_s1=dummy_s1 - umbr_mean(i_umbr)
      if(iwrap.eq.1)then
         if ( del_s1 .ge. 3.14d0 ) del_s1 = 6.28d0 - del_s1
         if ( del_s1 .le.-3.14d0 ) del_s1 = 6.28d0 + del_s1
      end if
      dummy_v=dexp(-(0.50d0*umbr_k(i_umbr)*del_s1*del_s1/kt))
      dummy_a(i_umbr)=dummy_a(i_umbr) +dum*dummy_v*prob(i_s1)
   enddo
enddo

CALL GlobSumR(dummy_a,dummy_a1,umbr_n)
do i_umbr=1,umbr_n
dummy_a(i_umbr)=dummy_a1(i_umbr)
end do

 cnvg = 0.0d0
 do i_umbr=1,umbr_n
 if(dummy_a(i_umbr).le.0.0d0)then
    dummy_a(i_umbr)=1.0d0
 else
    dummy_a(i_umbr) = 1.0d0/dummy_a(i_umbr)
 end if
 cnvg = cnvg + dabs(dlog(dummy_a(i_umbr))-dlog(a(i_umbr)))
 a(i_umbr) = dummy_a(i_umbr)
 enddo
 cnvg = kt*cnvg
 end subroutine
!**************************************************************************************************!



SUBROUTINE DistributeGrids(ncv,grid0,grid,rank,gleng1_min,gleng1_max)
IMPLICIT NONE
INTEGER :: ncv,rank
REAL*8 :: grid0(3,ncv), grid(3,ncv)
INTEGER :: i,ncpu,icpu,ngrids,ngrids_m,ngrids_o
INTEGER :: gleng1_min, gleng1_max

CALL Get_ncpu(ncpu)
CALL Get_cpuid(icpu)

write(*,*)'NCPUS',ncpu
DO i=1,ncv
  grid(1:3,i)=grid0(1:3,i)
END DO
rank=0

if(ncpu.eq.1) then
gleng1_min=1
gleng1_max=nint((grid0(2,1)-grid0(1,1))/grid0(3,1))+1
end if

if(icpu.eq.0)WRITE(*,'(3A12,3A16)') 'CPU','CV', 'GRID SIZE', 'GRID MIN', 'GRID MAX', 'GRID BIN'
CALL Sync_procs
IF(ncpu.GT.1)THEN
  ngrids=nint((grid0(2,1)-grid0(1,1))/grid0(3,1))+1
  write(*,*)'NEW_GRID',ngrids,icpu,ncpu
  ngrids_o=ngrids
  ngrids=ngrids/ncpu
  IF(icpu.eq.ncpu-1)THEN
    ngrids_m=ngrids+mod(ngrids_o,ncpu)
    grid(1,1)=grid(1,1)+DFLOAT(icpu*ngrids)*grid0(3,1)
    grid(2,1)=grid(1,1)+DFLOAT(ngrids_m-1)*grid0(3,1)
  ELSE
    ngrids_m=ngrids
    grid(1,1)=grid(1,1)+DFLOAT(icpu*ngrids)*grid0(3,1)
    grid(2,1)=grid(1,1)+DFLOAT(ngrids-1)*grid0(3,1)
  END IF
  CALL Sync_procs
  WRITE(*,'(3I12,3F16.6)') icpu, 1, ngrids_m, grid(1,1), grid(2,1), grid(3,1)
  rank=ngrids*icpu
  gleng1_min=ngrids*icpu+1
  gleng1_max=ngrids*(icpu+1)
  if(icpu.eq.ncpu-1) gleng1_max=ngrids*(icpu+1)+mod(ngrids_o,ncpu)
END IF
END
!****************************************************************************************!

SUBROUTINE WHAM_MPI_Init()
USE mpi_f08
IMPLICIT NONE
INTEGER :: i_err
call MPI_Init(i_err)
END
!****************************************************************************************!

SUBROUTINE WHAM_MPI_Finalize()
USE mpi_f08
IMPLICIT NONE
INTEGER :: i_err
call MPI_Finalize(i_err)
END
!****************************************************************************************!

SUBROUTINE get_ncpu(ncpu)
USE mpi_f08
IMPLICIT NONE
INTEGER :: ncpu, i_err
call MPI_Comm_size(MPI_COMM_WORLD,ncpu,i_err)
END
!****************************************************************************************!

SUBROUTINE get_cpuid(icpu)
USE mpi_f08
IMPLICIT NONE
INTEGER :: icpu, i_err
call MPI_Comm_rank(MPI_COMM_WORLD,icpu,i_err)
END
!****************************************************************************************!

SUBROUTINE IBcast(myint,leng)
USE mpi_f08
IMPLICIT NONE
INTEGER :: leng, myint(*), i_err
CALL MPI_Bcast(myint,leng,MPI_INTEGER,0,MPI_COMM_WORLD,i_err)
END
!****************************************************************************************!

SUBROUTINE RBcast(myreal,leng)
USE mpi_f08
IMPLICIT NONE
REAL*8 :: myreal(*)
INTEGER :: leng, i_err
CALL MPI_Bcast(myreal,leng,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,i_err)
END
!****************************************************************************************!

SUBROUTINE Sync_Procs
USE mpi_f08
IMPLICIT NONE
INTEGER :: i_err
call MPI_Barrier(MPI_COMM_WORLD,i_err)
END
!****************************************************************************************!

SUBROUTINE Set_Parent(parent)
USE mpi_f08
IMPLICIT NONE
LOGICAL :: parent
INTEGER :: icpu, i_err
parent=.false.
CALL MPI_Comm_rank(MPI_COMM_WORLD,icpu,i_err)
IF(icpu.eq.0)parent=.true.
END
!****************************************************************************************!

SUBROUTINE GlobSumR(myreal_in,myreal_out,leng)
USE mpi_f08
IMPLICIT NONE
REAL*8 :: myreal_in(*), myreal_out(*)
INTEGER :: leng,i_err
CALL MPI_Allreduce(myreal_in,myreal_out,leng,MPI_DOUBLE_PRECISION, &
                   MPI_SUM,MPI_COMM_WORLD,i_err)
END
!****************************************************************************************!

SUBROUTINE GlobSumI(myint_in,myint_out,leng)
USE mpi_f08
IMPLICIT NONE
INTEGER :: myint_in(*), myint_out(*)
INTEGER :: leng,i_err
call MPI_Allreduce(myint_in,myint_out,leng,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,i_err)
END
!****************************************************************************************!
