      Program rpEPR

! This is a program to resolution enhance and/or simulate EPR
! spectra using an interactive command driven interface

#ifdef ifort
      use ifposix
      use ifport
#endif
      use global
      use spectrum_parameters
      use spectrum_parameter_handlers
      use atom_types
      use init_atoms
      use input_output
      use radical_simulation

      IMPLICIT none

      interface
      SUBROUTINE CFFT1I (N, WSAVE, LENSAV, IER)
      INTEGER    N, LENSAV, IER
      REAL       WSAVE(LENSAV)
      end subroutine cfft1i
      end interface

      interface
      SUBROUTINE CFFT1F (N, INC, C, LENC, WSAVE, LENSAV, WORK, LENWRK, IER)
      INTEGER  N, INC, LENC, LENSAV, LENWRK, IER
      COMPLEX  C(LENC)
      REAL     WSAVE(LENSAV)     ,WORK(LENWRK)
      end subroutine cfft1f
      end interface

      interface
      SUBROUTINE CFFT1B (N, INC, C, LENC, WSAVE, LENSAV, WORK, LENWRK, IER)
      INTEGER  N, INC, LENC, LENSAV, LENWRK, IER
      COMPLEX       C(LENC)
      REAL     WSAVE(LENSAV)     ,WORK(LENWRK)
      end subroutine cfft1b
      end interface



      LOGICAL :: flag, header = .true.

      Character (len=40) specnum, num1, num2, num3, num4, num5, num6, &
     & num7, num8, num9, win_func, xrange, yrange
      CHARACTER (len=120) fname, specfile, parmfile,thisarg, refile, &
     & testfile, winstring
      character (len=256) :: prompt, simprompt, re_prompt, spec_prompt
      character (len=256) :: plotstring, rangestring

      CHARACTER(len=20) :: cmdans

      CHARACTER(len=5) :: my_number

      type(spec_parms), pointer :: specpt,scrpt
      type(isotope) :: atom(natoms)
      type(prompt_flag) :: pflag

      logical :: parmflag=.false., specflag=.false., reflag=.false.
      logical :: show_spec=.true., show_re=.true., show_sim=.true.
      logical :: labelflag=.false., simflag=.true., plotflag=.false.

      Real, dimension(10000) :: H, field
      real, dimension(10000) :: refx, refy, enhance_y, apod
      complex, dimension(10000) :: fid
      real :: wsave(33000),work(65536),xstep,tmax, big, small, khigh, klow
      real :: konst, tau, max1, max2, xlow, xhigh, ylow=-2000.0, yhigh=2000.0
      complex scr1(32768),myi
      real, parameter :: pi=3.1415
      real, external :: myrand

      integer :: numpts, ne1, ne2, ilist, I, j, k, l, nnuc(2), &
    &   totnuc, inuc, nspec, ierr, ylabel, kb, cutoff, lensave, &
    &   lenwork, lenspec, ipid


!    "gamma" is intentionally misspelled as "gama" to avoid stepping
!     on the FORTRAN keyword gamma.

      myi=cmplx(0,1)
      lensave = 33000
      lenwork = 655536

      call fill_atoms(atom)
      specpt => parms
      scrpt => parmsnew
      nspec = 1

      ! Change the window name to something unique
      call init_random_seed()
      write(my_number,'(I3)') int(100*myrand())
      winstring='wmctrl -r :ACTIVE: -T rpEPR-'//trim(adjustl(my_number))
      ierr=system(winstring)

      ! This is the command to take back focus from gnuplot
      winstring='wmctrl -R rpEPR-'//trim(adjustl(my_number))
      ierr=system(winstring)

      ! Parse the command line for possible simulation parameters, and experimental files
      ierr = iargc()
      specfile = ''
      if(ierr .gt. 0) then
        if(ierr .gt. 4) ierr=4
        do i = 1, ierr
          call getarg(i,thisarg)
          select case (trim(adjustl(thisarg)))
          case ("-P")
            call getarg(i+1, parmfile)
            parmflag=.true.
          case ("-E")
            call getarg(i+1, specfile)
            specflag=.true.
          case ("-?","-h","--help")
            write(*,*) "rpEPR takes up to 2 command line arguments"
            write(*,*) "-P parameter_file.par specifies a binary, simulation-parameter file from a previous session"
            write(*,*) "-E epr_file specifies a two-column text file containing an EPR spectrum"
            write(*,*) "Example: rpEPR -E my_epr_file -P some_simul_parms.par"
            stop
          end select
        enddo
      endif

        ! Get simulation parameters from file
        if(parmflag)then
         call read_parms(specpt,parmfile,ierr)
          if(ierr == -1) parmflag=.false.
          write(*,*) "Enter a base-name for the simulation files"
          read(*,*) fname
        end if

     ! Get the reference spectrum
     if(specflag)then
       call get_spectrum(specfile,refx,refy,numpts,ierr)
       if(ierr == 0)then
         if(parmflag)then
           specpt%origlow = refx(1)
           specpt%orighigh = refx(numpts)
         endif
       else
         write(*,*) "Error reading spectrum "//trim(adjustl(specfile))
         specflag=.false.
       endif
     endif

! Make a fifo, and fork a process to run gnuplot...intel compiler only
#ifdef ifort
      ! When compiled with gnu you must start gnuplot with a separate shell command
      ! before starting the simulation program
      ! The script rpEPR.sh should launch things properly
      ierr=system('rm plotcmds my_mouse')
      ierr=system('mkfifo plotcmds')
      ierr=system('mkfifo my_mouse')
      call pxffork(ipid,ierr)
      if(ipid == 0)then ! The Child
        ierr=system('gnuplot<plotcmds')
        call pxfexit(ierr)
      else ! The parent
      if(ierr .ne. 0)write(*,*) 'Error forking gnuplot process'
#endif
      open(20,file='plotcmds')

!     Set some default plot properties
      write(20,*) 'set border 1'
      write(20,*) 'unset tics'
      write(20,*) 'set xtics nomirror'
#ifdef ifort
      call pxffflush(20,ierr)
#endif


! Set up for the command loop

      simprompt='g, a, lw, hnu, H0,'// &
     &' intangs, isotope, nucspins1, nucspins2, run, show_sim'
      prompt='enhance, simulate, refspec, label, xrange, yrange, print, write, done'
      re_prompt='show_re'
      spec_prompt='show_spec, scale_spec'
      xrange = '[:]'
      yrange = '[-2000:2000]'

      write(*,*) ' '
      WRITE ( * ,  * ) "Basic commands: "//trim(prompt)
      if(parmflag)write(*,*) "Simulation commands: "//trim(simprompt)
      if(specflag)write(*,*) "Experimental spectrum commands: "//trim(spec_prompt)
      if(reflag)write(*,*) "Resolution Enhancement commands: "//trim(re_prompt)
      WRITE ( * , '(a)' ) 'Enter a command'



!  Begin the simulation

      interact: do 200

! If we already have files called base name, advance specnum until we
! have a new file name
      if(parmflag .and. simflag) then
      DO
        WRITE (specnum, * ) nspec
        specpt%fname = trim (fname) //trim (adjustl (specnum) )
        INQUIRE (file = specpt%fname, exist = flag)
        IF (flag) then
          nspec = nspec + 1
        ELSE
          exit
        ENDIF
      enddo

      WRITE ( * , * ) 'running: '//specpt%fname
      if(.not. specflag) numpts = 0

!      if(specpt%numelectrons == 2)then
!        call calculate(parms,H,field,numpts,1)
!      else
        call radsim(parms,H,field,numpts)
!      end if
      simflag=.false.


! Write the new spectrum to a file
      open(10,file=specpt%fname)
      do i = 1,numpts
        write(10,*) field(i), H(i)
      enddo
      close(10)

! Log the parameters for this spectrum
      call log_parameters(specpt,header=header)
      header = .false.


#ifdef ifort
      call pxffflush(20,ierr)
#endif

      endif ! End of block controlled by parmflag and simflag

! Plot up to 3 traces, eg: simul, experimental, Res. enhancement.
      write(20,*) 'unset label'
      rangestring = ' '//trim(xrange)//' '//trim(yrange)//' "'
      plotstring = 'plot'//rangestring
      plotflag=.false.
      if(parmflag .and. show_sim)then
        write(20,*) trim(adjustl(plotstring))//trim(specpt%fname)//'" with lines'
        plotstring='replot "'
        plotflag=.true.
      endif
      if(specflag .and. show_spec)then
        write(20,*) trim(adjustl(plotstring))//trim(specfile)//'" with lines'
        plotstring='replot "'
        plotflag=.true.
      endif
      if(reflag .and. show_re)then
        write(20,*) trim(adjustl(plotstring))//trim(refile)//'" with lines'
        plotstring='replot "'
        plotflag=.true.
      endif
#ifdef ifort
      call pxffflush(20,ierr)
#endif

      ! Get gnuplot to give us the current xrange
      if(plotflag)then
      open(11,file="my_mouse")
      write(20,*) 'set print "my_mouse"'
      write(20,*) 'print GPVAL_X_MIN, GPVAL_X_MAX'
#ifdef ifort
      call pxffflush(20,ierr)
#endif
      read(11,*) xlow, xhigh
      close(11)
      write(20,*) "set print"
#ifdef ifort
      call pxffflush(20,ierr)
#endif
      endif

! Add labels to the plot if called for
      if(labelflag .and. plotflag)then
        ylabel = 1850
        ! num1 is the position of the left edge
        write(num1,'(f10.1)') xlow
        ! Write the simul parameters
        if(parmflag .and. show_sim)then
          scrpt = specpt
          call parms_to_A(scrpt)
          num1 = adjustl(num1)
          write(num2,'(I8)') ylabel
          num2 = adjustl(num2)
          write(num3, '(f9.4)') scrpt%freq
          num3 = adjustl(num3)
          plotstring = 'set label "Frequency: '//trim(num3)//' GHz" at '&
          //trim(num1)//','//trim(num2)//' font "courier,6"'
          write(20,*) trim(adjustl(plotstring))
          ylabel = ylabel -100
          write(num2,'(I8)') ylabel
          num2 = adjustl(num2)
          write(num3, '(f8.5)') scrpt%g1(1)
          write(num4, '(f8.5)') scrpt%g1(2)
          write(num5, '(f8.5)') scrpt%g1(3)
          plotstring = 'set label "g1: '//trim(num3)//' '//trim(num4)//' '//trim(num5)//'" at '// &
          trim(num1)//','//trim(num2)//' font "courier,6"'
          write(20,*) trim(adjustl(plotstring))
#ifdef ifort
          call pxffflush(20,ierr)
#endif
          ylabel = ylabel -100
          write(num2,'(I8)') ylabel
          num2 = adjustl(num2)
          if(scrpt%numnuclei /= 0)then
            plotstring = 'set label "nuc   e     Ax     Ay     Az  alpha   beta  gamma" at '&
            //trim(num1)//','//trim(num2)//' \'
            write(20,*) trim(adjustl(plotstring))
            write(20,*) ' font "courier,6"'
            ylabel = ylabel -100
            write(num2,'(I8)') ylabel
            num2 = adjustl(num2)
          end if
#ifdef ifort
          call pxffflush(20,ierr)
#endif
          do i = 1,scrpt%numnuclei
            write(num3,'(f7.2)') scrpt%Atens(i,1)
            write(num4,'(f7.2)') scrpt%Atens(i,2)
            write(num5,'(f7.2)') scrpt%Atens(i,3)
            write(num6,'(f7.2)') scrpt%alpha(i)
            write(num7,'(f7.2)') scrpt%beta(i)
            write(num8,'(f7.2)') scrpt%gama(i)
            write(num9,'(I2)') scrpt%which_elec(i)
            plotstring = 'set label "'//scrpt%nucname(i)//trim(num9)&
            //trim(num3)//trim(num4)//trim(num5)//trim(num6)//trim(num7)&
            //trim(num8)//'" at '//trim(num1)//','//trim(num2)//' \'
            write(20,*) trim(adjustl(plotstring))
            write(20,*) ' font "courier,6"'
#ifdef ifort
            call pxffflush(20,ierr)
#endif
            ylabel = ylabel -100
            write(num2,'(I8)') ylabel
            num2 = adjustl(num2)
          end do
        endif
        ! Write the resolution enhancement parms
        if(reflag .and. show_re)then
          write(num2,'(I8)') ylabel
          num2 = adjustl(num2)
          write(num3,'(f5.1)') (-konst*xstep*kb)/(sqrt(3.0)*pi)
          write(num4,'(f6.3)') tmax
          plotstring= "set label 'RE parms: linewidth="&
          //trim(adjustl(num3))//' tmax='//trim(adjustl(num4))&
          //' window='//trim(win_func)//"'  \"
!          ' at '//trim(adjustl(num1))//','//trim(num2)//' font "courier,6"'
          write(20,*) trim(adjustl(plotstring))
!          print*, trim(plotstring)
!          write(20,*) 'font "courier,6"'
          write(20,*) ' at '//trim(adjustl(num1))//','//trim(num2)//' font "courier,6"'
#ifdef ifort
          call pxffflush(20,ierr)
#endif
        endif
      endif
      if(plotflag)then
        write(20,*) 'replot'
#ifdef ifort
        call pxffflush(20,ierr)
#endif
      endif

! The command loop
      comand: do 300

      ! Wait a second for gnuplot to finish, then grab focus back
      !call sleep(1)
      do i = 1,3
      ierr = system(winstring)
      enddo

      WRITE ( *, 99)
   99 FORMAT     (' command?>',$)
      READ ( *, * ) cmdans

      ! Ignore simul options if the parmflag is not set...ie. we don't have an active simul
      if(parmflag)then
      ! This select case handles simulation parameters
      select case (cmdans)
      case ('g')
      WRITE ( * , '(a,3(f8.4),a)' ) 'Current gx, gy, gz:', specpt%g1(1),specpt%g1(2), &
       specpt%g1(3) , ' New g values (3) ?'
      READ ( *, * ) specpt%g1(1), specpt%g1(2), specpt%g1(3)

      case ('lw')
      k = 1
      WRITE ( * , '(a,I1,a,3(f5.1),a)' ) 'Current linewidths:', specpt%lw(k,1) , &
      specpt%lw(k,2) , specpt%lw(k,3),' New linewidths (3) (Gauss) ?'
      READ ( *, * )  specpt%lw(k,1), specpt%lw(k,2), specpt%lw(k,3)

      case ('a')
      k = 1
      IF (nnuc(k) .gt.0) then
      do
      WRITE ( * , 400) 'Enter the number of a nucleus'
        DO l = 1, specpt%numnuclei
        IF ((specpt%which_elec(l) == k) .and. (specpt%geom_type(l) .le. maxatoms)) then
          WRITE(*,'(I2,a,a)') l, ' ', specpt%nucname(l)
        ENDIF
        enddo
        WRITE ( * , * ) 'Number ?'
        READ ( *, * ) inuc
        if(inuc .lt. 1 .or. inuc .gt. nnuc(k))then
          write(*,'(a,I2)') 'number must be between 1 and ',nnuc(k)
          cycle
        end if
        exit
        end do
        if(specpt%geom_type(inuc) == 0)then
          WRITE ( * ,  '(a,3(f7.1),3(f6.1))' ) 'Current hyperfine tensor and euler angles: ', &
          specpt%Atens(inuc, 1:3), specpt%alpha(inuc), specpt%beta(inuc), &
          specpt%gama(inuc)
        else if(specpt%geom_type(inuc) == 1)then
          write(*,'(a,f6.1,5(f7.1))') 'Current dihedral angle(phi), A0, A1, delta_Ax, delta_Ay, delta_Az', &
          specpt%Atens(inuc, 1:3), specpt%alpha(inuc), specpt%beta(inuc), &
          specpt%gama(inuc)
        end if
          WRITE ( * ,  '(a)' ) 'Enter new parameters for: '//&
    &      trim(specpt%nucname (inuc))//'(Gauss),(Degrees)'
          READ ( *, * ) specpt%Atens(inuc,1), specpt%Atens(inuc,2), specpt%Atens(inuc,3), &
           specpt%alpha(inuc), specpt%beta(inuc), specpt%gama(inuc)
      ELSE
        WRITE ( * , '(a,I2)' ) 'No nuclei coupled to electron'
      ENDIF

      case ('hnu')
      WRITE ( * ,  '(a,f7.4,a)' ) 'Current microwave frequency', specpt%freq,&
    &  ' new frequency(GHz) ?'
      READ ( *, * )  specpt%freq

      case ('H0', 'h0')
      WRITE ( * ,  '(a,2(f7.0),a)' ) 'Current field range',&
    &  specpt%origlow, specpt%orighigh, ' new field range(Gauss) ?'
      READ ( *, * ) specpt%origlow, specpt%orighigh

      case ('intangs')
      WRITE ( * , 402) 'Current integration angles', specpt%LT,&
    &  specpt%LP, 'new angles(degrees) ?'
      READ ( *, * ) specpt%LT, specpt%LP

      case ('run')
      nspec = nspec + 1
      WRITE (specnum, * ) nspec
      specpt%fname = trim (fname) //trim (adjustl (specnum) )
!     write(*,*) 'running: '//specpt%fname
      simflag=.true.
      exit comand

      case ('isotope')
      write(*,*) 'Current list of nuclei'
      do i = 1, specpt%numnuclei
        write(*,'(I2,a)') i, ' '//trim(specpt%nucname(i))
      end do
      write(*,'(a,I2,a)') 'Enter a number between 1 and ',specpt%numnuclei,' to choose a nucleus'
      read(*,*) inuc
      write(*,*) 'Available isotopes'
      l = specpt%what_isotope(inuc)
      do i = 1, natoms
        if(atom(l)%atomic_number == atom(i)%atomic_number) &
        write(*,'(I2,a)') i, ' '//trim(atom(i)%isoname)
      end do
      write(*,*) 'Enter the number of an isotope'
      read(*,*) i
      if(atom(l)%atomic_number /= atom(i)%atomic_number)then
        write(*,*) 'Not an isotope pair I recognise'
        cycle comand
      end if
      specpt%what_isotope(inuc) = i
      specpt%nucname(inuc) = atom(i)%isoname

      select case (specpt%geom_type(inuc))
      case (0) ! Ordinary nucleus
        specpt%Atens(inuc,1:3) = specpt%Atens(inuc,1:3)*(atom(i)%isognuc/atom(l)%isognuc)

      case (1) ! Primary beta nucleus
        specpt%Atens(inuc,2:3) = specpt%Atens(inuc,2:3)*(atom(i)%isognuc/atom(l)%isognuc)
        specpt%alpha(inuc) = specpt%alpha(inuc)*(atom(i)%isognuc/atom(l)%isognuc)
        specpt%beta(inuc) = specpt%beta(inuc)*(atom(i)%isognuc/atom(l)%isognuc)
        specpt%gama(inuc) = specpt%gama(inuc)*(atom(i)%isognuc/atom(l)%isognuc)

      end select

      case ('nucspins')
      header = .true.
      k = 1
      j=2
!      totnuc = specpt%numnuclei - nnuc(k)
!      scrpt = specpt
!      ne1 = nnuc(1)
!      ne2 = nnuc(2)
!      ilist = 0
!      if(k == 2)ilist = ne1
!      nnuc(k) = 0
      nnuc = 0
      do i=1,specpt%numnuclei
        if(specpt%which_elec(i) == 1) then
          nnuc(1) = nnuc(1) + 1
        else
          nnuc(2) = nnuc(2) + 1
        end if
      end do
      l=1
      do i = 1, specpt%numnuclei
        if(specpt%which_elec(i) == j) then
          scrpt%Atens(l,:) = specpt%Atens(i,:)
          scrpt%alpha(l) = specpt%alpha(i)
          scrpt%beta(l) = specpt%beta(i)
          scrpt%gama(l) = specpt%gama(i)
          scrpt%which_elec(l) = specpt%which_elec(i)
          scrpt%nucname(l) = specpt%nucname(i)
          scrpt%what_isotope(l) = specpt%what_isotope(i)
          scrpt%geom_type(l) = specpt%geom_type(i)
          l = l + 1
        end if
      end do
      specpt%geom_type = 0
      specpt%numnuclei = 0
      WRITE ( * , 400) 'Enter new nucleus set'
!        call spec_from_terminal(specpt,pflag)
      call nucspins(specpt,k)
      write (*,400) 'Enter hyperfine values for nuclei'
      call Avalues(specpt,k)
      do i = 1,nnuc(j)
        specpt%Atens(i+specpt%numnuclei,:) = scrpt%Atens(i,:)
        specpt%alpha(i+specpt%numnuclei) = scrpt%alpha(i)
        specpt%beta(i+specpt%numnuclei) = scrpt%beta(i)
        specpt%gama(i+specpt%numnuclei) = scrpt%gama(i)
        specpt%which_elec(i+specpt%numnuclei) = scrpt%which_elec(i)
        specpt%nucname(i+specpt%numnuclei) = scrpt%nucname(i)
        specpt%what_isotope(i+specpt%numnuclei) = scrpt%what_isotope(i)
        specpt%geom_type(i+specpt%numnuclei) = scrpt%geom_type(i)
      end do
      specpt%numnuclei = specpt%numnuclei + nnuc(j)
!      nnuc(k) = specpt%numnuclei
!      if(k == 1)then
!        specpt%numnuclei = specpt%numnuclei + ne2
!      else
!        specpt%numnuclei = specpt%numnuclei + ne1
!      end if
!      if(k == 1 )then !Copy the e2 nuclei after the new e1 nuclei
!        do i = 1,ne2
!          specpt%Atens(ilist+i,:) = scrpt%Atens(ne1+i,:)
!          specpt%alpha(ilist+i) = scrpt%alpha(ne1+i)
!          specpt%beta(ilist+i) = scrpt%beta(ne1+i)
!          specpt%gama(ilist+i) = scrpt%gama(ne1+i)
!        end do
!      end if
       end select
       endif

      ! A separate select case to handle non simulation parameter related commands
      select case (cmdans)
      case ("simulate")
       WRITE ( * , * ) 'Enter a base filename for the simulations'
       READ ( *, * ) fname
!       specpt%fname = trim(fname)
       pflag%fname = .false.

       !Dummy angles so spec_from_terminal only prompts for them in the 2 electron case
       specpt%LT = 10
       specpt%LP = 10
       specpt%numnuclei = 0
       pflag%intang = .false.
       pflag%field = .not. specflag
       call spec_from_terminal(specpt,pflag)
       !radsim will determine angle counts itself if given zeroes
       if(specpt%numelectrons == 1)then
         specpt%LT = 0
         specpt%LP = 0
       end if
       parmflag=.true.
       simflag=.true.

      ! Count how many nuclei are coupled to each electron
       nnuc = 0
       do i = 1,specpt%numnuclei
         nnuc(specpt%which_elec(i)) = nnuc(specpt%which_elec(i)) + 1
       end do
       exit comand

      case ('refspec', 'exp_spec')
   88 WRITE ( * , * ) 'Current reference spectrum "'//trim (specfile) //&
      '" new file name?, or "none" to skip'
      testfile = ''
      READ ( *, * ) testfile
      IF (trim (testfile) .ne.'none') then
        INQUIRE (file = trim (testfile), exist = flag)
        IF (.not.flag) then
          WRITE ( * , * ) trim (testfile) //' not found'
          GOTO 88
        ENDIf
        call get_spectrum(testfile,refx,refy,numpts,ierr)
        if(ierr == 0)then
          specfile=''
          specfile=testfile
          specpt%origlow = refx(1)
          specpt%orighigh = refx(numpts)
          write(*,'(a,2(f7.0))') 'New scan range: ',specpt%origlow,specpt%orighigh
          xrange = '[:]'
          xlow=specpt%origlow
          xhigh=specpt%orighigh
          specflag=.true.
          exit comand
        end if
      ENDIF

      case ('print')
        write(20,*) 'set term postscript color'
        write(20,*) 'set output "'//trim(specpt%fname)//'.ps"'
        write(20,*) 'set xlabel "Gauss"'
        write(20,*) 'replot'
#ifdef ifort
        call pxffflush(20,ierr)
#endif
        do i=1,5 ! give gnuplot up to 5 seconds to create the postscript file
          call sleep(1)
          inquire(file=trim(specpt%fname)//'.ps', exist=flag)
          if(flag)then
            ierr=system('lpr '//trim(specpt%fname)//'.ps')
            exit
          end if
        end do
        if(i == 6) write(*,*) 'Error accessing postscript file for printing'
        write(20,*) 'set term x11'
        write(20,*) 'set output'
#ifdef ifort
        call pxffflush(20,ierr)
#endif

      case ('write')
        ! Write an encapsulated postscript file
        if((parmflag .and. show_sim) .or. (specflag .and. show_spec)&
    &     .or. (reflag .and. show_re))then
          testfile = ''
          if(parmflag)then
            write(testfile,'(a)') trim(specpt%fname)//'.eps'
          else
            write(testfile,'(a)') trim(specfile)//'.eps'
          endif
          write(20,*) 'set term postscript color eps solid size 3.25, 3.25'
          write(20,*) 'set output "'//trim(adjustl(testfile))//'"'
          write(20,*) 'set xlabel "Gauss"'
          write(20,*) 'replot'
#ifdef ifort
          call pxffflush(20,ierr)
#endif
          write(20,*) 'set term wxt'
          write(20,*) 'set output'
#ifdef ifort
          call pxffflush(20,ierr)
#endif
          write(*,*) ' '
          write(*,*) 'Wrote figure to file: '//trim(testfile)
        else
          write(*,*) "No spectra shown, no figure written"
          write(*,*) "Toggle display on using: show_... and try again"
        end if

        ! Write a binary parameter file
        if(parmflag)then
         call write_parms(specpt,ierr)
        endif

      case ('label')
        labelflag=.not. labelflag
        exit comand

      case ('done')
        exit interact

      ! A longish section to perform a RE
      case ('enhance')
        ! Get a spectrum if we don't have one
        if(.not. specflag)then
          401 WRITE ( * ,  * ) 'Experimental spectrum, file name?'
          READ ( *, * ) specfile
          INQUIRE (file = trim (specfile), exist = flag)
          IF (.not.flag) then
            WRITE ( * , * ) trim (specfile) //' not found'
            GOTO 401
          ENDIF
          call get_spectrum(specfile,refx,refy,numpts,ierr)
          specflag=.true.
          xrange = '[:]'
          specpt%origlow = refx(1)
          specpt%orighigh = refx(numpts)
!          xlow=refx(1)
!          xhigh=refx(numpts)
        endif

        ! Get a filename for the RE
        write(*,*) "Resolution enhanced spectrum, file name?"
        read(*,*) refile

        !Begin the RE here
        xstep = (refx(numpts)-refx(1))/(numpts-1)
        ! Zerofill the spectrum
        kb = 1024
        do i = 1,5
          if(kb .gt. numpts)exit
          kb = kb * 2
        end do
        lenspec = kb
        scr1(1:kb)=cmplx(0.,0.)
        do i=1,numpts
          scr1(i+(kb-numpts)/2) = cmplx(refy(i),0)
        end do
        do i=1,kb/4
          scr1(i)=scr1(1+kb/4)
          scr1(i+numpts+kb/4)=scr1(numpts+kb/4)
        enddo
        ! Generate the FID
        call cfft1i(kb,wsave,lensave,ierr)
        call cfft1b(kb,1,scr1,lenspec,wsave,lensave,work,lenwork,ierr)
        open(11,file="abs(my_FID)",status="replace")
        scr1(1) = cmplx(0,0) ! This is a cheat to get rid of DC bias in the spectrum, makes the FID look better too!
        ! Display the FID in gnuplot, and get the user to set the cutoff using the mouse
        do i=1,kb/2
          write(11,*) i/(xstep*kb),abs(scr1(i))
        end do
        close(11)
        open(11,file="my_mouse")
        write(20,*) 'set print "my_mouse"'
        write(20,*) "set mouse"
        write(20,*) 'set title "Click the mouse at the point you want to cut off the FID"'
        write(20,*) 'plot "abs(my_FID)" with lines'
        write(20,*) "pause mouse any"
        write(20,*) "print MOUSE_X"
#ifdef ifort
        call pxffflush(20,ierr)
#endif
        read(11,*) tmax
        close(11)
        write(20,*) "unset mouse"
        write(20,*) "set print"
        write(20,*) "unset title"
#ifdef ifort
        call pxffflush(20,ierr)
#endif
        cutoff=1+(tmax*xstep*kb)

        ! Use a heuristic to generate an initial estimate of the lw (konst)
        big=0.
        small=huge(small)
        do i=1,cutoff/2
          if(big .lt. abs(scr1(i)))big=abs(scr1(i))
        enddo
        do i=1+cutoff/2,cutoff
          if(small .gt. abs(scr1(i)))small=abs(scr1(i))
        enddo
        if(small == 0) small = .01
        konst=log(small/big)/float(cutoff)
        khigh=konst
        klow=konst/10.
        do
          do i = 1, cutoff
            fid(i)=scr1(i)*exp(-konst*float(i))
          enddo
          max1=0.
          max2=0.
          do i=1,cutoff/2
            if(abs(fid(i)) .gt. max1) max1= abs(fid(i))
          end do
          do i=1+cutoff/2,cutoff
            if(abs(fid(i)) .gt. max2) max2= abs(fid(i))
          end do
          if(max1 .lt. max2)then
            khigh = konst
          else
            klow = konst
          endif
          konst = klow + (khigh-klow)/2.
          if(abs((khigh-klow)/min(klow,khigh)) .lt. .01)exit
        enddo
        konst = klow

        ! Refine the resolution enhancement interactively
        cmdans="go" ! A dummy cmdans causes us to fall through to default on the first trip
        win_func= "cosine"
        do !The interactive loop for the re

          ! Prompt user for new RE parameters
          !Wait a second for gnuplot to finish, then grab focus back
            !call sleep(1)
            do i = 1,3
            ierr = system(winstring)
            enddo
          if(trim(cmdans) == "?")then
            write(*,'("Res Enhance cmd?>",$)')
            read(*,*) cmdans
            select case (cmdans)
            case ("tm")
              write(*,'(a,f7.3,a)') "Current tm is ",tmax, " New tmax?"
              read(*,*) tmax
              cutoff=1+(tmax*xstep*kb)

            case ("lw")
              write(*,'(a,f6.1,a)') "Current linewidth is ",(-konst*xstep*kb)/(sqrt(3.0)*pi)," New linewidth?"
              read(*,*) konst
              konst=-konst*(sqrt(3.0)*pi)/(xstep*kb)

            case ("win")
              write(*,*) 'Current window function is "'//trim(adjustl(win_func))//'" New window function?'
              write(*,*) "The available windowing functions in order from conservative to aggressive  are:"
              write(*,*) "1 squared bartlett"
              write(*,*) "2 cosine"
              write(*,*) "3 gaussian"
              write(*,*) "4 raised cosine"
              write(*,*) "5 triangle"
              write(*,*) "6 welch"
              write(*,*) "7 square"
              write(*,*) "Enter the number of a windowing functions"
              read(*,*) j
              select case (j)
              case (1)
                win_func="squared bartlett"
              case (2)
                win_func="cosine"
              case (3)
                win_func="gaussian"
              case (4)
                win_func="raised cosine"
              case (5)
                win_func="triangle"
              case (6)
                win_func="welch"
              case (7)
                win_func="square"
              end select

            case ("done")
              exit

            case default
              write(*,*) ' '
              write(*,*) "We are adjusting your RE parameters just now"
              write(*,*) "The valid commands are: tm, lw, win, or done"
            end select
          endif

          apod(1:cutoff)=1.0
          apod(cutoff+1:kb)=0.

          do ! This will loop if a bad value of win_func is given

          ! Case construct to generate the chosen apodization function
          select case (trim(win_func))
          case ("square")
            exit

          case ("welch")
            do i = 1, cutoff
              tau = float(i-1)/(2.*float(cutoff-1))
              apod(i)=apod(i)*(1.-(2.*tau)**2)
            enddo
            exit

          case ("triangle")
            do i = 1, cutoff
              tau = float(i-1)/(2.*float(cutoff-1))
              apod(i)=apod(i)*(1-2*tau)
            enddo
            exit

          case ("raised cosine")
            do i = 1, cutoff
              tau = float(i-1)/(2.*float(cutoff-1))
              apod(i)=apod(i)*(0.54+0.46*cos(2*pi*tau))
            enddo
            exit

          case ("gaussian")
            do i = 1, cutoff
              tau = float(i-1)/(2.*float(cutoff-1))
              apod(i)=apod(i)*(exp(-.5*((5.0*tau)**2)))
            enddo
            exit

          case ("cosine")
            do i = 1, cutoff
              tau = float(i-1)/(2.*float(cutoff-1))
              apod(i)=apod(i)*(0.5+0.5*cos(2*pi*tau))
            enddo
            exit

          case ("squared bartlett")
            do i = 1, cutoff
              tau = float(i-1)/(2.*float(cutoff-1))
              apod(i)=apod(i)*(1-2*tau)**2
            enddo
            exit

          case default
            write(*,*) trim(win_func)//" is not one I recognize, try again"
            write(*,*) "The available windowing functions in order from conservative to aggressive  are:"
            write(*,*) "1 squared bartlett"
            write(*,*) "2 cosine"
            write(*,*) "3 gaussian"
            write(*,*) "4 raised cosine"
            write(*,*) "5 triangle"
            write(*,*) "6 welch"
            write(*,*) "7 square"
            write(*,'(a,$)') "Enter the number of a windowing function:"
            read(*,*) j
            select case (j)
            case (1)
              win_func="squared bartlett"
            case (2)
              win_func="cosine"
            case (3)
              win_func="gaussian"
            case (4)
              win_func="raised cosine"
            case (5)
              win_func="triangle"
            case (6)
              win_func="welch"
            case (7)
              win_func="square"
            end select
          end select
          enddo

          ! Multiply the apodization by the rising exponential
          do i = 1, cutoff
            apod(i)=exp(-konst*float(i))*apod(i)
            apod(kb-i+1) = apod(i)
          enddo

          ! Apply the apodization, and do the back transform
          do i = 1, kb
            fid(i)=scr1(i)*apod(i)
          enddo
          ! Transform back to the spectrum
          call cfft1f(kb,1,fid,lenspec,wsave,lensave,work,lenwork,ierr)

          ! Copy the real component
          do i = 1,numpts
            enhance_y(i) = numpts*real(fid(i+(kb-numpts)/2))
          end do

          !Scale to -2000 to 2000
          call scale_spec(enhance_y,numpts)


          ! Write the new RE to a file
          open(11,file=trim(adjustl(refile)),status="replace")
          do i = 1, numpts
            write(11,*) refx(i), enhance_y(i)
          enddo
          close(11)

          ! Plot the new RE with gnuplot
          rangestring = ' '//trim(xrange)//' '//trim(yrange)//' "'
          plotstring = 'plot'//trim(rangestring)
          write(20,*) trim(adjustl(plotstring))//trim(adjustl(refile))//'" with lines'
#ifdef ifort
          call pxffflush(20,ierr)
#endif
          cmdans="?"
        enddo
        reflag=.true.
        exit comand
        ! The RE section ends here

      case ('scale_spec')
        if(specflag)then
          call scale_spec(refy,numpts)
        else
          write(*,*) "No experimental spectrum loaded"
        endif
        ! Write the scaled spectrum to a new file
        specfile=trim(specfile)//'_scaled'
        open(11,file=trim(specfile),status='replace')
        do i = 1, numpts
          write(11,*) refx(i), refy(i)
        enddo
        close(11)
        exit comand

      case ('show_spec')
        show_spec = .not. show_spec
        exit comand

      case ('show_re')
        show_re = .not. show_re
        exit comand

      case ('show_sim')
        show_sim = .not. show_sim
        exit comand

      case ('xrange')
        write(*,*) "Current xrange "//trim(xrange)
        write(*,*) "Enter new xrange"
        read(*,*) xlow, xhigh
        write(num1,'(f10.1)') xlow
        write(num2,'(f10.1)') xhigh
        xrange = '['//trim(num1)//':'//trim(num2)//']'
        exit comand

      case ('yrange')
        write(*,*) "Current yrange "//trim(yrange)
        write(*,*) "Enter new yrange"
        read(*,*) ylow, yhigh
        write(num1,'(f10.1)') ylow
        write(num2,'(f10.1)') yhigh
        yrange = '['//trim(num1)//':'//trim(num2)//']'
        exit comand

      case default
        write(*,*) ' '
        WRITE (*,*) "Basic commands: "//trim(prompt)
        if(parmflag)write(*,*) "Simulation commands: "//trim(simprompt)
        if(specflag)write(*,*) "Experimental spectrum commands: "//trim(spec_prompt)
        if(reflag)write(*,*) "Resolution Enhancement commands: "//trim(re_prompt)
        write(*,*) ' '
        WRITE ( * , '(a)' ) 'Enter a command'

      end select
     !End of the comand input loop
  300 end do comand

      !End of the interactive loop
  200 end do interact

#ifdef ifort
      end if ! End fork if-block
#endif
      write(20,*) 'quit'
#ifdef ifort
      call pxffflush(20,ierr)
#endif
      close(20)
#ifdef ifort
      ierr=system('rm plotcmds')
      ierr=system('rm my_mouse')
      ierr=system('rm abs\(my_FID\)')
#endif

  400 FORMAT(a,I2)
  402 FORMAT(a,I4,I4)
      STOP
      END


      subroutine scale_spec(spec,n)
      ! This routine scales the argument spectrum to -2000 to 2000
      implicit none
      ! Arguments
      real, dimension(10000) :: spec
      integer :: n
      ! Locals
      real :: big, small, tau
      integer :: i

      small=huge(small)
      big=0.
      do i=1,n
        if(spec(i) .gt. big)big=spec(i)
        if(spec(i) .lt. small)small=spec(i)
      enddo
      tau=(big-small)/4000.
      do i = 1, n
        spec(i)=spec(i)/tau
        spec(i)=spec(i)-(big/tau - 2000.)
      enddo
      end subroutine scale_spec

         ! Random number wrapper
real FUNCTION myrand()
  real a
  CALL RANDOM_NUMBER(a)
  myrand=a
END FUNCTION myrand


! Get random seeds from /dev/urandom, linux specific
SUBROUTINE init_random_seed()
   INTEGER :: n
   INTEGER, DIMENSION(:), ALLOCATABLE :: seed
   integer, parameter :: randin=15

   CALL RANDOM_SEED(size = n)
   ALLOCATE(seed(n))

! nothing like a little ugly preprocessing to get around compiler incompatibilities
#ifdef ifort
   open(randin,file='/dev/urandom',form='binary') ! This works with ifort, other compilers????
#endif
#ifdef gnu
   open(randin,file='/dev/urandom',access='stream',form='unformatted') ! gfortran version
#endif
   read(randin) seed(1:n)
   close(randin)

   CALL RANDOM_SEED(PUT = seed)

   DEALLOCATE(seed)
END SUBROUTINE init_random_seed

