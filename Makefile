include make_header

obj/atoms_pt.o: source/atoms_pt.f90 
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c -o obj/atoms_pt.o source/atoms_pt.f90

obj/c1fgkb.o: fftpack5/c1fgkb.f 
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c -o obj/c1fgkb.o fftpack5/c1fgkb.f

obj/gaussian.o: source/gaussian.f90 
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c -o obj/gaussian.o source/gaussian.f90

obj/xerfft.o: fftpack5/xerfft.f 
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c -o obj/xerfft.o fftpack5/xerfft.f

obj/c1f5kf.o: fftpack5/c1f5kf.f 
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c -o obj/c1f5kf.o fftpack5/c1f5kf.f

obj/tables.o: fftpack5/tables.f 
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c -o obj/tables.o fftpack5/tables.f

obj/c1f2kf.o: fftpack5/c1f2kf.f 
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c -o obj/c1f2kf.o fftpack5/c1f2kf.f

obj/c1f5kb.o: fftpack5/c1f5kb.f 
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c -o obj/c1f5kb.o fftpack5/c1f5kb.f

obj/c1f4kb.o: fftpack5/c1f4kb.f 
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c -o obj/c1f4kb.o fftpack5/c1f4kb.f

obj/c1f3kb.o: fftpack5/c1f3kb.f 
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c -o obj/c1f3kb.o fftpack5/c1f3kb.f

obj/c1f2kb.o: fftpack5/c1f2kb.f 
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c -o obj/c1f2kb.o fftpack5/c1f2kb.f

obj/c1fgkf.o: fftpack5/c1fgkf.f 
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c -o obj/c1fgkf.o fftpack5/c1fgkf.f

obj/constants_pt.o: source/constants_pt.f90 
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c -o obj/constants_pt.o source/constants_pt.f90

obj/factor.o: fftpack5/factor.f 
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c -o obj/factor.o fftpack5/factor.f

obj/c1f3kf.o: fftpack5/c1f3kf.f 
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c -o obj/c1f3kf.o fftpack5/c1f3kf.f

obj/c1f4kf.o: fftpack5/c1f4kf.f 
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c -o obj/c1f4kf.o fftpack5/c1f4kf.f

obj/anneal_globals.o: source/anneal_globals.f90 obj/spectrum_parameters_pt.o obj/atoms_pt.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c -o obj/anneal_globals.o source/anneal_globals.f90

obj/cfft1f.o: fftpack5/cfft1f.f obj/xerfft.o obj/c1fm1f.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c -o obj/cfft1f.o fftpack5/cfft1f.f

obj/cfft1i.o: fftpack5/cfft1i.f obj/xerfft.o obj/mcfti1.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c -o obj/cfft1i.o fftpack5/cfft1i.f

obj/spectrum_parameters_pt.o: source/spectrum_parameters_pt.f90 obj/atoms_pt.o obj/constants_pt.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c -o obj/spectrum_parameters_pt.o source/spectrum_parameters_pt.f90

obj/mcfti1.o: fftpack5/mcfti1.f obj/factor.o obj/tables.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c -o obj/mcfti1.o fftpack5/mcfti1.f

obj/cfft1b.o: fftpack5/cfft1b.f obj/xerfft.o obj/c1fm1b.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c -o obj/cfft1b.o fftpack5/cfft1b.f

obj/input_output_routines_pt.o: source/input_output_routines_pt.f90 obj/anneal_globals.o obj/spectrum_parameters_pt.o obj/atoms_pt.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c -o obj/input_output_routines_pt.o source/input_output_routines_pt.f90

obj/convolute_fftpack.o: source/convolute_fftpack.f obj/cfft1i.o obj/cfft1f.o obj/cfft1b.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c -o obj/convolute_fftpack.o source/convolute_fftpack.f

obj/radsim_pt.o: source/radsim_pt.f90 obj/spectrum_parameters_pt.o obj/atoms_pt.o obj/gaussian.o obj/constants_pt.o obj/convolute_fftpack.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c -o obj/radsim_pt.o source/radsim_pt.f90

obj/c1fm1b.o: fftpack5/c1fm1b.f obj/c1f2kb.o obj/c1f3kb.o obj/c1f4kb.o obj/c1f5kb.o obj/c1fgkb.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c -o obj/c1fm1b.o fftpack5/c1fm1b.f

obj/c1fm1f.o: fftpack5/c1fm1f.f obj/c1f2kf.o obj/c1f3kf.o obj/c1f4kf.o obj/c1f5kf.o obj/c1fgkf.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c -o obj/c1fm1f.o fftpack5/c1fm1f.f

obj/rpEPR.o: source/rpEPR.F90 obj/anneal_globals.o obj/spectrum_parameters_pt.o obj/atoms_pt.o obj/input_output_routines_pt.o obj/radsim_pt.o obj/cfft1i.o obj/cfft1f.o obj/cfft1b.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c -o obj/rpEPR.o source/rpEPR.F90

OBJ_rpEPR = obj/rpEPR.o obj/cfft1b.o obj/c1fm1b.o obj/c1fgkb.o obj/c1f5kb.o obj/c1f4kb.o obj/c1f3kb.o obj/c1f2kb.o obj/xerfft.o obj/cfft1f.o obj/c1fm1f.o obj/c1fgkf.o obj/c1f5kf.o obj/c1f4kf.o obj/c1f3kf.o obj/c1f2kf.o obj/cfft1i.o obj/mcfti1.o obj/tables.o obj/factor.o obj/radsim_pt.o obj/convolute_fftpack.o obj/constants_pt.o obj/gaussian.o obj/atoms_pt.o obj/spectrum_parameters_pt.o obj/input_output_routines_pt.o obj/anneal_globals.o 

SRC = source/input_output_routines_pt.f90 fftpack5/c1f2kb.f fftpack5/c1f2kf.f fftpack5/c1f3kb.f fftpack5/c1f3kf.f fftpack5/c1f4kb.f fftpack5/c1f4kf.f fftpack5/c1f5kb.f fftpack5/c1f5kf.f fftpack5/c1fgkb.f fftpack5/c1fgkf.f fftpack5/c1fm1b.f fftpack5/c1fm1f.f fftpack5/cfft1b.f fftpack5/cfft1f.f fftpack5/cfft1i.f fftpack5/cfft2b.f fftpack5/cfft2f.f fftpack5/cfft2i.f fftpack5/cfftmb.f fftpack5/cfftmf.f fftpack5/cfftmi.f fftpack5/cmf2kb.f fftpack5/cmf2kf.f fftpack5/cmf3kb.f fftpack5/cmf3kf.f fftpack5/cmf4kb.f fftpack5/cmf4kf.f fftpack5/cmf5kb.f fftpack5/cmf5kf.f fftpack5/cmfgkb.f fftpack5/cmfgkf.f fftpack5/cmfm1b.f fftpack5/cmfm1f.f fftpack5/cosq1b.f fftpack5/cosq1f.f fftpack5/cosq1i.f fftpack5/cosqb1.f fftpack5/cosqf1.f fftpack5/cosqmb.f fftpack5/cosqmf.f fftpack5/cosqmi.f fftpack5/cost1b.f fftpack5/cost1f.f fftpack5/cost1i.f fftpack5/costb1.f fftpack5/costf1.f fftpack5/costmb.f fftpack5/costmf.f fftpack5/costmi.f fftpack5/factor.f fftpack5/mcfti1.f fftpack5/mcsqb1.f fftpack5/mcsqf1.f fftpack5/mcstb1.f fftpack5/mcstf1.f fftpack5/mradb2.f fftpack5/mradb3.f fftpack5/mradb4.f fftpack5/mradb5.f fftpack5/mradbg.f fftpack5/mradf2.f fftpack5/mradf3.f fftpack5/mradf4.f fftpack5/mradf5.f fftpack5/mradfg.f fftpack5/mrftb1.f fftpack5/mrftf1.f fftpack5/mrfti1.f fftpack5/msntb1.f fftpack5/msntf1.f fftpack5/r1f2kb.f fftpack5/r1f2kf.f fftpack5/r1f3kb.f fftpack5/r1f3kf.f fftpack5/r1f4kb.f fftpack5/r1f4kf.f fftpack5/r1f5kb.f fftpack5/r1f5kf.f fftpack5/r1fgkb.f fftpack5/r1fgkf.f fftpack5/rfft1b.f fftpack5/rfft1f.f fftpack5/rfft1i.f fftpack5/rfft2b.f fftpack5/rfft2f.f fftpack5/rfft2i.f fftpack5/rfftb1.f fftpack5/rfftf1.f fftpack5/rffti1.f fftpack5/rfftmb.f fftpack5/rfftmf.f fftpack5/rfftmi.f fftpack5/sinq1b.f fftpack5/sinq1f.f fftpack5/sinq1i.f fftpack5/sinqmb.f fftpack5/sinqmf.f fftpack5/sinqmi.f fftpack5/sint1b.f fftpack5/sint1f.f fftpack5/sint1i.f fftpack5/sintb1.f fftpack5/sintf1.f fftpack5/sintmb.f fftpack5/sintmf.f fftpack5/sintmi.f fftpack5/tables.f fftpack5/xercon.f fftpack5/xerfft.f source/atoms_pt.f90 source/convolute_fftpack.f source/gaussian.f90 source/anneal_globals.f90 source/rpEPR.F90 source/spectrum_parameters_pt.f90 source/constants_pt.f90 source/radsim_pt.f90 

bin/rpEPR: $(OBJ_rpEPR)
	$(LD) $(OBJ_rpEPR) -o bin/rpEPR $(LDFLAGS)

clean:
	rm obj/* bin/* mod/*

make:
	@mv Makefile Makefile.bak
	@./fort_deps_ctags.pl>Makefile

TAGS: $(SRC)
	@etags $(SRC)

tags: $(SRC)
	@ctags $(SRC)

