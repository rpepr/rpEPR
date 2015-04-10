module global ! Define some variables needed thoughout the annealing program
 use spectrum_parameters
 use atom_types, only: maxatoms
 public

type prompt_flag
  logical :: freq = .true.
  logical :: field = .true.
  logical :: ne = .true.
  logical :: fname = .true.
  logical :: zfs = .true.
  logical :: intang = .true.
  logical :: gtens(2) = .true.
  logical :: lw(2) = .true.
  logical :: nucspins(2) = .true.
  logical :: Atens(maxatoms) = .true.
end type prompt_flag

  integer, parameter :: maxparms = maxatoms*6 + 15 ! 15 plus 6 times maxatoms
  integer, parameter :: maxiso = maxatoms + 4 ! 4 plus maxatoms
  type(real8_pointer), dimension(maxparms) :: x_new
  type(spec_parms), target, save :: parms, parmslow, parmshigh, parmsnew !The parameters of the EPR simulation to anneal
  type(spec_parms), dimension(5), save :: isotopomer ! Parameters for isotopic derivatives of the main radical
  real, dimension(5,10000) :: isoref ! Spectra of the isotopic derivatives
  integer :: nD = -1, nE = -1, niso(maxiso) = -1, iso, num_isotopomers, isopts(5)
  real :: init_temperature
end module global
