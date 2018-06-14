!> @file      coalescence/src/globals.f90
module globals
    use constants
    implicit none
    character(len=99) :: &
        fileplacein,& !< location of input files
        fileplaceout,& !< location of output files
        inputs,& !< input file
        outputs !< output file with scalar variables
    logical :: &
	gelpoint,& !<gel point reached (t/f)
	help
    integer :: &
        fi1,fi2,fi3,fi4,fi5,& !<file indices
	number_of_bubbles,&
	number_of_coalesced_bubbles,&
        integrator,& !<integrator. 1=dlsode,2=dlsodes
        int_meth,& !< stiff or no? See MF for ODEPACK
        p,& !<number of internal nodes
        maxts,& !<maximum inner time steps between t and t+h
        its,& !<number of outer integration outer time steps
        teq,& !<temperature equation (index)
	last_pair,& !index of last added coalescence pair>
        xOHeq,& !<polyol conversion equation (index)
        xWeq !<water conversion equation (index)
    real(dp) :: &
	c_average,&
	porosity,&
        OH0,& !<initial concentration of polyol (don't set to zero)
        W0,& !<initial concentration of water (can be zero)
        NCO0,& !<initial concentration of isocyanate
        AOH,& !<frequential factor of gelling reaction
        EOH,& !<activation energy of gelling reaction
        AW,& !<frequential factor of blowing reaction
        EW,& !<activation energy of blowing reaction
        dHOH,& !<gelling reaction enthalpy
        dHW,& !<blowing reaction enthalpy
        time,& !<time (s)
        rel_tol,& !<relative tolerance
        abs_tol,& !<absolute tolerance
        rhop,& !<polymer density
        cp,& !<heat capacity of the reaction mixture
        cppol,& !<heat capacity of polymer
        rhobl,& !<density of liquid physical blowing agent
        eta,& !<viscosity
        pamb,& !<ambient pressure (in the liquid)
        sigma,& !<interfacial tension
        timestep,& !<timestep
	tstart,&
	tend,&
        gelpointconv,& !<conversion of polyol at gel point
	KH,& !<Henry constant
	Mbl,& !<blowing agent molar mass
	dHv,& !<evaporation heat of blowing agent
        cpblg,& !<heat capacity of blowing agent in gas phase
        cpbll,& !<heat capacity of blowing agent in liquid phase
	D,&
	t0,&
	c0
    logical, dimension(:), allocatable :: &	
	exist_bubble
    real(dp), dimension(:), allocatable :: &
        y,& !<state variables   
	R0,&
	P0,&
	c_r,&
	pair,&
	pair0,&
	z_coord,&
	y_coord,&
        x_coord
end module globals
