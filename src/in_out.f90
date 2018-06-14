!> @file      coalescence/src/in_out.f90
module in_out
    use constants
    use globals
    use ioutils, only:newunit,str
    implicit none
    private
    integer :: ii
    public set_paths,read_inputs,save_integration_header,&
        save_integration_step,save_integration_close,save_growth_of_bubble
contains
!********************************BEGINNING*************************************
!> Sets paths to all files.
!!
!! Holds the names and paths of all files.
subroutine set_paths
    fileplacein='./inputs/'
    fileplaceout='./'
    inputs='input.json'
    outputs='outputs.csv'   
    inputs=TRIM(ADJUSTL(fileplacein))//TRIM(ADJUSTL(inputs))
    outputs=TRIM(ADJUSTL(fileplaceout))//TRIM(ADJUSTL(outputs))
end subroutine set_paths
!***********************************END****************************************
!> Sets input values from the input file.
!!
!! Quantites are saved in global variables.
subroutine read_inputs
    use fson
    use fson_value_m, only: fson_value_get
    character(len=1024) :: strval
    type(fson_value), pointer :: json_data
    write(*,*) 'loading input file ',TRIM(inputs)
    json_data => fson_parse(inputs)
    call fson_get(json_data, "coalescence.integrator", strval)
    if (strval=="dlsode") then
        integrator=1
    elseif (strval=="dlsodes") then
        integrator=2
    else
        stop 'unknown integrator'
    endif
    call fson_get(json_data, "coalescence.method", strval)
    if (strval=="nonstiff") then
        if (integrator==1 .or. integrator==2) then
            int_meth=10
        endif
    elseif (strval=="stiff") then
        if (integrator==1 .or. integrator==2) then
            if (integrator==1) then
                int_meth=22
            elseif (integrator==2) then
                int_meth=222
            endif
        endif
    else
        stop 'method can be either stiff or nonstiff'
    endif
    call fson_get(json_data, "coalescence.number_of_bubbles", number_of_bubbles)
    call fson_get(json_data, "coalescence.initialTime", tstart)
    call fson_get(json_data, "coalescence.finalTime", tend)

    call fson_get(json_data, "coalescence.outerTimeSteps", its)
    call fson_get(json_data, "coalescence.maxInnerTimeSteps", maxts)
    call fson_get(json_data, "coalescence.relativeTolerance", rel_tol)
    call fson_get(json_data, "coalescence.absoluteTolerance", abs_tol)

    call fson_get(json_data, "physicalProperties.pressure", pamb)
    call fson_get(&
        json_data, "physicalProperties.blowingAgents.CO2.molarMass", Mbl)
    call fson_get(json_data, "physicalProperties.polymer.heatCapacity", cppol)

    call fson_get(&
        json_data, &
        "physicalProperties.blowingAgents.CO2.heatCapacityInLiquidPhase", &
        cpbll)
    call fson_get(&
        json_data, &
        "physicalProperties.blowingAgents.CO2.heatCapacityInGaseousPhase", &
        cpblg)
    call fson_get(&
        json_data, &
        "physicalProperties.blowingAgents.CO2.evaporationHeat", dHv)
    call fson_get(json_data, "initialConditions.temperature", T0)
   
 
    call fson_get(json_data, "initialConditions.concentrations.water", W0)
    call fson_get(json_data, "kinetics.gelPoint", gelpointconv)
    call fson_get(&
            json_data, "initialConditions.concentrations.polyol", OH0)
    call fson_get(&
            json_data, "initialConditions.concentrations.isocyanate", NCO0)
    call fson_get(&
            json_data, "kinetics.gellingReaction.frequentialFactor", AOH)
    call fson_get(&
            json_data, "kinetics.gellingReaction.activationEnergy", EOH)
    call fson_get(&
            json_data, "kinetics.blowingReaction.frequentialFactor", AW)
    call fson_get(&
            json_data, "kinetics.blowingReaction.activationEnergy", EW)
    call fson_get(&
        json_data, &
        "initialConditions.concentrations.blowingAgents.CO2", c0)
    call fson_get(&
            json_data, "kinetics.gellingReaction.reactionEnthalpy", dHOH)
    call fson_get(&
            json_data, "kinetics.blowingReaction.reactionEnthalpy", dHW)

    call fson_get(json_data, "physicalProperties.polymer.density", rhop)
    call fson_get(json_data, "physicalProperties.surfaceTension", sigma)
    call fson_get(&
            json_data, &
            "physicalProperties.blowingAgents.CO2.diffusivity", D)
    call fson_get(&
            json_data, &
            "physicalProperties.blowingAgents.CO2.solubility", KH)
    call fson_destroy(json_data)
    write(*,*) 'done: inputs loaded'
    write(*,*)
end subroutine read_inputs
!***********************************END****************************************


!********************************BEGINNING*************************************
!> Opens output files and writes a header.
!!
!! Names of variables are written on the first line in the file.
!! Soubroutine is called before the integration commences.
subroutine save_integration_header
	ii=0
    open (unit=newunit(fi1), file = outputs)
    write(fi1,*) 't', ',', 'radius1', ',', 'radius2', ',', 'radius3', ',', 'viscosity', ',', 'number_of_coalesced_bubbles'
    !open (unit=newunit(fi2), file = outputs_GR)
    !write(fi2,'(1000A23)') '#GrowthRate1', 'GrowthRate2', 'temperature', &
    !    'bubbleRadius','c1','c2','p1','p2'
    !open (unit=newunit(fi4), file = outputs_c)
    !open (unit=newunit(fi5), file = outputs_drain)
    !write(fi5,'(1000A24)') '#time','radius','porosity','viscosity'
end subroutine save_integration_header
!***********************************END****************************************


!********************************BEGINNING*************************************
!> Writes an integration step to output file.
!!
!! Subroutine is called once per time step.
subroutine save_integration_step(iout)
	integer, intent(in) :: iout !< number of the time step
	integer :: i,ii
	character*13 :: filename

	ii=iout/35

	if(ii .LT. 10) then
		write(filename,'("TEST",I1,".TXT")')ii
	else
		write(filename,'("TEST",I2,".TXT")')ii
	endif
	open (unit=newunit(ii), file = filename)
    	write(ii,*) 'x coord',', ', 'y coord', ', ', 'z coord', ', ', 'radius'
	do i=1,number_of_bubbles
		if (exist_bubble(i)) then
			write(ii,*) x_coord(i),', ', y_coord(i), ', ', z_coord(i), ', ', y(i)
		endif
	enddo
	close(unit=(ii))

   ! write(fi2,"(1000es23.15)") grrate, temp, radius, avconc, pressure
   ! write(fi4,"(1000es23.15)") (Y(fceq+i+1),i=0,ngas*p,ngas)
   ! write(fi5,"(1000es24.15e3)") time,radius,porosity,eta
    ! save arrays, which are preserved for future use
   ! etat(iout,1)=time
   ! etat(iout,2)=eta
   ! port(iout,1)=time
   ! port(iout,2)=porosity
   ! init_bub_rad(iout,1)=time
   ! init_bub_rad(iout,2)=radius
end subroutine save_integration_step

subroutine save_growth_of_bubble(t)
	real(dp), intent(in) :: t

	write(fi1,*) t, ',', y(25), ',', y(20), ',', y(125), ',', eta, ',' , number_of_coalesced_bubbles

end subroutine save_growth_of_bubble



subroutine save_integration_close(iout)
    integer, intent(in) :: iout !< number of the time step
    integer :: i
    real(dp), dimension(:,:), allocatable :: matr
    open (unit=newunit(fi2), file = 'diameter.csv')
    	write(fi2,*)  'diameter'
	do i=1,number_of_bubbles
		if (exist_bubble(i)) then
			write(iout,*) 2*y(i)
		endif
	enddo
	close(unit=(fi2))
    close(fi1)
  !  close(fi2)
  !  close(fi4)
  ! close(fi5)
end subroutine save_integration_close
end module in_out
