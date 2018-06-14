!> @file	 coalescence/src/integration.f90
module integration
	use constants
	use globals
	use list
	use in_out, only: save_integration_header,save_integration_step,save_integration_close,save_growth_of_bubble
	use model, only:odesystem, calc_conc_and_eta,coalesce
	implicit none
	private
        integer,parameter :: seed = 86456 !111111
	integer :: meth, itmeth, iatol, jpretype, igstype, maxl, ier!, itask
        !integer(c_long) :: ipar(1), mu, ml, iouts(25), maxtss
	real(dp) :: rout(10), rpar(1), delt!, tout, atol, rtol, t
	integer :: iout, iopt, istate, itask, itol, liw, lrw, nnz, lenrat, neq, mf
        real(dp) :: jac,tout,rtol,atol,t
        real(dp), dimension(:), allocatable :: rwork!y
        integer, dimension(:), allocatable :: iwork
	abstract interface
		subroutine sub(neq,t,y,ydot)
 	    		use constants
            		integer, intent(in) :: neq
            		real(dp), intent(in) :: t
	    		real(dp), intent(in) :: y(neq)
            		real(dp), intent(out) :: ydot(neq)
		end subroutine sub
	end interface
	procedure (sub), pointer :: sub_ptr => odesystem
	public integrate
contains

subroutine set_integrator
	!write(*,*) integrator
	!stop
        if (integrator==1 .or. integrator==2) then
        mf=int_meth
        select case(integrator)
        case(1)
            select case(mf)
            case(10)
                allocate(rwork(20+16*neq),iwork(20))
            case(22)
                allocate(rwork(22+9*neq+neq**2),iwork(20+neq))
            case default
                stop 'unknown mf'
            end select
        case(2)
            select case(mf)
            case(10)
                allocate(rwork(20+16*neq),iwork(30))
            case(222)
                nnz=neq**10 !not sure, smaller numbers make problems for low p
                lenrat=2 !depends on dp
                allocate(rwork(int(20+(2+1._dp/lenrat)*nnz+(11+9._dp/lenrat)*neq)),&
                    iwork(30))
            case default
                stop 'unknown mf'
            end select
        case default
            stop 'unknown integrator'
        end select
        itask = 4
        istate = 1
        iopt = 1
        rwork(1)=tend+epsilon(tend)*1.0e3_dp
        rwork(5:10)=0
        iwork(5:10)=0
        lrw = size(rwork)
        liw = size(iwork)
        iwork(6)=maxts
        itol = 1 !don't change, or you must declare atol as atol(neq)
        rtol=rel_tol
        atol=abs_tol
	endif
end subroutine set_integrator

subroutine set_initial_conditions
	integer :: i,j
	gelpoint=.false.
        timestep=(tend-tstart)/its
	t = tstart
    	tout = t+timestep
        neq=2*number_of_bubbles+3
    	allocate(y(neq))
        allocate(x_coord(number_of_bubbles))
	allocate(y_coord(number_of_bubbles))
	allocate(z_coord(number_of_bubbles))
	allocate(R0(number_of_bubbles))
	allocate(P0(number_of_bubbles))
	allocate(c_r(number_of_bubbles))
	allocate(pair0(number_of_bubbles))
	allocate(pair(number_of_bubbles))
	allocate(exist_bubble(number_of_bubbles))
	last_pair=1
    	y=0
	teq=2*number_of_bubbles+1
	xOHeq=teq+1
	xWeq=xOHeq+1
	y(teq) = T0   !temperature
	y(xOHeq) = 0.0_dp        !xOH
    	y(xWeq) = 0.0_dp        !xW
	do i=1,number_of_bubbles    !pozor este treba skontrolovat aby sa neprekryvali!!!!!!
		x_coord(i)=rand()*0.93e-4_dp
                y_coord(i)=rand()*0.93e-4_dp
		z_coord(i)=rand()*0.93e-4_dp
                R0(i)= rand()*(8.1e-8_dp-8.0e-8_dp)+8.0e-8_dp!2*sigma/(c0/KH-pamb)!rand()*1e-5_dp !rand()*(5.0e-7_dp-1.0e-7)+1.0e-7_dp
		!write(*,*) R0(i)
		P0(i)= (2*sigma/R0(i)+pamb) !c0/KH!(2*sigma/R0(i)+pamb)  !c0/KH+pamb !pozor to kontrolovat ako to ma byt s pavlovym input!!!!
		!c_r(i)=0.0_dp
		!write(*,*) R0(i), P0(i)
		exist_bubble(i)=.TRUE.
	enddo
	!stop
	!write(*,*) c0/KH+pamb, (2*sigma/R0(1)+pamb)
	!stop
        do i=1,number_of_bubbles
        	y(i) = R0(i) !radius
        enddo
	do i=number_of_bubbles+1,2*number_of_bubbles
		pair0(i-number_of_bubbles)=P0(i-number_of_bubbles)
        	y(i) = 0.0_dp
		!write(*,*) y(i)
	enddo
	number_of_coalesced_bubbles=0
        
end subroutine set_initial_conditions

subroutine integrate
	integer :: i,j,min_idx,max_idx, number_of_existing_bubbles
	integer, dimension(:), allocatable :: hist
	real(dp) :: min_rad,max_rad, deltat_max,min_R
	real(dp) :: ydot(neq) !< derivatives of integrated variables
	call set_initial_conditions
	call set_integrator
	write(*,*) 'integrating...'
	!write(*,*) t, tout
	!stop
	!write(*,*) t,tout
	deltat_max=0.0111_dp !0.11_dp
	call save_integration_header
	do iout = 1,its
		!write(*,*) 'iteration', iout
		call calc_conc_and_eta(iout)
		!write(*,*) c_average
		!stop
		!do i=1,number_of_bubbles
		!	if (exist_bubble(i)) then
		!		write(*,*) i, t, y(i), y(i+number_of_bubbles)-pamb, y(teq), y(xOHeq), y(xWeq),eta, exist_bubble(i)!, (y(i+number_of_bubbles)*y(i)**3/y(teq)-P0(i)*R0(i)**3/T0)
		!	endif
!
		!enddo
		
		
		select case (integrator)
        		case(1)
				!call odesystem(neq,t,y,ydot)
				!write(*,*) t, tout
				call dlsode (sub_ptr, neq, y, t, tout, itol, rtol, atol, itask, &
                                                     istate, iopt, rwork, lrw, iwork, liw, jac, mf)
        		case(2)
				istate=1
				!write(*,*) t, tout
				
        			call dlsodes(sub_ptr, neq, y, t, tout, itol, rtol, atol, itask, &
                				istate, iopt, rwork, lrw, iwork, liw, jac, mf)
				istate=1
				
				!write(*,*) istate
        		case default
        		    stop 'unknown integrator'
        	end select
		!write(*,*) P0(1), y(1+number_of_bubbles), pair(1), pair0(1), pamb, y(1+number_of_bubbles)-pair(1)-pamb

       ! if (gelpoint) then
       !     write(*,'(2x,A,f7.3,A)') 'gel point reached at time t = ',time,' s'
       !     write(*,'(2x,A,f6.1,A)') 'temperature at gel point T = ',temp,' K'
       !     write(*,'(2x,A,f5.3)') 'conversion at gel point X = ',conv
       !     exit
       ! endif
	!write(*,*) 'iterace ', iout
	min_rad=10e10
	max_rad=0.0_dp
	number_of_existing_bubbles=0
		
		if((mod(iout,10000)==0) .or. (iout==1)) then
			do j=1,number_of_bubbles
				if ((exist_bubble(j)) .and. (y(j)<min_rad)) then
					min_rad=y(j)
					min_idx=j
				endif
				if ((exist_bubble(j)) .and. (y(j)>max_rad)) then
					max_rad=y(j)
					max_idx=j
				endif
				if (exist_bubble(j)) then
					number_of_existing_bubbles=number_of_existing_bubbles+1
				endif
			enddo
			write(*,*) 'iteration', iout
			write(*,*) 'time:', t, 'n_coal_bub:', number_of_coalesced_bubbles, 'number_of_existing_bubbles:', number_of_existing_bubbles
			write(*,*) 'min_idx:', min_idx, ' min R:', min_rad, 'max_idx', max_idx, ' max R:', max_rad, ' porosity:', porosity, ' last_pair', last_pair
			write(*,*) 'pressure:', pair(1)+y(1+number_of_bubbles), 'pressure_initial:', P0(1), 'c_r:', c_r(1), 'c_average:', c_average
		endif
		!if((max_rad-min_rad)>1.0e-6_dp) then
	!		help=.true.
!			do j=1,number_of_bubbles
!			write(*,*) y(j), exist_bubble(j),j
!			enddo
!			pause(100)
!		endif
		!help=.false.
	call coalesce(t,iout)
	
	!write(*,*) y


	!if(mod(iout,35)==0) then
	!if(iout==1000) then 
	!	call save_integration_step(iout)
		!stop
	!endif	
	call save_growth_of_bubble(t)
	!write(*,*) t,tout
	if (timestep<deltat_max) then
		timestep=timestep*1.00001_dp
	else 
		timestep=deltat_max
	endif
	tout = t+timestep
	!write(*,*) t,tout
	!stop
	if (gelpoint) then
		do j=1,number_of_bubbles
				if ((exist_bubble(j)) .and. (y(j)<min_rad)) then
					min_rad=y(j)
					min_idx=j
				endif
				if ((exist_bubble(j)) .and. (y(j)>max_rad)) then
					max_rad=y(j)
					max_idx=j
				endif
				if (exist_bubble(j)) then
					number_of_existing_bubbles=number_of_existing_bubbles+1
				endif
		enddo
		write(*,*) 'gelpoint reached'
		write(*,*) 'time:', t, 'number of coalesced bubbles', number_of_coalesced_bubbles, 'max_rad', max_rad, 'min_rad', min_rad, 'porosity', porosity, 'viscosity', eta
		exit
	endif
	!if (porosity>0.8) then
	!	write(*,*) 'close-packing of equal spheres'
	!	write(*,*) 'time:', t, 'number of coalesced bubbles', number_of_coalesced_bubbles, 'max_rad', max_rad, 'min_rad', min_rad, 'porosity', porosity, 'viscosity', eta
	!	exit
	!endif
	!stop
	end do
	
	call save_integration_close(iout)
	allocate(hist(11))
	hist=0
	min_R=50e-6_dp
	do i=1,number_of_bubbles
		if(exist_bubble(i) .and. (2*y(i)<=min_R)) then
			hist(1)=hist(1)+1
		else if((exist_bubble(i) .and. (2*y(i)<=2*min_R)) .and. (2*y(i)>min_R)) then
			hist(2)=hist(2)+1
		else if((exist_bubble(i) .and. (2*y(i)<=3*min_R)) .and. (2*y(i)>2*min_R)) then
			hist(3)=hist(3)+1
		else if((exist_bubble(i) .and. (2*y(i)<=4*min_R)) .and. (2*y(i)>3*min_R)) then
			hist(4)=hist(4)+1
		else if((exist_bubble(i) .and. (2*y(i)<=5*min_R)) .and. (2*y(i)>4*min_R)) then
			hist(5)=hist(5)+1
		else if((exist_bubble(i) .and. (2*y(i)<=6*min_R)) .and. (2*y(i)>5*min_R)) then
			hist(6)=hist(6)+1
		else if((exist_bubble(i) .and. (2*y(i)<=7*min_R)) .and. (2*y(i)>6*min_R)) then
			hist(7)=hist(7)+1
		else if((exist_bubble(i) .and. (2*y(i)<=8*min_R)) .and. (2*y(i)>7*min_R)) then
			hist(8)=hist(8)+1
		else if((exist_bubble(i) .and. (2*y(i)<=9*min_R)) .and. (2*y(i)>8*min_R)) then
			hist(9)=hist(9)+1
		else if((exist_bubble(i) .and. (2*y(i)<=10*min_R)) .and. (2*y(i)>9*min_R)) then
			hist(10)=hist(10)+1
		else
			hist(11)=hist(11)+1
		endif			
	enddo
	write (*,*) 'histogram:', hist/number_of_existing_bubbles
	write(*,*) 'done: integration'
end subroutine integrate



end module integration
