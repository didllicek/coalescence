!>@file      coalescence/src/model.f90
module model
	use constants
	use globals
	implicit none
	type coalescence_pair
        	integer bub1_idx,bub2_idx    
        	real time
		logical did_coalesce
        end type coalescence_pair
	public calc_conc_and_eta,odesystem,coalesce
	type(coalescence_pair) pairs(800**2)

contains

subroutine calc_conc_and_eta(iout)
	integer, intent(in) :: iout
	real(dp) :: Aeta,Eeta,Cg,AA,B !viscosity model constants

	integer :: i
	real(dp) :: dxOH,dxW,v0,vgas
	vgas=0.0_dp

	!c_generated=c_generated+W0*ydot(xWeq)*1e4*timestep	
	!c_average=c0+c_generated	

	do i=1,number_of_bubbles
		if (exist_bubble(i)) then
			vgas=vgas+4*pi*y(i)**3/3
		endif
		v0=(0.93e-4_dp)**3-4*pi*R0(i)**3/3 !attention
		pair(i)=pair0(i)*R0(i)**3/y(i)**3*y(teq)/T0 
	enddo
	porosity=vgas/(v0+vgas)

	c_average=c0+W0*y(xWeq)!*1.2!/10.0_dp!/20!*1e4
	
	do i=1,number_of_bubbles
		c_average=c_average-4*pi*y(i)**3*((y(i+number_of_bubbles))/3/Rg/y(teq))/v0
	enddo
	do i=1,number_of_bubbles
		c_r(i)=KH*(y(i+number_of_bubbles))
		!write(*,*) pair(i), pair0(i), c_r(i), c_average
		!stop
	enddo
	!write(*,*) 'pair:', pair(1), 'pressure of CO2:', y(1+number_of_bubbles), 'pair0:', pair0(1), 'c_r:', c_r(i), 'c_average:', c_average
	
	
	if (y(teq)>500) then
        print*, 'temperature > 500', y(teq)
        return
	endif
	
	if (.not. gelpoint) then
		Aeta=4.1e-8_dp
	        Eeta=38.3e3_dp
        	AA=4.0_dp
            	B=-2.0_dp	
            	eta=Aeta*exp(Eeta/(Rg*y(teq)))*(gelpointconv/(gelpointconv-y(xOHeq)))**(AA+B*y(xOHeq))
		!write(*,*) gelpointconv-y(xOHeq), y(xOHeq), eta
	end if
	!write(*,*) gelpointconv
	if (eta>1.0e10_dp) eta=1.0e10_dp
	if (y(xOHeq)>gelpointconv*0.995_dp) then
            gelpoint=.true.
        endif
	cp=cppol!+c_average*Mbl*cpbll/rhop
	!write(*,*) cp
	!stop
	!dxOH=(AOH*exp(-EOH/Rg/y(teq))*(1-y(xOHeq))*(NCO0-2*W0*y(xWeq)-OH0*y(xOHeq)))
	!write(*,*) dxOH
	!y(xOHeq) = y(xOHeq)+ (AOH*exp(-EOH/Rg/y(teq))*(1-y(xOHeq))*(NCO0-2*W0*y(xWeq)-OH0*y(xOHeq)))*timestep !polyol conversion
	!write(*,*) ydot(xOHeq)
        !if (W0>1e-3) then
            ! water conversion
            ! ydot(xWeq) = AW*exp(-EW/Rg/temp)*(1-y(xWeq))*&
            !     (NCO0-2*W0*y(xWeq)-OH0*y(xOHeq)) 2nd order
	 !dxW=(AW*exp(-EW/Rg/y(teq))*(1-y(xWeq)))
	 !write(*,*) dxW
         !y(xWeq) = y(xWeq)+ (AW*exp(-EW/Rg/y(teq))*(1-y(xWeq)))*timestep !1st order

	 !   write(*,*) ydot(xWeq)
        !endif
        !temperature (enthalpy balance)
	!y(teq) = y(teq)+ (-dHOH*OH0/(rhop*cp)*dxOH-dHW*W0/(rhop*cp)*dxW)*timestep  
	!ydot(teq) = ydot(teq) - dHv(i)*12*pi*Mbl(i)*D(i)*radius**4/&      pripadne vyparna entalpia 
         !       (rhop*cp*Vsh)*(y(fceq+i-1)-KH(i)*y(fpeq+i-1))/(dz(1)/2
	!write(*,*) 'som tu'
        
	

end subroutine calc_conc_and_eta

!! evolution of pressure and radius of bubbles
subroutine  odesystem (neq, t, y, ydot)
	integer, intent(in) :: neq !< number of equations
	real(dp), intent(in) :: t !< time
	real(dp), intent(in) :: y(neq) !< integrated variables
	real(dp), intent(out) :: ydot(neq) !< derivatives of integrated variables
	integer :: i,j
        real(dp) :: min_der,max_der,min_rad,max_rad
	ydot=0


	ydot(xOHeq) = AOH*exp(-EOH/Rg/y(teq))*(1-y(xOHeq))*(NCO0-2*W0*y(xWeq)-OH0*y(xOHeq)) !polyol conversion

        ydot(xWeq) = AW*exp(-EW/Rg/y(teq))*(1-y(xWeq)) !1st order

	ydot(teq) = -dHOH*OH0/(rhop*cp)*ydot(xOHeq)-dHW*W0/(rhop*cp)*ydot(xWeq)
	
        do i=1,number_of_bubbles
		
		if ((exist_bubble(i))) then! .and. ((y(teq)-T0) .gt. 0.0_dp)) then
			ydot(i)=y(i)/(4*eta)*(y(i+number_of_bubbles)+pair(i)-pamb-2*sigma/y(i))
			!write(*,*) c_r, c_average
                	if ((c_average-c_r(i)) .gt. 0.0_dp) then
				ydot(i+number_of_bubbles)=6.0_dp*D*(Rg*y(teq)**2*(c_average-c_r(i))**2*y(i))/((y(i+number_of_bubbles)+pair(i))*y(i)**3/y(teq)-P0(i)*R0(i)**3/T0)&
						-3*(y(i+number_of_bubbles)+pair(i))/y(i)*ydot(i)+(y(i+number_of_bubbles)+pair(i))/y(teq)*ydot(teq)

                	else if ((c_average-c_r(i)) .lt. 0.0_dp) then
				ydot(i+number_of_bubbles)=-6.0_dp*D*(Rg*y(teq)**2*(c_average-c_r(i))**2*y(i))/((y(i+number_of_bubbles)+pair(i))*y(i)**3/y(teq)-P0(i)*R0(i)**3/T0)&
						-3*(y(i+number_of_bubbles)+pair(i))/y(i)*ydot(i)+(y(i+number_of_bubbles)+pair(i))/y(teq)*ydot(teq)	
			else 
				ydot(i+number_of_bubbles)=-3*(y(i+number_of_bubbles)+pair(i))/y(i)*ydot(i)+(y(i+number_of_bubbles)+pair(i))/y(teq)*ydot(teq)	
				!write(*,*) 'here I am'		
				
			endif
			
			
		else
			ydot(i)=0.0_dp
			ydot(i+number_of_bubbles)=0.0_dp
		endif 
		if ((y(i)<1e-8_dp) .or. (ydot(i)<0)) then
			ydot(i+number_of_bubbles)=0.0_dp
			ydot(i)=0.0_dp
			exist_bubble(i)=.false.
		endif
        enddo
	!write(*,*) ydot(1+number_of_bubbles)
	
  
end subroutine odesystem


subroutine coalesce(t,iout)
	real(dp) :: t
	integer :: i,j,idx1,idx2,bub1_idx,bub2_idx,k,iout
        real(dp) :: edge,r,distance,press_new,tc,nair_1,nair_2,nCO2_1,nCO2_2,xair,xCO2
	logical :: already_in_list,dc
	type(coalescence_pair) :: pair_bb


	do i=1,last_pair-1
		idx1=pairs(i)%bub1_idx
		idx2=pairs(i)%bub2_idx
		if (exist_bubble(idx1) .and. exist_bubble(idx2)) then
			tc=pairs(i)%time
			dc=pairs(i)%did_coalesce		
			if (((tc < t) .and. (last_pair>1)) .and. (.not. dc)) then		
				!write(*,*) tc	
				!stop
				pairs(i)%did_coalesce=.TRUE.
				exist_bubble(idx2)=.FALSE.
				x_coord(idx1)=(x_coord(idx1)+x_coord(idx2))/2
				y_coord(idx1)=(y_coord(idx1)+y_coord(idx2))/2
				z_coord(idx1)=(z_coord(idx1)+z_coord(idx2))/2
				write(*,*) 'b1:', idx1, 'povodne:', 'R1:', y(idx1), 'P1:', y(idx1+number_of_bubbles), 'pair1:', pair(idx1)
				write(*,*) 'b2:', idx2, 'povodne:', 'R2:', y(idx2), 'P2:', y(idx2+number_of_bubbles), 'pair2:', pair(idx2)
				nair_1=pair(idx1)*4*pi*y(idx1)**3/3/Rg/y(teq)
				nair_2=pair(idx2)*4*pi*y(idx2)**3/3/Rg/y(teq)
				nCO2_1=y(idx1+number_of_bubbles)*4*pi*y(idx1)**3/3/Rg/y(teq)
				nCO2_2=y(idx2+number_of_bubbles)*4*pi*y(idx2)**3/3/Rg/y(teq)
				xair=(nair_1+nair_2)/(nair_1+nair_2+nCO2_1+nCO2_2)
				xCO2=(nCO2_1+nCO2_2)/(nair_1+nair_2+nCO2_1+nCO2_2)
				press_new=2.0_dp/(1.0_dp/(y(idx1+number_of_bubbles)+pair(idx1))+1.0_dp/(y(idx2+number_of_bubbles)+pair(idx2)))
				!write(*,*) xair, xCO2, press_new
				y(idx1)=(3.0_dp/4/pi/(press_new)*&
				(4*pi/3.0*y(idx1)**3*(y(idx1+number_of_bubbles)+pair(idx1))+&
				4*pi/3.0*y(idx2)**3*(y(idx2+number_of_bubbles)+pair(idx2))))**(1./3)
				pair(idx1)=xair*press_new			
				y(idx1+number_of_bubbles)=xCO2*press_new
				write(*,*) 'new:', idx1, 'R:', y(idx1), 'P:', y(idx1+number_of_bubbles), 'pair:', pair(idx1)
				!write(*,*) i, idx1, idx2
				number_of_coalesced_bubbles=number_of_coalesced_bubbles+1
				pairs(i)%did_coalesce=.true.
				!do j=1,last_pair
				!	if(.not.(pairs(j)%did_coalesce)) then				
				!		if (pairs(j)%bub1_idx==idx2) then
				!			bub1_idx=pairs(j)%bub1_idx
				!			distance=sqrt((x_coord(idx1)-x_coord(bub1_idx))**2+(y_coord(idx1)-y_coord(bub1_idx))**2&
				!				+(z_coord(idx1)-z_coord(bub1_idx))**2)	
				!			if(distance<(y(idx1)+y(bub1_idx))) then
				!				pairs(j)%bub1_idx=idx1
				!				pairs(j)%time=t+6*eta*500000/&
				!				abs(y(idx1+number_of_bubbles)+pair(idx1+number_of_bubbles)-y(bub1_idx+number_of_bubbles)-&
				!					pair(bub1_idx+number_of_bubbles))	
				!			else 
				!				pairs(j)%did_coalesce=.true.
				!			endif
				!		endif
				!		if (pairs(j)%bub2_idx==idx2) then
				!			bub2_idx=pairs(j)%bub2_idx
				!			distance=sqrt((x_coord(idx1)-x_coord(bub2_idx))**2+(y_coord(idx1)-y_coord(bub2_idx))**2+&
				!					(z_coord(idx1)-z_coord(bub2_idx))**2)	
				!			if(distance<(y(idx1)+y(bub2_idx))) then
				!				pairs(j)%bub2_idx=idx1
				!				pairs(j)%time=t+6*eta*500000/&
				!				abs(y(idx1+number_of_bubbles)+pair(idx1+number_of_bubbles)-y(bub2_idx+number_of_bubbles)-&
				!						pair(bub2_idx+number_of_bubbles))	
				!			else 
				!				pairs(j)%did_coalesce=.true.
				!			endif						
				!		endif	
				!	endif			
				!enddo
							
			endif

		else
			pairs(i)%did_coalesce=.TRUE.
		endif
	enddo

	do i=1,number_of_bubbles
		do j=1,number_of_bubbles
			already_in_list=.false.
			do k=1,last_pair
				if (((pairs(k)%bub1_idx==i) .and. (pairs(k)%bub2_idx==j)) .or. ((pairs(k)%bub1_idx==j) .and. (pairs(k)%bub2_idx==i))) then
					already_in_list=.true.
				endif
			end do
			if ((i>j) .and. (.not. (already_in_list))) then
				if((exist_bubble(i)) .and. (exist_bubble(j))) then
					distance=sqrt((x_coord(i)-x_coord(j))**2+(y_coord(i)-y_coord(j))**2+(z_coord(i)-z_coord(j))**2)
					if(distance<y(i)+y(j)) then
						if (y(i)>=y(j)) then
							pair_bb%time=t+6*eta*500000/&
								abs(y(i+number_of_bubbles)+pair(i+number_of_bubbles)-y(j+number_of_bubbles)-pair(j+number_of_bubbles))
							pair_bb%bub1_idx=i
							pair_bb%bub2_idx=j	
							pair_bb%did_coalesce=.false.			
						else
							pair_bb%time=t+6*eta*500000/&
								abs(y(i+number_of_bubbles)+pair(i+number_of_bubbles)-y(j+number_of_bubbles)-pair(j+number_of_bubbles))
							pair_bb%bub1_idx=j
							pair_bb%bub2_idx=i	
							pair_bb%did_coalesce=.false.				
						endif	
						if (pair_bb%time<65.0_dp) then  !pozor miesto 70 by malo byt tend
							pairs(last_pair)=pair_bb
							!write(*,*) pairs(last_pair)%time, pairs(last_pair)%bub1_idx, pairs(last_pair)%bub2_idx
							!write(*,*) 'b1', i, 'b2', j, 'time', t, 'time coal:', pairs(last_pair)%time
							last_pair=last_pair+1	
						endif						
					endif		
				endif	
			 endif
		enddo
	enddo
	!if(mod(iout,200)==0) then
	!	do i=1,last_pair
	!		idx1=pairs(i)%bub1_idx
	!		idx2=pairs(i)%bub2_idx
	!		write(*,*) i, idx1, idx2, pairs(i)%did_coalesce, pairs(i)%time
	!	enddo
	!endif
end subroutine coalesce
end module model
