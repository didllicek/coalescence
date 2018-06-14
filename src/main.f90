program coalescence
	use integration,only:integrate
	use in_out,only:set_paths,read_inputs
	implicit none
	call set_paths	
	call read_inputs	
	call integrate
end program
