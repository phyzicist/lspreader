[Title]
simulation_title "<TITLE>"
;
[Control]
;Courant limit
courant_multiplier 0.9
;Time-advance
 time_limit <TIMELIM_NS> ; # How many ns to run?
;time_limit 0.2e-6 
;number_of_steps 2

 time_step_ns <TIMESTEP_NS> ; 0.2e-4 is 1/16th 1um laser light; 0.1e-4 is 1/30th
;Restarts
 restart_interval_ns 75000e-6 ; probably much longer than max dump time
 maximum_restart_dump_time 11.5 ;in hours
;Parallel Processing
 balance_interval_ns 0
 load_balance_flag OFF
;Field Solution and Modification
 time_bias_coefficient 0
 time_bias_iterations 1
;Implicit Field Algorithm
 error_current_filtering_parameter 0.95
 implicit_iterations 10
 implicit_tolerance 1.e-5
;Matrix Solution Algorithm
 preconditioner JACOBI
 linear_solution GMRES
;Fluid Physics Algorithm
 fluid_electron_streaming_factor 0.1
 fluid_ion_streaming_factor 0.01 ;Tony insists this is 0.01 instead of 0.005
 flux_limit_fraction 0.2
;(Diagnostic Output) Flags
 dump_current_density_flag OFF
 dump_fields_flag <FLDDUMP>
 dump_scalars_flag <SCLDUMP>
 dump_velocities_flag OFF
 dump_particles_flag OFF
 dump_number_densities_flag <SCLDUMP> ; Densities in scalars?
 ;dump_time_zero_flag ON ; dump the results of the 'zeroth' time step...does it actually start?
 extract_photons_flag OFF
 dump_particles_flag OFF
;(Diagnostic Output) Dump Intervals
 dump_interval <SKIP_T>
 dump_steps
1 
end ; Always include the first time step as a dump.
 spatial_skip_x <SKIP_X>
 spatial_skip_y 1
 spatial_skip_z <SKIP_Z>
 probe_interval 1
;(Diagnostic Output) Formats
 photon_output_format ASCII
 target_output_format ASCII
 use_its_format_flag OFF
 print_region_flag OFF
;(Diagnostic Output) Movie Controls
 particle_movie_interval_ns <PMOV_NS> ; (no dumps if equal to 1.e+9) was 0.2e-6
 particle_movie_components Q X Y Z VX VY VZ XI YI ZI
;Numerical Checks and Reports
 domain_boundary_check ON
 report_timing_flag ON
 dump_timing_flag ON
;
[Grid]
;
grid1
xmin             <XMIN_CM>
xmax              <XMAX_CM>
x-cells           <XCELLS>
;
zmin             <ZMIN_CM>
zmax              <ZMAX_CM>
z-cells           <ZCELLS>
;
[Regions]
;
region1
xmin <XMIN_CM>
xmax  <XMAX_CM>

zmin <ZMIN_CM>
zmax  <ZMAX_CM>

number_of_domains <NDOMS>
split_direction ZSPLIT 
number_of_cells AUTO
;

[Boundaries]
;back this is the laser
outlet
from <XMIN_CM> -0.0000 <ZMIN_CM>
to   <XMIN_CM>  0.0000  <ZMAX_CM>
phase_velocity 1.0
drive_model LASER
reference_point <FOCX> 0 0 ; focal point position
;direction 0 0 0
;magnitude 1.0
;wavelength 0.8e-4 ; 800 nm
;spotsize 2.26e-4 ;these replace the laser analytic function
components 0 0 1
phases 0 0 0 ; polarization 1.1781
temporal_function 1
analytic_function 2
time_delay 0.0
;front
outlet
from  <XMAX_CM> -0.0000 <ZMIN_CM>
to    <XMAX_CM>  0.0000  <ZMAX_CM>
phase_velocity 1.0
drive_model NONE
;right
outlet
from <XMIN_CM> -0.0000 <ZMAX_CM>
to    <XMAX_CM>  0.0000 <ZMAX_CM>
phase_velocity 1.0
drive_model NONE
;left
outlet
from <XMIN_CM> -0.0020 <ZMIN_CM> ; Why is Y 0.002?
to    <XMAX_CM>  0.0020 <ZMIN_CM>
phase_velocity 1.0
drive_model NONE

;;;;;;;;;;;;;;;;
;; species
;;;;;;;;;;;;;;;;
[Materials]
 ;
 material oxygen
 atomic_number 8
 atomic_weight 16
 ;ionization_potential 35.1; eV
 ionization_potential 13.6; eV
 specific_heat 4.186; J/gK
 thermal_conductivity 0.0058; (W / mK)/100 

[Particle Species]
species1 ; neutral O
charge 0
mass 29392
atomic_number 8
migrant_species_flag off
implicit_species_flag on
particle_motion_flag on
particle_forces_option PRIMARY
transverse_weighting_flag on
particle_kinematics_option STANDARD
scattering_flag off
implicit_filtering_parameter 0.1
selection_ratio 0.05
;
species2 ; O+
charge +1
mass 29391
atomic_number 8
migrant_species_flag off
implicit_species_flag on
particle_motion_flag on
particle_forces_option PRIMARY
transverse_weighting_flag on
particle_kinematics_option STANDARD
scattering_flag off
implicit_filtering_parameter 0.1
selection_ratio 0.01
;
species3 ; O++
charge +2
mass 29390
atomic_number 8
migrant_species_flag off
implicit_species_flag on
particle_motion_flag on
particle_forces_option PRIMARY
transverse_weighting_flag on
particle_kinematics_option STANDARD
scattering_flag off
implicit_filtering_parameter 0.1
selection_ratio 0.01
;
species4 ; O 3+
charge +3
mass 29389
atomic_number 8
migrant_species_flag off
implicit_species_flag on
particle_motion_flag on
particle_forces_option PRIMARY
transverse_weighting_flag on
particle_kinematics_option STANDARD
scattering_flag off
implicit_filtering_parameter 0.1
selection_ratio 0.01
;
species5 ; O 4+
charge +4
mass 29388
atomic_number 8
migrant_species_flag off
implicit_species_flag on
particle_motion_flag on
particle_forces_option PRIMARY
transverse_weighting_flag on
particle_kinematics_option STANDARD
scattering_flag off
implicit_filtering_parameter 0.1
selection_ratio 0.01
;
species6 ; O 5+
charge +5
mass 29387
atomic_number 8
migrant_species_flag off
implicit_species_flag on
particle_motion_flag on
particle_forces_option PRIMARY
transverse_weighting_flag on
particle_kinematics_option STANDARD
scattering_flag off
implicit_filtering_parameter 0.1
selection_ratio 0.01
;
species7 ; O 6+
charge +6
mass 29386
atomic_number 8
migrant_species_flag off
implicit_species_flag on
particle_motion_flag on
particle_forces_option PRIMARY
transverse_weighting_flag on
particle_kinematics_option STANDARD
scattering_flag off
implicit_filtering_parameter 0.1
selection_ratio 0.01
;
species8 ; O 7+
charge +7
mass 29385
atomic_number 8
migrant_species_flag off
implicit_species_flag on
particle_motion_flag on
particle_forces_option PRIMARY
transverse_weighting_flag on
particle_kinematics_option STANDARD
scattering_flag off
implicit_filtering_parameter 0.1
selection_ratio 0.01
;
species9 ; O 8+
charge +8
mass 29384
atomic_number 8
migrant_species_flag off
implicit_species_flag on
particle_motion_flag on
particle_forces_option PRIMARY
transverse_weighting_flag on
particle_kinematics_option STANDARD
scattering_flag off
implicit_filtering_parameter 0.1
selection_ratio 0.01
;
species10 ; kinetic electrons
charge -1
mass 1.0
migrant_species_flag off
implicit_species_flag on
particle_motion_flag on
particle_forces_option PRIMARY
transverse_weighting_flag on
particle_kinematics_option STANDARD
scattering_flag on
implicit_filtering_parameter 0.1
selection_ratio 0.01
;
species11 ; protons
charge +1
mass 1836
atomic_number 1
migrant_species_flag off
implicit_species_flag on
particle_motion_flag on
particle_forces_option PRIMARY
transverse_weighting_flag on
particle_kinematics_option STANDARD
scattering_flag off
implicit_filtering_parameter 0.1
selection_ratio 0.01

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

[Particle Creation]

;; initial states ;;

plasma ; O+
from <PC_XMIN_CM> -0.0000 <PC_ZMIN_CM>
to    <PC_XMAX_CM>  0.0000  <PC_ZMAX_CM>
species 2
movie_tag 3
unbound off
discrete_numbers 3 3 3
density_function 6
reference_point 0 0 0
density_flags 1 0 1
momentum_flags 0 0 0
thermal_energy 1
movie_fraction 0.000
;
plasma ; e-
from <PC_XMIN_CM> -0.0000 <PC_ZMIN_CM>
to    <PC_XMAX_CM>  0.0000  <PC_ZMAX_CM>
species 10
movie_tag 3
unbound off
discrete_numbers 3 3 3
density_function 5
reference_point 0 0 0
density_flags 1 0 1
momentum_flags 0 0 0
thermal_energy 1
movie_fraction 0.050
;
plasma ; p
from <PC_XMIN_CM> -0.0000 <PC_ZMIN_CM>
to    <PC_XMAX_CM>  0.0000  <PC_ZMAX_CM>
species 11
movie_tag 3
unbound off
discrete_numbers 3 3 3
density_function 7
reference_point 0 0 0
density_flags 1 0 1
momentum_flags 0 0 0
thermal_energy 1
movie_fraction 0.000

;; ionization states ;; ;

higherstate              ; O -> O+
from <XMIN_CM> -0.0000 <ZMIN_CM>
to    <XMAX_CM>  0.0000  <ZMAX_CM>
interval 1
species 1
ion_species 2
movie_tag 5
electron_species 10
movie_tag 3
ionization_potential 13.6
cross_sections
  0.0
  0.0
  0.0
  0.0
  0.0
  0.0
  0.0
  0.0
  0.0
  0.0
  0.0
end
movie_fraction 0.0
;
higherstate              ; O+ -> O++
from <XMIN_CM> -0.0000 <ZMIN_CM>
to    <XMAX_CM>  0.0000  <ZMAX_CM>
interval 1
species 2
ion_species 3
movie_tag 5
electron_species 10
movie_tag 3
ionization_potential 35.1
cross_sections
  0.0
  0.0
  0.0
  0.0
  0.0
  0.0
  0.0
  0.0
  0.0
  0.0
  0.0
end
movie_fraction 0.0
;
higherstate              ; O++ -> O 3+
from <XMIN_CM> -0.0000 <ZMIN_CM>
to    <XMAX_CM>  0.0000  <ZMAX_CM>
interval 1
species 3
ion_species 4
movie_tag 5
electron_species 10
movie_tag 3
ionization_potential 54.9
cross_sections
  0.0
  0.0
  0.0
  0.0
  0.0
  0.0
  0.0
  0.0
  0.0
  0.0
  0.0
end
movie_fraction 0.0
;
higherstate              ; O 3+ -> O 4+
from <XMIN_CM> -0.0000 <ZMIN_CM>
to    <XMAX_CM>  0.0000  <ZMAX_CM>
interval 1
species 4
ion_species 5
movie_tag 5
electron_species 10
movie_tag 3
ionization_potential 77.4
cross_sections
  0.0
  0.0
  0.0
  0.0
  0.0
  0.0
  0.0
  0.0
  0.0
  0.0
  0.0
end
movie_fraction 0.0
;
higherstate              ; O 4+ -> O 5+
from <XMIN_CM> -0.0000 <ZMIN_CM>
to    <XMAX_CM>  0.0000  <ZMAX_CM>
interval 1
species 5
ion_species 6
movie_tag 5
electron_species 10
movie_tag 3
ionization_potential 113.9
cross_sections
  0.0
  0.0
  0.0
  0.0
  0.0
  0.0
  0.0
  0.0
  0.0
  0.0
  0.0
end
movie_fraction 0.0
;
;
higherstate              ; O 5+ -> O 6+
from <XMIN_CM> -0.0000 <ZMIN_CM>
to    <XMAX_CM>  0.0000  <ZMAX_CM>
interval 1
species 6
ion_species 7
movie_tag 5
electron_species 10
movie_tag 3
ionization_potential 138.1
cross_sections
  0.0
  0.0
  0.0
  0.0
  0.0
  0.0
  0.0
  0.0
  0.0
  0.0
  0.0
end
movie_fraction 0.0
;
higherstate              ; O 6+ -> O 7+
from <XMIN_CM> -0.0000 <ZMIN_CM>
to    <XMAX_CM>  0.0000  <ZMAX_CM>
interval 1
species 7
ion_species 8
movie_tag 5
electron_species 10
movie_tag 3
ionization_potential 739.3
cross_sections
  0.0
  0.0
  0.0
  0.0
  0.0
  0.0
  0.0
  0.0
  0.0
  0.0
  0.0
end
movie_fraction 0.0
;
higherstate              ; O 7+ -> O 8+
from <XMIN_CM> -0.0000 <ZMIN_CM>
to    <XMAX_CM>  0.0000  <ZMAX_CM>
interval 1
species 8
ion_species 9
movie_tag 5
electron_species 10
movie_tag 3
ionization_potential 871.4
cross_sections
  0.0
  0.0
  0.0
  0.0
  0.0
  0.0
  0.0
  0.0
  0.0
  0.0
  0.0
end
movie_fraction 0.0

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

[Particle Extraction]
;
extract1
species 10
direction X
maximum_number 1000000000
start_time 0.0
stop_time 1
at <XMIN_CM> 0 0
;
extract2
species 10
direction X
maximum_number 1000000000
start_time 0.0
stop_time 1
at  <XMAX_CM> 0 0
;
extract3
species 10
direction Z
maximum_number 1000000000
start_time 0.0
stop_time 1
at 0 0 <ZMAX_CM>
;
extract4
species 10
direction Z
maximum_number 1000000000
start_time 0.0
stop_time 1
at 0 0 <ZMIN_CM>
;

;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Functions
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

[Functions]
function1 ; laser temporal function
type 30
data_file <LAS_DAT>
independent_variable_multiplier <LAS_INDEP> ; =2xFWHM,  30 fs FWHM pulse
;dependent_variable_multiplier 1.736e8  ; = Emax in kV/cm units, 1.736e8 => 4*10^19 W/cm^2
;dependent_variable_multiplier 8.68e5  ; = Emax in kV/cm units, 8.68e5 => 10^15 W/cm^2
;dependent_variable_multiplier 2.75e7  ; = Emax in kV/cm units, 2.75e7 => 10^18 W/cm^2
dependent_variable_multiplier <LAS_DEP> ; = Emax in kV/cm units, 4.763e7 => 3 x 10^18 W/cm^2
;dependent_variable_multiplier 2.02e7  ; = Emax in kV/cm units, 2.02e7 => 5.4 x 10^17 W/cm^2

function2 ;laser analytic function for lsp v10
type 19   ; \lambda spotsize
coefficients <WLEN_CM> <SPOT_CM> end
;coefficients 0.8e-4 2.174e-4 end

;;
function5 ; electrons
type 40
data_file <DENS_DAT>
independent_variable_multiplier 1.0
dependent_variable_multiplier 1.0
;;
function6 ; Oxygen
type 40
data_file <DENS_DAT>
independent_variable_multiplier 1.0
dependent_variable_multiplier 0.33
;;
function7 ; Protons
type 40
data_file <DENS_DAT>
independent_variable_multiplier 1.0
dependent_variable_multiplier 0.67
;
;
[Probes]
probe1 ; 
point E X
at <XMIN_CM> 0 0
;
probe2 ; 
point E Y
at <XMIN_CM> 0 0
;                                        ;
probe3 ; 
point E Z
at <XMIN_CM> 0 0
;
probe4 ;
point B X
at <XMIN_CM> 0 0
;                                        
probe5 ;
point B Y
at <XMIN_CM> 0 0
                                        ;
probe6 ;
point B Z
at <XMIN_CM> 0 0

probe7
performance cpu_time
;
probe8
energy net_energy
;
