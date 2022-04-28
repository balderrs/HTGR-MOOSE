# ==============================================================================
# Blowdown of a 1/20th Scaled-down HTGR UI
# ------------------------------------------------------------------------------
# Idaho Falls, INL, July 23, 2021 8:00 AM
# Author(s): Silvino A. Balderrama Prieto, Guillaume L. Giudicelli
# ------------------------------------------------------------------------------
# Description: The following model simulates the blowdown of a scaled-down HTGR
# as a result of a break in the RPV. The porous media compressible finite volume
# KT is implemented to solve for the transient depressurization.
# ==============================================================================
# MODEL PARAMETERS
# ------------------------------------------------------------------------------
ambient_pressure = 101325 # Ambient pressure [Pa]
ambient_temperature = ${fparse (22.4+273.15)} # Initial temperature [K]
# Containment ------------------------------------------------------------------
T_Cav = ${fparse (65.65+273.15)} # Initial temperature [K]
P_Cav = 101325 # Initial pressure [Pa]
# RPV --------------------------------------------------------------------------
T_RPV = ${fparse 197.7+273.15} # Initial temperature [K]
P_RPV = ${fparse 168.296*6894.75729} # Initial pressure [Pa] (153.6psig)
# ------------------------------------------------------------------------------
# Global Parameters
# ------------------------------------------------------------------------------
[GlobalParams]
  fp = fp
  gravity = '0 -9.81 0'
[]
# ------------------------------------------------------------------------------
# Mesh: Import mesh file
# ------------------------------------------------------------------------------
[Mesh]
  [fmg]
    type = FileMeshGenerator
    file = './test.e'
  []
[]
# ------------------------------------------------------------------------------
# Variables
# ------------------------------------------------------------------------------
[Variables]
  [rho]
    type = MooseVariableFVReal
    block = '1 2'
  []
  [rhou]
    type = MooseVariableFVReal
    block = '1 2'
  []
  [rhov]
    type = MooseVariableFVReal
    block = '1 2'
  []
  [rhow]
    type = MooseVariableFVReal
    block = '1 2'
  []
  [rho_et]
    type = MooseVariableFVReal
    block = '1 2'
  []
[]
# ------------------------------------------------------------------------------
# FVKernels
# ------------------------------------------------------------------------------
[FVKernels]
  # ------ Conservation of Mass --------
  [mass_time]
    type = FVTimeKernel
    variable = rho
  []
  [mass_advection]
    type = PCNSFVKT
    variable = rho
    eqn = "mass"
  []
  #------ Conservation of Momentum (x-component) -------
  [momentum_x_time]
    type = FVTimeKernel
    variable = rhou
  []
  [momentum_x_advection]
    type = PCNSFVKT
    variable = rhou
    momentum_component = x
    eqn = "momentum"
  []
  [viscosity_x]
    type = FVOrthogonalDiffusion
    variable = rhou
    coeff = 'mu_air'
  []
  #------ Conservation of Momentum (y-component) -------
  [momentum_y_time]
    type = FVTimeKernel
    variable = rhov
  []
  [momentum_y_advection]
    type = PCNSFVKT
    variable = rhov
    momentum_component = y
    eqn = "momentum"
  []
  [viscosity_y]
    type = FVOrthogonalDiffusion
    variable = rhov
    coeff = 'mu_air'
  []
  [y_gravity] # Gravity only acts on the y-comp.
    type = PCNSFVMomentumGravity
    variable = rhov
    momentum_component = y
  []
  #------ Conservation of Momentum (z-component) -------
  [momentum_z_time]
    type = FVTimeKernel
    variable = rhow
  []
  [momentum_z_advection]
    type = PCNSFVKT
    variable = rhow
    momentum_component = z
    eqn = "momentum"
  []
  [viscosity_z]
    type = FVOrthogonalDiffusion
    variable = rhow
    coeff = 'mu_air'
  []
  # ------ Conservation of Energy  (Fluid Regions)-------
  [fluid_energy_time]
    type = FVTimeKernel
    variable = rho_et
  []
  [fluid_energy_advection]
    type = PCNSFVKT
    variable = rho_et
    eqn = "energy"
  []
  [fluid_conduction]
    type = FVOrthogonalDiffusion
    variable = rho_et
    coeff = 'k_air'
  []
[]
# ------------------------------------------------------------------------------
# AuxVariables
# ------------------------------------------------------------------------------
[AuxVariables]
  [Ma]
    order = CONSTANT
    family = MONOMIAL
    fv = true
    block = '1 2'
  []
  [pressure]
    order = CONSTANT
    family = MONOMIAL
    fv = true
    block = '1 2'
  []
  [U]
    order = CONSTANT
    family = MONOMIAL
    fv = true
    block = '1 2'
  []
  [temperature]
    order = CONSTANT
    family = MONOMIAL
    fv = true
    block = '1 2'
  []
  [porosity]
    type = MooseVariableFVReal
  []
[]
# ------------------------------------------------------------------------------
# AuxKernels
# ------------------------------------------------------------------------------
[AuxKernels]
  [Ma_aux]
    type = NSMachAux
    variable = Ma
    fluid_properties = fp
    use_material_properties = true
  []
  [p_aux]
    type = ADMaterialRealAux
    variable = pressure
    property = pressure
  []
  [U_aux]
    type = ADMaterialRealAux
    variable = U
    property = speed
  []
  [temperature_aux]
    type = ADMaterialRealAux
    variable = temperature
    property = T_fluid
  []
  [porosity]
    type = MaterialRealAux
    variable = porosity
    property = porosity
    execute_on = 'TIMESTEP_END'
  []
[]
# ------------------------------------------------------------------------------
# Initial Conditions
# ------------------------------------------------------------------------------
[ICs]
  [RPV_P_IC] # Initial pressure in RPV
    type = NSFunctionInitialCondition
    initial_pressure = ${P_RPV}
    initial_temperature = ${T_RPV}
    initial_velocity = '0 -1.0e-12 0'
    fluid_properties = fp
    variable = 'pressure'
    block = 1
  []
  [Cavity_P_IC] # Initial pressure in cavity
    type = NSFunctionInitialCondition
    initial_pressure = ${P_Cav}
    initial_temperature = ${T_Cav}
    initial_velocity = '0 -1.0e-12 0'
    fluid_properties = fp
    variable = 'pressure'
    block = 2
  []
  [RPV_T_IC] # Initial temperature in RPV
    # type = NSInitialCondition
    type = NSFunctionInitialCondition
    initial_pressure = ${P_RPV}
    initial_temperature = ${T_RPV}
    initial_velocity = '0 -1.0e-12 0'
    fluid_properties = fp
    variable = 'temperature'
    block = 1
  []
  [Cavity_T_IC] # Initial temperature in cavity
    type = NSFunctionInitialCondition
    initial_pressure = ${P_Cav}
    initial_temperature = ${T_Cav}
    initial_velocity = '0 -1.0e-12 0'
    fluid_properties = fp
    variable = 'temperature'
    block = 2
  []
  [RPV_rho_IC]
    type = NSFunctionInitialCondition
    initial_pressure = ${P_RPV}
    initial_temperature = ${T_RPV}
    initial_velocity = '0 -1.0e-12 0'
    fluid_properties = fp
    variable = 'rho'
    block = 1
  []
  [Cavity_rho_IC]
    type = NSFunctionInitialCondition
    initial_pressure = ${P_Cav}
    initial_temperature = ${T_Cav}
    initial_velocity = '0 -1.0e-12 0'
    fluid_properties = fp
    variable = 'rho'
    block = 2
  []
  [RPV_rhou_IC]
    type = NSFunctionInitialCondition
    initial_pressure = ${P_RPV}
    initial_temperature = ${T_RPV}
    initial_velocity = '0 -1.0e-12 0'
    fluid_properties = fp
    variable = 'rhou'
    block = 1
  []
  [Cavity_rhou_IC]
    type = NSFunctionInitialCondition
    initial_pressure = ${P_Cav}
    initial_temperature = ${T_Cav}
    initial_velocity = '0 -1.0e-12 0'
    fluid_properties = fp
    variable = 'rhou'
    block = 2
  []
  [RPV_rhov_IC]
    type = NSFunctionInitialCondition
    initial_pressure = ${P_RPV}
    initial_temperature = ${T_RPV}
    initial_velocity = '0 -1.0e-12 0'
    fluid_properties = fp
    variable = 'rhov'
    block = 1
  []
  [Cavity_rhov_IC]
    type = NSFunctionInitialCondition
    initial_pressure = ${P_Cav}
    initial_temperature = ${T_Cav}
    initial_velocity = '0 -1.0e-12 0'
    fluid_properties = fp
    variable = 'rhov'
    block = 2
  []
  [RPV_rhow_IC]
    type = NSFunctionInitialCondition
    initial_pressure = ${P_RPV}
    initial_temperature = ${T_RPV}
    initial_velocity = '0 -1.0e-12 0'
    fluid_properties = fp
    variable = 'rhow'
    block = 1
  []
  [Cavity_rhow_IC]
    type = NSFunctionInitialCondition
    initial_pressure = ${P_Cav}
    initial_temperature = ${T_Cav}
    initial_velocity = '0 -1.0e-12 0'
    fluid_properties = fp
    variable = 'rhow'
    block = 2
  []
  [RPV_rhoE_IC]
    type = NSFunctionInitialCondition
    initial_pressure = ${P_RPV}
    initial_temperature = ${T_RPV}
    initial_velocity = '0 -1.0e-12 0'
    fluid_properties = fp
    variable = 'rho_et'
    block = 1
  []
  [Cavity_rhoE_IC]
    type = NSFunctionInitialCondition
    initial_pressure = ${P_Cav}
    initial_temperature = ${T_Cav}
    initial_velocity = '0 -1.0e-12 0'
    fluid_properties = fp
    variable = 'rho_et'
    block = 2
  []
[]
# ------------------------------------------------------------------------------
# FVBCs
# ------------------------------------------------------------------------------
[FVBCs]
  [x_pressure_walls]
    type = PCNSFVImplicitMomentumPressureBC
    variable = rhou
    momentum_component = x
    boundary = 'RPVIW RPVOW PCV ContIW ContSP RPVSP'
  []
  [y_pressure_walls]
    type = PCNSFVImplicitMomentumPressureBC
    variable = rhov
    momentum_component = y
    boundary = 'RPVIW RPVOW PCV ContIW ContSP RPVSP'
  []
  [z_pressure_walls]
    type = PCNSFVImplicitMomentumPressureBC
    variable = rhow
    momentum_component = z
    boundary = 'RPVIW RPVOW PCV ContIW ContSP RPVSP'
  []
  [shear_x_walls]
    type = FVOrthogonalBoundaryDiffusion
    function = 0
    variable = rhou
    coeff = 'mu_air'
    diffusing_quantity = 'vel_x'
    boundary = 'RPVIW RPVOW PCV ContIW ContSP RPVSP'
  []
  [shear_y_walls]
    type = FVOrthogonalBoundaryDiffusion
    function = 0
    variable = rhov
    coeff = 'mu_air'
    diffusing_quantity = 'vel_y'
    boundary = 'RPVIW RPVOW PCV ContIW ContSP RPVSP'
  []
  [shear_z_walls]
    type = FVOrthogonalBoundaryDiffusion
    function = 0
    variable = rhow
    coeff = 'mu_air'
    diffusing_quantity = 'vel_z'
    boundary = 'RPVIW RPVOW PCV ContIW ContSP RPVSP'
  []
  # Venting section outlet
  [rho_outlet]
    type = PCNSFVStrongBC
    boundary = 'Outlet'
    variable = rho
    pressure = ${ambient_pressure}
    eqn = 'mass'
  []
  [rhou_outlet]
    type = PCNSFVStrongBC
    boundary = 'Outlet'
    variable = rhou
    pressure = ${ambient_pressure}
    eqn = 'momentum'
    momentum_component = x
  []
  [rhov_outlet]
    type = PCNSFVStrongBC
    boundary = 'Outlet'
    variable = rhov
    pressure = ${ambient_pressure}
    eqn = 'momentum'
    momentum_component = y
  []
  [rhow_outlet]
    type = PCNSFVStrongBC
    boundary = 'Outlet'
    variable = rhow
    pressure = ${ambient_pressure}
    eqn = 'momentum'
    momentum_component = z
  []
  [rhoet_outlet]
    type = PCNSFVStrongBC
    boundary = 'Outlet'
    variable = rho_et
    pressure = ${ambient_pressure}
    eqn = 'energy'
  []
  # [HT_Walls]
  #   type = FVNeumannBC
  #   boundary = 'RPVIW'
  #   variable = rho_et
  #   value = ${fparse 540.0/1.577917964}
  # []
  [RPV_walls]
    type = FVThermalResistanceBC
    geometry = cartesian
    T_ambient = ${T_RPV}
    thermal_conductivities = 16.3 # Thermal conductivity SS-316 [W/m-K]
    conduction_thicknesses = 0.0127 # Containment building thickness [m]
    htc = 'SS_HTC'
    emissivity = 0 # emissivity of SS-316 [-]
    boundary = 'RPVIW'
    variable = temperature
  []
  [cold_walls]
    type = FVThermalResistanceBC
    geometry = cartesian
    T_ambient = ${ambient_temperature}
    thermal_conductivities = 45 # Thermal conductivity carbon steel [W/m-K]
    conduction_thicknesses = 0.0047625 # Containment building thickness [m]
    htc = 'CS_HTC'
    emissivity = 0 # emissivity of carbon steel [-]
    boundary = 'ContIW'
    variable = temperature
  []
[]
# ------------------------------------------------------------------------------
# Materials and fluid properties
# ------------------------------------------------------------------------------
[Modules]
  [FluidProperties]
    [fp]
      type = IdealGasFluidProperties
      gamma = 1.4 # Heat capacity ratio [-]
      molar_mass = 28.97e-3 # [kg/mol] molar mass of air
    []
  []
[]
[Materials]
  [var_mat]
    type = PorousConservedVarMaterial
    rho = rho
    rho_et = rho_et
    superficial_rhou = rhou
    superficial_rhov = rhov
    superficial_rhow = rhow
    porosity = porosity
    block = '1 2'
  []
  [air_func]
    type = ADGenericFunctionMaterial
    prop_names = 'mu_air k_air cp_air'
    prop_values = '21.74e-6 0.03162 1010.7' # Heat capacity of air [J/kg-K]'
    block = '1 2'
  []
  [HTC_f_1]
    type = ADGenericConstantMaterial
    prop_names = 'CS_HTC'
    prop_values = 25 # airat transfer coefficient [W/m2-K]
    block = '1 2'
  []
  [HTC_f_2]
    type = ADGenericConstantMaterial
    prop_names = 'SS_HTC'
    prop_values = 50 # Heat transfer coefficient [W/m2-K]
    block = '1 2'
  []
  [porosity]
    type = GenericConstantMaterial
    prop_names = 'porosity'
    prop_values = 1
  []
[]
# ------------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------------
[Functions]
  [Cav_temp_func]
    type = ParsedFunction
    # Change of temperature in the y axis
    value = 309.15+(74.0216223*y)
  []
[]
# ------------------------------------------------------------------------------
# Preconditioning
# ------------------------------------------------------------------------------
[Preconditioning]
  [preco_SMP]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type'
    petsc_options_value = 'lu'
  []
[]
# ------------------------------------------------------------------------------
# Executioner
# ------------------------------------------------------------------------------
[Executioner]
  type = Transient
  # Time stepping parameters
  end_time = 25
  dtmax = 1
  dt = 0.0025
  dtmin = 1e-6
  # Solver parameters
  l_tol = 1e-6
  l_max_its = 50
  nl_max_its = 20
  # nl_rel_tol = 1e-8
  # nl_abs_tol = 1e-7
  automatic_scaling = true
  # steady_state_detection = false
  # steady_state_tolerance = 1e-10
[]
# ------------------------------------------------------------------------------
# Outputs
# ------------------------------------------------------------------------------
[Outputs]
  checkpoint = true
  exodus = true # Export exodus file
  csv = true # Export csv file with temp. and vel. values
  interval = 1  # only output every 40 timesteps
  print_nonlinear_converged_reason = true
  print_linear_residuals = true
  print_linear_converged_reason = true
[]
[Debug]
  show_var_residual_norms = true
[]
[Postprocessors]
  [TC01]
    type = PointValue
    point = '-0.060325 0.898525 0.23495'
    variable = temperature
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [TC02]
    type = PointValue
    point = '-0.460375 0.898525 0.23495'
    variable = temperature
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [TC03]
    type = PointValue
    point = '-0.060325 0.492125 0.23495'
    variable = temperature
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [TC04]
    type = PointValue
    point = '-0.460375 0.492125 0.23495'
    variable = temperature
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [TC05]
    type = PointValue
    point = '-0.060325 0.098425 0.23495'
    variable = temperature
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [TC06]
    type = PointValue
    point = '-0.460375 0.098425 0.23495'
    variable = temperature
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [TC20]
    type = PointValue
    point = '-0.028575 1.2573 0'
    variable = temperature
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [TC21]
    type = PointValue
    point = '-0.26035 1.2573 0.23098'
    variable = temperature
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [TC22]
    type = PointValue
    point = '-0.49133 1.2573 0'
    variable = temperature
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [TC23]
    type = PointValue
    point = '-0.62071 0.4457 0'
    variable = temperature
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [TC24]
    type = PointValue
    point = '-0.62071 0.21724 0'
    variable = temperature
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [TC26]
    type = PointValue
    point = '-1.171575 1.4224 0'
    variable = temperature
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [TC27]
    type = PointValue
    point = '-1.171575 0.4572 0.23495'
    variable = temperature
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [TC28]
    type = PointValue
    point = '-1.171575 -0.3556 0.23495'
    variable = temperature
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [TC31]
    type = PointValue
    point = '-1.222375 -0.30734 0'
    variable = temperature
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [TC32]
    type = PointValue
    point = '-1.222375 -0.30734 0.0508'
    variable = temperature
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [TC33]
    type = PointValue
    point = '-1.222375 -0.30734 0.1016'
    variable = temperature
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [VT08]
    type = PointValue
    point = '-0.58261 0.4457 0'
    variable = U
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [VT10]
    type = PointValue
    point = '-0.58261 0.21724 0'
    variable = U
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [VT13]
    type = PointValue
    point = '-1.222375 -0.30734 0.0254'
    variable = U
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [VT14]
    type = PointValue
    point = '-1.222375 -0.30734 0.0762'
    variable = U
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [Pressure_RPV]
    type = ElementAverageValue
    variable = pressure
    block = 1
  []
  [Pressure_Cav]
    type = ElementAverageValue
    variable = pressure
    block = 2
  []
  [Temperature_RPV]
    type = ElementAverageValue
    variable = temperature
    block = 1
  []
  [Temperature_Cav]
    type = ElementAverageValue
    variable = temperature
    block = 2
  []
[]
