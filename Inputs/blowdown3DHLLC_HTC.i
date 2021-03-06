# ==============================================================================
# Blowdown of a 1/20th Scaled-down HTGR UI
# ------------------------------------------------------------------------------
# Idaho Falls, INL, July 23, 2021 8:00 AM
# Author(s): Silvino A. Balderrama Prieto, Guillaume L. Giudicelli
# ------------------------------------------------------------------------------
# Description: The following model simulates the blowdown of a scaled-down HTGR
# as a result of a break in the RPV. The porous media compressible finite volume
# HLLC is implemented to solve for the transient depressurization.
# ==============================================================================
# MODEL PARAMETERS
# ------------------------------------------------------------------------------
ambient_pressure = 101325 # Ambient pressure [Pa]
ambient_temperature = ${fparse (14.77+273.15)} # Initial temperature [K]
# Containment ------------------------------------------------------------------
T_Cav = ${fparse (37.2+273.15)} # Initial temperature [K]
P_Cav = 101325 # Initial pressure [Pa]

# RPV --------------------------------------------------------------------------
T_RPV = ${fparse 125.7718+273.15} # Initial temperature [K]
P_RPV = ${fparse 161.31415*6894.75729} # Initial pressure [Pa]
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
    file = './HTGR3D_B01_V12.e'
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
  # ------ Conservation of Mass -------
  [mass_time]
    type = FVTimeKernel
    variable = rho
  []
  [mass_advection]
    type = PCNSFVMassHLLC
    variable = rho
  []
  #------ Conservation of Momentum (x-component) -------
  [momentum_x_time]
    type = FVTimeKernel
    variable = rhou
  []
  [momentum_x_advection]
    type = PCNSFVMomentumHLLC
    variable = rhou
    momentum_component = x
  []
  [viscosity_x]
    type = FVOrthogonalDiffusion
    variable = rhou
    coeff = 'muHe'
  []
  #------ Conservation of Momentum (y-component) -------
  [momentum_y_time]
    type = FVTimeKernel
    variable = rhov
  []
  [momentum_y_advection]
    type = PCNSFVMomentumHLLC
    variable = rhov
    momentum_component = y
  []
  [viscosity_y]
    type = FVOrthogonalDiffusion
    variable = rhov
    coeff = 'muHe'
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
    coeff = 'muHe'
  []
  # ------ Conservation of Energy  (Fluid Regions)-------
  [fluid_energy_time]
    type = FVTimeKernel
    variable = rho_et
  []
  [fluid_energy_advection]
    type = PCNSFVFluidEnergyHLLC
    variable = rho_et
  []
  [fluid_conduction]
    type = FVOrthogonalDiffusion
    variable = rho_et
    coeff = 'kHe'
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
    type = NSInitialCondition
    initial_pressure = ${P_RPV}
    initial_temperature = ${T_RPV}
    initial_velocity = '0 -1.0e-12 0'
    fluid_properties = fp
    variable = 'pressure'
    block = 1
  []
  [Cavity_P_IC] # Initial pressure in cavity
    type = NSInitialCondition
    initial_pressure = ${P_Cav}
    initial_temperature = ${T_Cav}
    initial_velocity = '0 -1.0e-12 0'
    fluid_properties = fp
    variable = 'pressure'
    block = 2
  []
  [RPV_T_IC] # Initial temperature in RPV
    type = NSInitialCondition
    initial_pressure = ${P_RPV}
    initial_temperature = ${T_RPV}
    initial_velocity = '0 -1.0e-12 0'
    fluid_properties = fp
    variable = 'temperature'
    block = 1
  []
  [Cavity_T_IC] # Initial temperature in cavity
    type = NSInitialCondition
    initial_pressure = ${P_Cav}
    initial_temperature = ${T_Cav}
    initial_velocity = '0 -1.0e-12 0'
    fluid_properties = fp
    variable = 'temperature'
    block = 2
  []
  [RPV_rho_IC]
    type = NSInitialCondition
    initial_pressure = ${P_RPV}
    initial_temperature = ${T_RPV}
    initial_velocity = '0 -1.0e-12 0'
    fluid_properties = fp
    variable = 'rho'
    block = 1
  []
  [Cavity_rho_IC]
    type = NSInitialCondition
    initial_pressure = ${P_Cav}
    initial_temperature = ${T_Cav}
    initial_velocity = '0 -1.0e-12 0'
    fluid_properties = fp
    variable = 'rho'
    block = 2
  []
  [RPV_rhou_IC]
    type = NSInitialCondition
    initial_pressure = ${P_RPV}
    initial_temperature = ${T_RPV}
    initial_velocity = '0 -1.0e-12 0'
    fluid_properties = fp
    variable = 'rhou'
    block = 1
  []
  [Cavity_rhou_IC]
    type = NSInitialCondition
    initial_pressure = ${P_Cav}
    initial_temperature = ${T_Cav}
    initial_velocity = '0 -1.0e-12 0'
    fluid_properties = fp
    variable = 'rhou'
    block = 2
  []
  [RPV_rhov_IC]
    type = NSInitialCondition
    initial_pressure = ${P_RPV}
    initial_temperature = ${T_RPV}
    initial_velocity = '0 -1.0e-12 0'
    fluid_properties = fp
    variable = 'rhov'
    block = 1
  []
  [Cavity_rhov_IC]
    type = NSInitialCondition
    initial_pressure = ${P_Cav}
    initial_temperature = ${T_Cav}
    initial_velocity = '0 -1.0e-12 0'
    fluid_properties = fp
    variable = 'rhov'
    block = 2
  []
  [RPV_rhow_IC]
    type = NSInitialCondition
    initial_pressure = ${P_RPV}
    initial_temperature = ${T_RPV}
    initial_velocity = '0 -1.0e-12 0'
    fluid_properties = fp
    variable = 'rhow'
    block = 1
  []
  [Cavity_rhow_IC]
    type = NSInitialCondition
    initial_pressure = ${P_Cav}
    initial_temperature = ${T_Cav}
    initial_velocity = '0 -1.0e-12 0'
    fluid_properties = fp
    variable = 'rhow'
    block = 2
  []
  [RPV_rhoE_IC]
    type = NSInitialCondition
    initial_pressure = ${P_RPV}
    initial_temperature = ${T_RPV}
    initial_velocity = '0 -1.0e-12 0'
    fluid_properties = fp
    variable = 'rho_et'
    block = 1
  []
  [Cavity_rhoE_IC]
    type = NSInitialCondition
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
    coeff = 'muHe'
    diffusing_quantity = 'vel_x'
    boundary = 'RPVIW RPVOW PCV ContIW ContSP RPVSP'
  []
  [shear_y_walls]
    type = FVOrthogonalBoundaryDiffusion
    function = 0
    variable = rhov
    coeff = 'muHe'
    diffusing_quantity = 'vel_y'
    boundary = 'RPVIW RPVOW PCV ContIW ContSP RPVSP'
  []
  [shear_z_walls]
    type = FVOrthogonalBoundaryDiffusion
    function = 0
    variable = rhow
    coeff = 'muHe'
    diffusing_quantity = 'vel_z'
    boundary = 'RPVIW RPVOW PCV ContIW ContSP RPVSP'
  []
  # Venting section outlet
  [rho_outlet]
    type = PCNSFVStrongBC
    boundary = 'Outlet'
    variable = rho
    T_fluid = ${ambient_temperature}
    pressure = ${ambient_pressure}
    eqn = 'mass'
  []
  [rhou_outlet]
    type = PCNSFVStrongBC
    boundary = 'Outlet'
    variable = rhou
    T_fluid = ${ambient_temperature}
    pressure = ${ambient_pressure}
    eqn = 'momentum'
    momentum_component = x
  []
  [rhov_outlet]
    type = PCNSFVStrongBC
    boundary = 'Outlet'
    variable = rhov
    T_fluid = ${ambient_temperature}
    pressure = ${ambient_pressure}
    eqn = 'momentum'
    momentum_component = y
  []
  [rhow_outlet]
    type = PCNSFVStrongBC
    boundary = 'Outlet'
    variable = rhow
    T_fluid = ${ambient_temperature}
    pressure = ${ambient_pressure}
    eqn = 'momentum'
    momentum_component = z
  []
  [rhoet_outlet]
    type = PCNSFVStrongBC
    boundary = 'Outlet'
    variable = rho_et
    T_fluid = ${ambient_temperature}
    pressure = ${ambient_pressure}
    eqn = 'energy'
  []
  [cold_walls]
    type = FVThermalResistanceBC
    geometry = cartesian
    T_ambient = ${ambient_temperature}
    thermal_conductivities = 45 # Thermal conductivity carbon steel [W/m-K]
    conduction_thicknesses = 0.0047625 # Containment building thickness [m]
    htc = 'HTC'
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
      gamma = ${fparse 5/3.} # Heat capacity ratio [-]
      molar_mass = 4.002602e-3 # [kg/mol] molar mass of helium
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
    porosity = porosity
    block = '1 2'
  []
  [he_func]
    type = ADGenericFunctionMaterial
    prop_names = 'muHe kHe cpHe'
    prop_values = 'muHe kHe 5190.' # Heat capacity of He [J/kg-K]'
    block = '1 2'
  []
  [fHTC]
    type = ADGenericConstantMaterial
    prop_names = 'HTC'
    prop_values = 15 # Heat transfer coefficient [W/m2-K]
    block = '1 2'
  []
  [sound_speed]
    type = SoundspeedMat
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
  [kHe]
    type = PiecewiseLinear # x in Kelvin and y in W/m-K
    x = '250	270	290	310	330	350	370	390	410	430	450	470	490	510	530	550	570	590	610	630	650'
    y = '0.13754	0.14503	0.15236	0.15955	0.1666	0.17353	0.18034	0.18705	0.19366	0.20017	0.20659	0.21294	0.2192	0.22539	0.2315	0.23755	0.24354	0.24946	0.25533	0.26113	0.26689'
  []
  [muHe]
    type = PiecewiseLinear # x in Kelvin and y in Pa-s
    x = '250	270	290	310	330	350	370	390	410	430	450	470	490	510	530	550	570	590	610	630	650'
    y = '1.76e-05	1.85e-05	1.95e-05	2.04e-05	2.13e-05	2.22e-05	2.30e-05	2.39e-05	2.47e-05	2.55e-05	2.64e-05	2.72e-05	2.80e-05	2.88e-05	2.95e-05	3.03e-05	3.11e-05	3.18e-05	3.26e-05	3.33e-05	3.41e-05'
  []
  # [temp_func]
  #   type = ParsedFunction
  #   # Bottom at -0.446088m, top at 1.439685 m
  #   value = 288.75 + 64.16 * (y+0.446088)/(1.885773)
  # []
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
  end_time = 150
  dtmax = 1
  dt = 0.005
  dtmin = 1e-5
  l_tol = 1e-6
  l_max_its = 50
  nl_max_its = 25
  # nl_rel_tol = 1e-8
  # nl_abs_tol = 1e-7
  automatic_scaling = true
  # steady_state_detection = false
  # steady_state_tolerance = 1e-10
  # [TimeStepper]
  #   type = PostprocessorDT
  #   postprocessor = cfl_dt
  #   dt = 0.1
    # cutback_factor_at_failure = 0.5
    # reset_dt = true
  # []
[]
# ------------------------------------------------------------------------------
# Outputs
# ------------------------------------------------------------------------------
[Outputs]
  checkpoint = true
  exodus = true # Export exodus file
  csv = true # Export csv file with temp. and vel. values
  interval = 40  # only output every 40 timesteps
[]
[Postprocessors]
  # [cfl_dt]
  #   type = ADCFLTimeStepSize
  #   c_names = 'sound_speed'
  #   vel_names = 'speed'
  # []
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
    point = '-1.222375 1.235075 0'
    variable = temperature
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [TC32]
    type = PointValue
    point = '-1.222375 1.235075 0.0508'
    variable = temperature
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [TC33]
    type = PointValue
    point = '-1.222375 1.235075 0.1016'
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
    point = '-1.222375 1.235075 0.0254'
    variable = U
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [VT14]
    type = PointValue
    point = '-1.222375 1.235075 0.0762'
    variable = U
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [Pressure_RPV]
    type = ElementAverageValue
    variable = pressure
    block = 1
  []
  [PT-01]
    type = PointValue
    variable = pressure
    point = '-0.260355 0.02 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]
