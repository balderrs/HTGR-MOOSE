# ==============================================================================
# 2D Blowdown of a 1/20th Scaled-down HTGR UI
# ------------------------------------------------------------------------------
# Idaho Falls, INL, July 23, 2021 8:00 AM
# Author(s): Silvino A. Balderrama Prieto, Guillaume L. Giudicelli
# ------------------------------------------------------------------------------
# Description: The following model simulates the blowdown of a 2D scaled-down HTGR
# as a result of a break in the RPV. The porous media compressible finite volume
# KT is implemented to solve for the transient depressurization.
# ==============================================================================
# MODEL PARAMETERS

# Ambient ----------------------------------------------------------------------
ambient_pressure = 101325 # Ambient pressure [Pa]
ambient_temperature = ${fparse (16.2+273.15)} # Ambient temperature [K]

# Containment ------------------------------------------------------------------
# T_Cav = ${fparse (19+273.15)} # Initial temperature [K]
P_Cav = 101325 # Initial pressure [Pa]

# RPV --------------------------------------------------------------------------
T_RPV = ${fparse 125.7718+273.15} # Initial temperature [K]
P_RPV = ${fparse 146.61825*6894.75729} # Initial pressure [Pa]
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
    file = './HTGR2D_B01_V12.e'
  []
[]
# ------------------------------------------------------------------------------
# Variables
# ------------------------------------------------------------------------------
[Variables]
  [pressure]
    type = MooseVariableFVReal
    block = '1 2'
  []
  [temperature]
    type = MooseVariableFVReal
    block = '1 2'
  []
  [sup_vel_x]
    type = MooseVariableFVReal
    block = '1 2'
  []
  [sup_vel_y]
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
    mat_prop_time_derivative = 'dsuperficial_rho_dt'
    variable = pressure
  []
  [mass_advection]
    type = PCNSFVKT
    variable = pressure
    eqn = "mass"
  []
  #------ Conservation of Momentum (x-component) -------
  [momentum_x_time]
    type = FVTimeKernel
    mat_prop_time_derivative = 'dsuperficial_rhou_dt'
    variable = sup_vel_x
  []
  [momentum_x_advection]
    type = PCNSFVKT
    variable = sup_vel_x
    momentum_component = x
    eqn = "momentum"
  []
  [viscosity_x]
    type = FVOrthogonalDiffusion
    variable = sup_vel_x
    diffusing_quantity = vel_x
    coeff = 'muHe'
  []
  #------ Conservation of Momentum (y-component) -------
  [momentum_y_time]
    type = FVTimeKernel
    mat_prop_time_derivative = 'dsuperficial_rhov_dt'
    variable = sup_vel_y
  []
  [momentum_y_advection]
    type = PCNSFVKT
    variable = sup_vel_y
    momentum_component = y
    eqn = "momentum"
  []
  [viscosity_y]
    type = FVOrthogonalDiffusion
    variable = sup_vel_y
    diffusing_quantity = vel_y
    coeff = 'muHe'
  []
  [y_gravity] # Gravity only acts on the y-comp.
    type = PCNSFVMomentumGravity
    variable = sup_vel_y
    gravity = '0 -9.81 0'
    momentum_component = y
  []
  # ------ Conservation of Energy  (Fluid Regions)-------
  [fluid_energy_time]
    type = FVTimeKernel
    mat_prop_time_derivative = 'dsuperficial_rho_et_dt'
    variable = temperature
  []
  [fluid_energy_advection]
    type = PCNSFVKT
    variable = temperature
    eqn = "energy"
  []
  [fluid_conduction]
    type = FVOrthogonalDiffusion
    variable = temperature
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
  []
  [U]
    type = MooseVariableFVReal
  []
  [vel_x]
    type = MooseVariableFVReal
  []
  [vel_y]
    type = MooseVariableFVReal
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
  [U_aux]
    type = ADMaterialRealAux
    variable = U
    property = speed
  []
  [vel_x]
    type = ADMaterialRealAux
    variable = vel_x
    property = vel_x
    execute_on = 'TIMESTEP_END'
  []
  [vel_y]
    type = ADMaterialRealAux
    variable = vel_y
    property = vel_y
    execute_on = 'TIMESTEP_END'
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
  [RPV_T_IC]
    type = ConstantIC
    variable = temperature
    value = ${T_RPV}
    block = 1
  []
  [Cavity_T_IC]
    type = FunctionIC
    variable = temperature
    function = 'temp_func'
    block = 2
  []
  [RPV_P_IC]
    type = ConstantIC
    variable = pressure
    value = ${P_RPV}
    block = 1
  []
  [Cavity_P_IC]
    type = ConstantIC
    variable = pressure
    value = ${P_Cav}
    block = 2
  []
[]
# ------------------------------------------------------------------------------
# FVBCs
# ------------------------------------------------------------------------------
[FVBCs]
  # Walls
  [pressure_x_walls]
    type = PCNSFVImplicitMomentumPressureBC
    variable = sup_vel_x
    momentum_component = x
    boundary = 'RPVIW RPVOW PCU ContIW'
  []
  [pressure_y_walls]
    type = PCNSFVImplicitMomentumPressureBC
    variable = sup_vel_y
    momentum_component = y
    boundary = 'RPVIW RPVOW PCU ContIW'
  []
  [shear_x_walls]
    type = FVOrthogonalBoundaryDiffusion
    function = 0
    variable = sup_vel_x
    diffusing_quantity = vel_x
    coeff = 'muHe'
    boundary = 'RPVIW RPVOW PCU ContIW'
  []
  [shear_y_walls]
    type = FVOrthogonalBoundaryDiffusion
    function = 0
    variable = sup_vel_y
    diffusing_quantity = vel_y
    coeff = 'muHe'
    boundary = 'RPVIW RPVOW PCU ContIW'
  []
  # Venting section outlet
  [rho_outlet]
    type = PCNSFVStrongBC
    boundary = 'Outlet'
    variable = pressure
    T_fluid = ${ambient_temperature}
    pressure = ${ambient_pressure}
    eqn = 'mass'
  []
  [rhou_outlet]
    type = PCNSFVStrongBC
    boundary = 'Outlet'
    variable = sup_vel_x
    T_fluid = ${ambient_temperature}
    pressure = ${ambient_pressure}
    eqn = 'momentum'
    momentum_component = x
  []
  [rhov_outlet]
    type = PCNSFVStrongBC
    boundary = 'Outlet'
    variable = sup_vel_y
    T_fluid = ${ambient_temperature}
    pressure = ${ambient_pressure}
    eqn = 'momentum'
    momentum_component = y
  []
  [rho_et_outlet]
    type = PCNSFVStrongBC
    boundary = 'Outlet'
    variable = temperature
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
    type = PorousMixedVarMaterial
    pressure = pressure
    T_fluid = temperature
    superficial_rhou = sup_vel_x
    superficial_rhov = sup_vel_y
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
  [temp_func]
    type = PiecewiseLinear
    x = '-446 98.425 1275 1475'
    y = '291.15 298.45 399 401'
    axis = y
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
  end_time = 150
  l_tol = 1e-7
  automatic_scaling = true
  [./TimeStepper]
    type = PostprocessorDT
    postprocessor = cfl_dt
  [../]
  steady_state_detection = false
  steady_state_tolerance = 1e-10
[]
# ------------------------------------------------------------------------------
# Outputs
# ------------------------------------------------------------------------------
[Outputs]
  exodus = true # Export exodus file
  csv = true # Export csv file with temp. and vel. values
  interval = 50  # only output every 10 timesteps
[]
[Postprocessors]
  [cfl_dt]
    type = ADCFLTimeStepSize
    c_names = 'sound_speed'
    vel_names = 'speed'
  []
  [TC01]
    type = PointValue
    point = '-0.068262 0.827087 0'
    variable = temperature
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [TC02]
    type = PointValue
    point = '-0.509587 0.827087 0'
    variable = temperature
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [TC03]
    type = PointValue
    point = '-0.068262 0.496887 0'
    variable = temperature
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [TC04]
    type = PointValue
    point = '-0.509587 0.496887 0'
    variable = temperature
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [TC06]
    type = PointValue
    point = '-0.509587 0.103187 0'
    variable = temperature
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [TC20]
    type = PointValue
    point = '-0.0341315 1.267 0'
    variable = temperature
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [TC21]
    type = PointValue
    point = '-0.2651125 1.267 0'
    variable = temperature
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [TC22]
    type = PointValue
    point = '-0.4960945 1.267 0'
    variable = temperature
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [TC23]
    type = PointValue
    point = '-0.582613 0.437833 0'
    variable = temperature
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [TC24]
    type = PointValue
    point = '-0.582613 0.234633 0'
    variable = temperature
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [TC25]
    type = PointValue
    point = '-0.847638 1.427163 0'
    variable = temperature
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [TC26]
    type = PointValue
    point = '-1.035138 1.427163 0'
    variable = temperature
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [TC31]
    type = PointValue
    point = '-1.227138 1.2398375 0'
    variable = temperature
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [VT08]
    type = PointValue
    point = '-0.611187 0.437833 0'
    variable = U
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [VT10]
    type = PointValue
    point = '-0.611187 0.234633 0'
    variable = U
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [VT13]
    type = PointValue
    point = '-1.252537 1.239837 0'
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
    point = '-0.2651125 0.023 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]
