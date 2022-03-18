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
ambient_temperature = ${fparse (14.77+273.15)} # Initial temperature [K]
# Containment ------------------------------------------------------------------
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
  [pressure]
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
  [T_fluid]
    type = MooseVariableFVReal
    block = '1 2'
  []
  [Mass_Fraction]
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
    type = FVMatPropTimeKernel
    variable = pressure
    mat_prop_time_derivative = 'dsuperficial_rho_dt'
  []
  [mass_advection]
    type = GasMixPCNSFVKT
    variable = pressure
    eqn = "mass"
  []
  #------ Conservation of Momentum (x-component) -------
  [momentum_x_time]
    type = FVMatPropTimeKernel
    variable = rhou
    mat_prop_time_derivative = 'dsuperficial_rhou_dt'
  []
  [momentum_x_advection]
    type = GasMixPCNSFVKT
    variable = rhou
    momentum_component = x
    eqn = "momentum"
  []
  # [viscosity_x]
  #   type = FVOrthogonalDiffusion
  #   variable = rhou
  #   coeff = 'muHe'
  # []
  #------ Conservation of Momentum (y-component) -------
  [momentum_y_time]
    type = FVMatPropTimeKernel
    variable = rhov
    mat_prop_time_derivative = 'dsuperficial_rhov_dt'
  []
  [momentum_y_advection]
    type = GasMixPCNSFVKT
    variable = rhov
    momentum_component = y
    eqn = "momentum"
  []
  # [viscosity_y]
  #   type = FVOrthogonalDiffusion
  #   variable = rhov
  #   coeff = 'muHe'
  # []
  [y_gravity] # Gravity only acts on the y-comp.
    type = PCNSFVMomentumGravity
    variable = rhov
    momentum_component = y
  []
  #------ Conservation of Momentum (z-component) -------
  [momentum_z_time]
    type = FVMatPropTimeKernel
    variable = rhow
    mat_prop_time_derivative = 'dsuperficial_rhow_dt'
  []
  [momentum_z_advection]
    type = GasMixPCNSFVKT
    variable = rhow
    momentum_component = z
    eqn = "momentum"
  []
  # [viscosity_z]
  #   type = FVOrthogonalDiffusion
  #   variable = rhow
  #   coeff = 'muHe'
  # []
  # ------ Conservation of Energy  (Fluid Regions)-------
  [fluid_energy_time]
    type = FVMatPropTimeKernel
    variable = T_fluid
    mat_prop_time_derivative = 'dsuperficial_rhoet_dt'
  []
  [fluid_energy_advection]
    type = GasMixPCNSFVKT
    variable = T_fluid
    eqn = "energy"
  []
  # [fluid_conduction]
  #   type = FVOrthogonalDiffusion
  #   variable = T_fluid
  #   coeff = k
  # []
  # --------------- Mass fraction advection ---------------
  [mass_frac_time]
    type = FVMatPropTimeKernel
    variable = Mass_Fraction
    mat_prop_time_derivative = 'drho_f_dt'
  []
  [MF_Advection]
    type = GasMixPCNSFVKT
    variable = Mass_Fraction
    eqn = scalar
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
  [rho]
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
  [rho_aux]
    type = ADMaterialRealAux
    variable = rho
    property = rho
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
    fluid_properties = fp_helium
    variable = 'pressure'
    block = 1
  []
  [Cavity_P_IC] # Initial pressure in cavity
    type = NSFunctionInitialCondition
    initial_pressure = ${P_Cav}
    initial_temperature = 'Cav_temp_func'
    initial_velocity = '0 -1.0e-12 0'
    fluid_properties = fp_air
    variable = 'pressure'
    block = 2
  []
  [RPV_T_IC] # Initial temperature in RPV
    # type = NSInitialCondition
    type = NSFunctionInitialCondition
    initial_pressure = ${P_RPV}
    initial_temperature = ${T_RPV}
    initial_velocity = '0 -1.0e-12 0'
    fluid_properties = fp_helium
    variable = 'temperature'
    block = 1
  []
  [Cavity_T_IC] # Initial temperature in cavity
    type = NSFunctionInitialCondition
    initial_pressure = ${P_Cav}
    initial_temperature = 'Cav_temp_func'
    initial_velocity = '0 -1.0e-12 0'
    fluid_properties = fp_air
    variable = 'temperature'
    block = 2
  []
  [RPV_rho_IC]
    type = NSFunctionInitialCondition
    initial_pressure = ${P_RPV}
    initial_temperature = ${T_RPV}
    initial_velocity = '0 -1.0e-12 0'
    fluid_properties = fp_helium
    variable = 'rho'
    block = 1
  []
  [Cavity_rho_IC]
    type = NSFunctionInitialCondition
    initial_pressure = ${P_Cav}
    initial_temperature = 'Cav_temp_func'
    initial_velocity = '0 -1.0e-12 0'
    fluid_properties = fp_air
    variable = 'rho'
    block = 2
  []
  [RPV_rhou_IC]
    type = NSFunctionInitialCondition
    initial_pressure = ${P_RPV}
    initial_temperature = ${T_RPV}
    initial_velocity = '0 -1.0e-12 0'
    fluid_properties = fp_helium
    variable = 'rhou'
    block = 1

  []
  [Cavity_rhou_IC]
    type = NSFunctionInitialCondition
    initial_pressure = ${P_Cav}
    initial_temperature = 'Cav_temp_func'
    initial_velocity = '0 -1.0e-12 0'
    fluid_properties = fp_air
    variable = 'rhou'
    block = 2
  []
  [RPV_rhov_IC]
    type = NSFunctionInitialCondition
    initial_pressure = ${P_RPV}
    initial_temperature = ${T_RPV}
    initial_velocity = '0 -1.0e-12 0'
    fluid_properties = fp_helium
    variable = 'rhov'
    block = 1
  []
  [Cavity_rhov_IC]
    type = NSFunctionInitialCondition
    initial_pressure = ${P_Cav}
    initial_temperature = 'Cav_temp_func'
    initial_velocity = '0 -1.0e-12 0'
    fluid_properties = fp_air
    variable = 'rhov'
    block = 2
  []
  [RPV_rhow_IC]
    type = NSFunctionInitialCondition
    initial_pressure = ${P_RPV}
    initial_temperature = ${T_RPV}
    initial_velocity = '0 -1.0e-12 0'
    fluid_properties = fp_helium
    variable = 'rhow'
    block = 1
  []
  [Cavity_rhow_IC]
    type = NSFunctionInitialCondition
    initial_pressure = ${P_Cav}
    initial_temperature = 'Cav_temp_func'
    initial_velocity = '0 -1.0e-12 0'
    fluid_properties = fp_air
    variable = 'rhow'
    block = 2
  []
  # [RPV_rhoE_IC]
  #   type = NSFunctionInitialCondition
  #   initial_pressure = ${P_RPV}
  #   initial_temperature = ${T_RPV}
  #   initial_velocity = '0 -1.0e-12 0'
  #   fluid_properties = fp_helium
  #   variable = 'rho_et'
  #   block = 1
  # []
  # [Cavity_rhoE_IC]
  #   type = NSFunctionInitialCondition
  #   initial_pressure = ${P_Cav}
  #   initial_temperature = 'Cav_temp_func'
  #   initial_velocity = '0 -1.0e-12 0'
  #   fluid_properties = fp_air
  #   variable = 'rho_et'
  #   block = 2
  # []
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
  # [shear_x_walls]
  #   type = FVOrthogonalBoundaryDiffusion
  #   function = 0
  #   variable = rhou
  #   coeff = 'muHe'
  #   diffusing_quantity = 'vel_x'
  #   boundary = 'RPVIW RPVOW PCV ContIW ContSP RPVSP'
  # []
  # [shear_y_walls]
  #   type = FVOrthogonalBoundaryDiffusion
  #   function = 0
  #   variable = rhov
  #   coeff = 'muHe'
  #   diffusing_quantity = 'vel_y'
  #   boundary = 'RPVIW RPVOW PCV ContIW ContSP RPVSP'
  # []
  # [shear_z_walls]
  #   type = FVOrthogonalBoundaryDiffusion
  #   function = 0
  #   variable = rhow
  #   coeff = 'muHe'
  #   diffusing_quantity = 'vel_z'
  #   boundary = 'RPVIW RPVOW PCV ContIW ContSP RPVSP'
  # []
  # Venting section outlet
  [rho_outlet]
    type = GasMixPCNSFVStrongBC
    boundary = 'Outlet'
    variable = pressure
    pressure = ${ambient_pressure}
    eqn = 'mass'
  []
  [rhou_outlet]
    type = GasMixPCNSFVStrongBC
    boundary = 'Outlet'
    variable = rhou
    pressure = ${ambient_pressure}
    eqn = 'momentum'
    momentum_component = x
  []
  [rhov_outlet]
    type = GasMixPCNSFVStrongBC
    boundary = 'Outlet'
    variable = rhov
    pressure = ${ambient_pressure}
    eqn = 'momentum'
    momentum_component = y
  []
  [rhow_outlet]
    type = GasMixPCNSFVStrongBC
    boundary = 'Outlet'
    variable = rhow
    pressure = ${ambient_pressure}
    eqn = 'momentum'
    momentum_component = z
  []
  [rhoet_outlet]
    type = GasMixPCNSFVStrongBC
    boundary = 'Outlet'
    variable = T_fluid
    pressure = ${ambient_pressure}
    eqn = 'energy'
  []
  [Mass_Fraction_Outlet]
    type = GasMixPCNSFVStrongBC
    boundary = 'Outlet'
    variable = Mass_Fraction
    pressure = ${ambient_pressure}
    eqn = 'scalar'
  []
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
    [fp_helium]
      type = IdealGasFluidProperties
      gamma = ${fparse 5/3.} # Heat capacity ratio [-]
      molar_mass = 4.002602e-3 # [kg/mol] molar mass of helium
      # k = 0.13754
      # mu = 1.76e-05
    []
    [fp_air]
      type = IdealGasFluidProperties
      gamma = ${fparse 7/5.} # Heat capacity ratio [-]
      molar_mass = 28.97e-3 # [kg/mol] molar mass of air
      # k = 0.02241
      # mu = 1.606e-05
    []
    [fp]
      type = GasMixPHFluidProperties
      fp_primary = fp_helium
      fp_secondary = 'fp_air'
    []
  []
[]
[Materials]
  [var_mat]
    type = GasMixPorousMixedVarMaterial
    fp = fp
    pressure = pressure
    T_fluid = temperature
    superficial_rhou = rhou
    superficial_rhov = rhov
    superficial_rhow = rhow
    secondary_fraction = Mass_Fraction
    porosity = porosity
    block = '1 2'
  []
  # [var_mat]
  #   type = PorousConservedVarMaterial
  #   rho = rho
  #   rho_et = rho_et
  #   superficial_rhou = rhou
  #   superficial_rhov = rhov
  #   superficial_rhow = rhow
  #   porosity = porosity
  #   block = '1 2'
  # []
  # [he_func]
  #   type = ADGenericFunctionMaterial
  #   prop_names = 'muHe kHe cpHe'
  #   prop_values = 'muHe kHe 5190.' # Heat capacity of He [J/kg-K]'
  #   block = '1 2'
  # []
  # [air_func]
  #   type = ADGenericFunctionMaterial
  #   prop_names = 'muAir kAir cpHe'
  #   prop_values = 'muAir kAir 1007.48' # Average heat capacity of air [J/kg-K]
  #   block = '1 2'
  # []
  [HTC_f_1]
    type = ADGenericConstantMaterial
    prop_names = 'CS_HTC'
    prop_values = 25 # Heat transfer coefficient [W/m2-K]
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
  # [kHe]
  #   type = PiecewiseLinear # x in Kelvin and y in W/m-K
  #   x = '250	270	290	310	330	350	370	390	410	430	450	470	490	510	530	550	570	590	610	630	650'
  #   y = '0.13754	0.14503	0.15236	0.15955	0.1666	0.17353	0.18034	0.18705	0.19366	0.20017	0.20659	0.21294	0.2192	0.22539	0.2315	0.23755	0.24354	0.24946	0.25533	0.26113	0.26689'
  # []
  # [muHe]
  #   type = PiecewiseLinear # x in Kelvin and y in Pa-s
  #   x = '250	270	290	310	330	350	370	390	410	430	450	470	490	510	530	550	570	590	610	630	650'
  #   y = '1.76e-05	1.85e-05	1.95e-05	2.04e-05	2.13e-05	2.22e-05	2.30e-05	2.39e-05	2.47e-05	2.55e-05	2.64e-05	2.72e-05	2.80e-05	2.88e-05	2.95e-05	3.03e-05	3.11e-05	3.18e-05	3.26e-05	3.33e-05	3.41e-05'
  # []
  # [kAir]
  #   type = PiecewiseLinear # x in Kelvin and y in W/m-K
  #   x = '250 300 350 400 450 500 550 600 650'
  #   y = '0.02241 0.02623 0.02984 0.03328 0.03656 0.03971 0.04277 0.04573 0.04863'
  # []
  # [muAir]
  #   type = PiecewiseLinear # x in Kelvin and y in Pa-s
  #   x = '250 300 350 400 450 500 550 600 650'
  #   y = '1.606e-05 1.857e-05 2.09e-05 2.31e-05 2.517e-05 2.713e-05 2.902e-05 3.082e-05 3.257e-05'
  # []
  [RPV_temp_func]
    type = ParsedFunction
    # Change of temperature in the y axis
    value = 272.619607+(23.78578*y)
  []
  [Cav_temp_func]
    type = ParsedFunction
    # Change of temperature in the y axis
    value = 292.9447+(23.785781*y)
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
  end_time = 25
  dtmax = 1
  dt = 0.0025
  dtmin = 1e-6
  l_tol = 1e-6
  l_max_its = 50
  nl_max_its = 25
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
  interval = 40  # only output every 40 timesteps
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
