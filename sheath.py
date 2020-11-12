##  Hi, man! Never give up!  ##
import math as m
Te        = 2.50e+01
Ti        = 2.50e+01
Density0  = 6.05e+17
B0        = 5.00e-01
Theta     = 6.00e+01
Z         = 1.00e+00
A         = 1.00e+00
MiMe      = 1836
PPC       = 1.00e+03
Dx        = 5.00e-01
XCELL     = 512
Nperod    = 5.00e-01
Dump      = 100

### --------------- Calculate the plasma paramters ---------------###
FunderCharge    = 1.602176565e-19
EleMass         = 9.109382616e-31
LightSpeed      = 299792458
Epsilon         = 8.854187818e-12
Wpe             = m.sqrt((Density0*FunderCharge**2)/(EleMass*Epsilon))
Wpi             = m.sqrt((Density0*FunderCharge**2)/(EleMass*MiMe*Epsilon*A))
Wce             = FunderCharge*B0/EleMass
Wci             = Z*FunderCharge*B0/EleMass/MiMe/A
Vte             = m.sqrt(Te*FunderCharge/EleMass)
Vti             = m.sqrt(Te*FunderCharge/EleMass/MiMe)
Re              = Vte/Wce
Ri              = Vti/Wci
De              = Vte/Wpe
Di              = Vti/Wpi
### -------------- Normalization -----------------###
Te              = Te/511.0e3
Ti              = Ti/511.0e3
N0              = Density0/(Epsilon*EleMass*Wpe**2/FunderCharge**2) 
B0              = B0/(EleMass*Wpe/FunderCharge)
Dx              = Dx*De/(LightSpeed/Wpe)
Tsim           = int(Nperod*2*m.pi/Wci*Wpe)+10

### ------------- Setting for Smilei ----------------###

Main(
    geometry = "1Dcartesian",
    interpolation_order = 2,
    timestep_over_CFL = 0.95,
    simulation_time = Tsim,
    cell_length  = [Dx],
    number_of_cells  = [XCELL],
    number_of_patches = [64],
    EM_boundary_conditions = [ ["reflective"] ],
    print_every = 1,
    random_seed = smilei_mpi_rank
)
fp = trapezoidal(1., xvacuum=0.0, xplateau=Dx*XCELL)
fm = trapezoidal(1., xvacuum=0.0, xplateau=(XCELL-500)*Dx)

Species(
    name = "ion",
    position_initialization = "regular",
    momentum_initialization = "mj",
    particles_per_cell = PPC, 
    c_part_max = 1.0,
    mass = MiMe*A,
    charge = Z,
    number_density = N0,
    mean_velocity = [0.0, 0.0, 0.0],
    temperature = [Ti],
    pusher = "boris",
    boundary_conditions = [["reflective", "reflective"]],
)
Species(
    name = "ion2",
    position_initialization = "regular",
    momentum_initialization = "mj",
    particles_per_cell = PPC, 
    c_part_max = 1.0,
    mass = MiMe*A,
    charge = Z,
    number_density = 0,
    mean_velocity = [0.0, 0.0, 0.0],
    temperature = [Ti],
    pusher = "boris",
    boundary_conditions = [["reflective", "reflective"]],
)
Species(
    name                    = "eon",
    position_initialization = "ion",
    momentum_initialization = "mj",
    particles_per_cell      = PPC, 
    c_part_max              = 1.0,
    mass                    = 1.0,
    charge                  = -1.0,
    number_density          = N0,
    mean_velocity           = [0.0, 0.0, 0.0],
    temperature             = [Te],
    pusher                  = "boris",
    boundary_conditions     = [["remove", "reflective"]],
)
ParticleInjector(
    name                    = "IonInjector",
    species                 = "ion2",
    box_side                = "xmax",
    position_initialization = "species",
    momentum_initialization = "mj",
    mean_velocity           = [-Vti/LightSpeed, 0.0, 0.0],
    temperature             = [Ti],
    number_density          = 1,
)
DiagScalar(
    every                   = 1,
    vars                    = [],
    precision               = 10
)
DiagFields(
    every                   = Dump,
    time_average            = 1,
    fields                  = [],
)
DiagParticleBinning(
    deposited_quantity      = "weight",
    every                   = Dump,
    time_average            = 1,
    species                 = ["ion"],
    axes                    = [ ["x", -0.01, 1.01, 25]],
)
DiagParticleBinning(
    deposited_quantity      = "weight",
    every                   = Dump,
    time_average            = 1,
    species                 = ["eon"],
    axes                    = [ ["x", -0.01, 1.01, 25]],
)
