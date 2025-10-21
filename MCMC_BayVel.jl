# ========================
# Activate Project Environment
# ========================
using Pkg
pathToYourDirectory::String = "pathToYourDirectory"
BayVelPath::String = pathToYourDirectory * "/BayVel"
Pkg.activate(pathToYourDirectory) # Activates the Julia environment in the specified directory, ensuring reproducibility.


# ========================
# Install Required Julia Packages (DECOMMENT JUST THE FIRST TIME)
# ========================
# Pkg.add("LanguageServer")
# Pkg.add("Revise")
# Pkg.add("Distributions")
# Pkg.add("Random")
# Pkg.add("LinearAlgebra")
# Pkg.add("PDMats")
# Pkg.add("StatsBase")
# Pkg.add("RCall")
# Pkg.add("ToggleableAsserts")
# Pkg.add("ProgressMeter")
# Pkg.add("LogExpFunctions")
# Pkg.develop(path = BayVelPath)

# ========================
# Load Required Julia Packages
# ========================
using LanguageServer           # Support for language tools like IDE integration
using Revise                   # Automatically track code changes without restarting session
using Distributions, Random    # Tools for probability distributions and random number generation
using LinearAlgebra, PDMats, StatsBase  # Basic matrix and statistical operations
using RCall                    # Interface for calling R from Julia
using ToggleableAsserts        # Asserts that can be turned off (e.g., for performance)
using ProgressMeter            # For showing progress bars
using LogExpFunctions          # Numerically stable exponential/log functions
using BayVel                   # Custom package of BayVel



# ========================
# Define Simulation Configuration
# ========================
typeSIM::String = "sim"                 # Dataset or tissue type (possible values: "sim", "DentateGyrus", "Pancreas")
typeSW::String = "SW1"                  # Switching structure type (typeSW = "SW1", "SW10")
typeT::String = "T1"                    # Subgroup structure (typeT = "T1", "T2", "T3" for simulated data, typeT = "T1", "T2" for real data)
typeD::String = "D1"                    # Likelihood distribution ("D1" for Poisson, "D4" for Negative Binomial)

model::BayVel.groupSubgroup = BayVel.groupSubgroup()
n_genes::Int64 = 2000                   # number of simulated genes

# Pass variables to R using RCall
@rput typeSIM
@rput typeSW
@rput typeT
@rput typeD
@rput n_genes
@rput pathToYourDirectory

# ========================
# Load Data from R Depending on Simulation Type
# ========================
R"""
real <- FALSE  # Flag for real data vs simulated

nameSim <- paste0(typeSW, "-", typeT, "-", typeD)
# Load DentateGyrus or Pancreas dataset
if(startsWith(typeSIM, "DentateGyrus") | startsWith(typeSIM, "Pancreas")){
  path <- paste0(pathToYourDirectory, "/real data/", typeSIM, "/filter/")
  real <- TRUE
  n_clusters <- ifelse(typeT == "T1", 1, "mult")
  nameSimREAL <- paste0(typeSIM, n_genes, "genes_", n_clusters, "subgr.RData")
  load(paste0(path, nameSimREAL)) # load real data
}else{
  path <- paste0(pathToYourDirectory, "/simulations/")
  load(paste0(path, nameSim, "/", nameSim, "_", n_genes, ".RData")) # Load simulated data
}

# Set switching time labels depending on SW type
if(typeSW == "SW1"){
  typeCellT0_off <- rep(1, length(typeCell))
}else{
  typeCellT0_off <- typeCell
}

# Format simulated data matrices
if(!real){
  attr(u0_off_real , "dim") <- c(n_typeT0_off, n_genes)
  attr(s0_off_real , "dim") <- c(n_typeT0_off, n_genes)
}

# Number of cells
n_cells = dim(Y_u)[1]

# Print some diagnostic info
print("# ----------------------------")
print(paste0("Type of simulation: ", typeSIM))
print(paste0("Name of the simulation: ", nameSim))
print(paste0("Number of cells: ", n_cells))
print(paste0("Number of genes: ", n_genes))
print(paste0("Number of subgroups: ", length(unique(subtypeCell))))
print(paste0("Number of groups: ", length(unique(typeCell))))
print(paste0("Number of switching clusters: ", length(unique(typeCellT0_off))))
print("# ----------------------------")
"""

# ========================
# Retrieve Data from R into Julia
# ========================
@rget real        # Is the dataset real or simulated
@rget typeSIM     # Simulation type
@rget nameSim     # Simulation identifier
@rget subtypeCell # Subclass labels for each cell
@rget typeCell    # Class labels for each cell
@rget typeCellT0_off # Switching labels
@rget Y_u         # Unspliced gene expression
@rget Y_s         # Spliced gene expression
@rput n_genes     # Number of genes


################################################
#### Parameter Simulation Setup
################################################

# Flags to indicate which parameters should be simulated (used for code debugging):
simSS::Bool   = true   # Simulate steady-state parameters
simT0::Bool   = true   # Simulate switching time parameters
simTau::Bool  = true   # Simulate subgroups time parameters
simEta::Bool  = true   # Simulate overdispersion parameters 
simCatt::Bool = true   # Simulate capture efficiency effects

# Based on dynamics type, selectively disable simulation of certain parameters
if (typeD == "D1") || (typeD == "D3")
    simCatt = false
end
if (typeD == "D1") || (typeD == "D2")
    simEta = false
end

# Send simulation flags to R
@rput simSS simT0 simTau simEta simCatt

################################################
#### Random Seed Setup
################################################

rngseed::Int64 = 1234  # Set the seed for reproducibility
@rput rngseed
Random.seed!(rngseed)


################################################
#### MCMC Configuration
################################################
molt::Int64 = 5  # Multiplier for MCMC steps
mcmc = (
    iter = 50000 * molt,     # Total number of iterations
    thin = 5 * molt,         # Thinning factor
    burnin = 40000 * molt    # Burn-in period
)

# Extract individual MCMC settings
mcmcIter     = deepcopy(mcmc.iter)
mcmcBurnin::Int64 = deepcopy(mcmc.burnin)
mcmcThin::Int64   = deepcopy(mcmc.thin)

# Compute how many samples to retain
SampleToSave::Int64 = Int64(trunc((mcmc.iter - mcmc.burnin) / mcmc.thin))

# Send MCMC configuration to R
@rput mcmcIter mcmcBurnin mcmcThin

# Adaptation parameters (used for adaptation of the variance parameter)
adaptType::String = "A/(B+n)"                   # Adaptation strategy
stepVect::Vector{Float64} = [2500.0, 20000.0]   # Values used for adaptation


################################################
#### Data Structure 
################################################
# Number of initial switching point cell types
n_typeT0_off::Int64 = maximum(typeCellT0_off)

# Number of group and subgroups
n_typeC::Int64 = maximum(typeCell)
n_subtypeC::Int64 = maximum(subtypeCell)

# Number of genes and cells in the dataset
n_genes::Int64 = size(Y_u, 2)
n_cells::Int64 = size(Y_u, 1)

################################################
#### Load (for simulated data) or Initialize (for real data) True Parameters
################################################
if(!real)
    @rget alpha_real;       # Transcription rate
    @rget beta_real;        # Splicing rate
    @rget gamma_real;       # Degradation rate
    @rget t_real;           # Subgroup time (per subgroup and gene)
    @rget tau_real;         # t - t0_on
    @rget t0_off_real;      # Off switching time (OFF-state), per group type and gene
    @rget t0_on_real;       # On switching time (ON-state), per gene, usually fixed at zero (identifiability)
    @rget u0_off_real;      # Unspliced mRNA at t0_off
    @rget s0_off_real;      # Spliced mRNA at t0_off
    @rget u0_on_real;       # Unspliced mRNA at t0_on 
    @rget s0_on_real;       # Spliced mRNA at t0_on 
    @rget k_real;           # State indicator (0 = off, 1 = on)
    @rget eta_real;         # Gene-specific overdispersion 
    @rget catt_real;        # Cell-specific capture efficiency

    k_real[:, :] .= k_real ./ 2 # Normalize the switching matrix if needed
else 
    # If using real data, allocate empty placeholders for parameters
    alpha_real = Matrix{Float64}(undef, n_genes, 2)
    beta_real  = ones(Float64, n_genes)
    gamma_real = Vector{Float64}(undef, n_genes)
    t_real     = Matrix{Float64}(undef, n_cells, n_genes)
    tau_real   = Matrix{Float64}(undef, n_cells, n_genes)
    t0_off_real = Matrix{Float64}(undef, n_typeT0_off, n_genes)
    t0_on_real  = Vector{Float64}(undef, n_genes)
    u0_off_real = Matrix{Float64}(undef, n_typeT0_off, n_genes)
    s0_off_real = Matrix{Float64}(undef, n_typeT0_off, n_genes)
    u0_on_real  = Vector{Float64}(undef, n_genes)
    s0_on_real  = Vector{Float64}(undef, n_genes)
    k_real      = Matrix{Int64}(undef, 2, n_genes)
    eta_real    = Vector{Float64}(undef, n_genes)
    catt_real   = Vector{Float64}(undef, n_cells)

    # Set defaults for eta and catt depending on dynamics type
    if (typeD == "D1") || (typeD == "D3")
        catt_real[:] .= 1.0
    end
    if (typeD == "D1") || (typeD == "D2")
        eta_real[:] .= 0.0
    end
end

################################################
#### Initialize Parameters for Inference
################################################
# Vectors for steady state initial values
initUoff   = Vector{Float64}(undef, n_genes)  # For unspliced RNA steady-state OFF
initSoff   = Vector{Float64}(undef, n_genes)  # For spliced RNA steady-state OFF
initDiffU  = Vector{Float64}(undef, n_genes)  # Initial differences for unspliced steady state coordinates
initBeta   = Vector{Float64}(undef, n_genes)  # For splicing rate beta

# Matrices for initial conditions at t0_off (per group)
initT0_off  = Matrix{Float64}(undef, n_typeT0_off, n_genes)
initU0_off  = Matrix{Float64}(undef, n_typeT0_off, n_genes)
initS0_off  = Matrix{Float64}(undef, n_typeT0_off, n_genes)

# Optional placeholders for initial conditions at t0_on
initU0_on = Matrix{Float64}(undef, 1, n_genes)
initS0_on = Matrix{Float64}(undef, 1, n_genes)

# Cell-by-gene latent time initialization
initTStar_withM = Matrix{Float64}(undef, size(Y_u))
initTStar       = Matrix{Float64}(undef, size(Y_u))
initTau         = Matrix{Float64}(undef, size(Y_u))

# Initialization of cell-by-gene branch indicators, overdispersion and capture efficiency
initK     = Matrix{Int64}(undef, size(Y_u))
initEta   = Vector{Float64}(undef, n_genes)
initCatt  = Vector{Float64}(undef, size(Y_u, 1))

# Initialize parameters based on truth + simulation flags
initParam!(
    initUoff, initSoff, initDiffU, initBeta,
    initT0_off, initU0_off, initS0_off, initU0_on, initS0_on,
    initTStar_withM, initTStar, initTau, initK, initEta, initCatt,
    alpha_real, beta_real, gamma_real, t0_off_real, u0_off_real, s0_off_real,
    t_real, tau_real, Int64.(k_real), eta_real, catt_real,
    simSS, simT0, simTau, simEta, simCatt,
    Int64.(subtypeCell), Int64.(typeCell), Int64.(typeCellT0_off),
    model
)

################################################
#### Sanity Checks for Initialization
################################################

# Check that k is binary (0 or 1)
@toggled_assert all((initK .== 0) .| (initK .== 1)) "Invalid initialization: K must be binary (0 or 1)"

# If steady-state is not simulated, check that initializations match real values
@toggled_assert simSS | (!simSS & all(initUoff[:] .== (alpha_real[:,1] ./ beta_real[:]))) "Mismatch in U_off initialization"
@toggled_assert simSS | (!simSS & all(initSoff[:] .== (alpha_real[:,1] ./ gamma_real[:]))) "Mismatch in S_off initialization"
@toggled_assert simSS | (!simSS & all(initDiffU[:] .== ((alpha_real[:,2] .- alpha_real[:,1]) ./ beta_real[:]))) "Mismatch in DiffU initialization"
@toggled_assert simSS | (!simSS & all(initBeta[:] .== beta_real[:])) "Mismatch in Beta initialization"

# If t0 is not simulated, check that it matches real values
@toggled_assert simT0 | (!simT0 & all(initT0_off .== t0_off_real)) "Mismatch in t0_off initialization"

# If neither SS nor T0 are simulated, check initial conditions at t0
@toggled_assert (simT0 | simSS) | (!simT0 & !simSS & all(round.(initU0_off; digits=10) .== round.(u0_off_real; digits=10))) "Mismatch in u0_off initialization"
@toggled_assert (simT0 | simSS) | (!simT0 & !simSS & all(round.(initS0_off; digits=10) .== round.(s0_off_real; digits=10))) "Mismatch in s0_off initialization"

# If Tau is not simulated, check if time is correctly initialized
@toggled_assert simTau | (!simTau & all(initTStar .== t_real)) "Mismatch in TStar initialization"

# Ensure that for cells with switching state = 1, tau ≤ t0_off
for g in 1:size(initK, 2)
    for sty in 1:n_subtypeC
        cellSubTy = findall(Int64.(subtypeCell) .== sty)
        tyT0_off = Int64(typeCellT0_off[cellSubTy[1]])
        if initK[cellSubTy[1], g] == 1
            @toggled_assert initTau[cellSubTy[1], g] <= initT0_off[tyT0_off, g] "Tau > t0_off for gene $g"
        end
    end
end

# If eta or catt are not simulated, check for match with real values
@toggled_assert simEta  | (!simEta  & all(initEta[:] .== eta_real[:])) "Mismatch in Eta initialization"
@toggled_assert simCatt | (!simCatt & all(initCatt[:] .== catt_real[:])) "Mismatch in Catt initialization"

# Disable all toggleable assertions from this point onward (to improve performance)
toggle(false)




################################################
#### Prior Distributions for Steady-State Parameters
################################################
# Maximum unscaled RNA count (used for scaling the priors of the steady states)
a::Float64 = 3000  # Note: could be computed with quantiles of the data

# Hyperparameters for priors
par1::Float64 = 1.0
par2::Float64 = 1.0
par3::Float64 = 1.0

# Set priors for steady state
priorsLogU_SS_on   = BayVel.Log_Beta(par1 + par2, par3, a)
priorsLogS_SS_on   = BayVel.Log_Beta(par1 + par2, par3, a)
priorsLogU_SS_off  = BayVel.Log_Beta(par1, par2, 1.0)
priorsLogBeta      = BayVel.Log_Beta(par1, par2, 1.0)

################################################
#### Priors for Time Parameters
################################################
# Define prior distribution over subgroups time (t) with mass at zero (t = 0)
p_zero::Float64 = 0.01
priorsT = BayVel.UniformWithMass(p_zero)

# Define prior over log of t0_off using exponential distribution
priorsLogT0_off = BayVel.Log_Exponential(1 / initBeta[1]) 

# Ensure simulation flags are re-synced to R (for later MCMC usage)
simSS = true
simT0 = true
simTau = true
@rput simSS simT0 simTau

################################################
#### Run MCMC Inference
################################################

# Call the core MCMC sampling function, passing:
# - the model
# - observed data (unspliced/spliced matrices)
# - priors
# - initialization
# - control flags for simulation
LogSS_chain, SS_Star_chain, LogT0_off_chain, T0_off_chain, phi_chain, TStar_chain, TStar_withM_chain, Tau_chain, k_chain, LogEta_chain, LogitCatt_chain, acceptRateSS, acceptRateT0_off, acceptRateTStar, acceptRateEta, acceptRateCatt = MCMC(
    model,
    Int64.(Y_u), # unspliced
    Int64.(Y_s), # spliced
    Int64.(typeCell); # typeCell
    mcmc = mcmc,
    priorsLogU_SS_on = priorsLogU_SS_on, 
    priorsLogS_SS_on = priorsLogS_SS_on,
    priorsLogU_SS_off = priorsLogU_SS_off, 
    priorsLogBeta = priorsLogBeta, 
    priorsLogT0_off = priorsLogT0_off,
    priorsT = priorsT,
    priorsEta = Truncated(Normal(0.0, 10000.0), 0.0, Inf), 
    priorsCatt = Uniform(0.0,1.0),
    initUoff = initUoff,
    initSoff = initSoff,
    initDiffU = initDiffU,
    initBeta = initBeta,
    initT0_off = initT0_off,
    initTStar = initTStar,    
    initTau = initTau,
    initEta = initEta,
    initCatt = initCatt,
    simSS = simSS,
    simT0 = simT0,
    simTau = simTau, 
    simEta = simEta,
    simCatt = simCatt,
    alphaTarget = 0.25,
    stepVect = stepVect,
    typeCellT0_off = Int64.(typeCellT0_off),
    subtypeCell = Int64.(subtypeCell)
);



################################################
#### Send Results to R
################################################
# Transfer all posterior chains and metadata to R
@rput LogSS_chain SS_Star_chain LogT0_off_chain T0_off_chain phi_chain
@rput TStar_chain TStar_withM_chain Tau_chain k_chain
@rput LogEta_chain LogitCatt_chain
@rput mcmcIter mcmcBurnin mcmcThin
@rput acceptRateSS acceptRateT0_off acceptRateTStar acceptRateEta acceptRateCatt
@rput adaptType stepVect

# Record prior distribution types as strings (for reproducibility/debugging)
typePriorT0_off::String = string(typeof(priorsLogT0_off))
typePriorT::String = string(typeof(priorsT))
@rput typePriorT0_off typePriorT


R"""
if(typeSIM == "sim"){
    dir.create(paste0(pathToYourDirectory,"/simulations/", nameSim), showWarnings = FALSE)
    dir.create(paste0(pathToYourDirectory,"/simulations/", nameSim, "/output/"), showWarnings = FALSE)
    save.image(file = paste0(pathToYourDirectory, "/simulations/", nameSim, "/output/res_", typeSIM, "_", nameSim, "_", n_genes, "genes_", mcmcIter, ".RData"))
}else{
    dir.create(paste0(pathToYourDirectory, "/real data/", typeSIM,  "/filter/", nameSim), showWarnings = FALSE)
    dir.create(paste0(pathToYourDirectory, "/real data/", typeSIM,  "/filter/", nameSim, "/output"), showWarnings = FALSE)

    save.image(file = paste0(pathToYourDirectory, "/real data/", typeSIM, "/filter/", nameSim, "/output/res_", typeSIM, "_", nameSim, "_", n_genes, "genes_", mcmcIter, ".RData"))
}
"""
