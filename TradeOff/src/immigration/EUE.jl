using TradeOff

"""
    Function to calculate the amount of free energy dissipated
    calculate_D(
    ΔG0 :: Float64, Standard Gibb's free energy change when 1 mole of reaction occurs
    R :: Float64, Gas constant
    T :: Float64, Temperature
    S::Float64, Substrate concentration
    P::Float64, Product concentration
    η::Float64, Amount of free energy obtained from given reaction
    ΔGATP :: Float64)

TBW
"""
function calculate_D(
    ΔG0 :: Float64,
    R :: Float64,
    T :: Float64,
    S::Float64,
    P::Float64,
    η::Float64,
    ΔGATP :: Float64)

    D = ΔG0 + (R * T * log1p(Q(S, P))) + (η * ΔGATP)
    return (D)
end 

"""
    Function to calculate the change in free energy of the reaction.
    calculate_ΔGT(
    ΔG0 :: Float64,
    R :: Float64,
    T :: Float64,
    S::Float64,
    P::Float64
    )

TBW
"""
function calculate_ΔGT(
    ΔG0 :: Float64,
    R :: Float64,
    T :: Float64,
    S::Float64,
    P::Float64
    )
    ΔGT = ΔG0 + (R * T * log1p(Q(S, P)))

    return (ΔGT)
end 

"""
    Function that calculates the thermodynamic factor of a reaction
    calculate_θ(
    S::Float64, Substrate concentration
    P::Float64, Product concentration
    T::Float64, Temperature
    η::Float64, Amount of free energy obtained from given reaction (units of number of ATP molecules)
    ΔG0::Float64 Gibbs free-enrgy per mole of ATP
    )
    Keq is a function to find the equilibrium constant defined in the TradeOff package
    Q is a function to find the reaction quotient Q, in the case of 1 to 1 stoichiometry also defined in the TradeOff package
"""
function calculate_θ(
    S::Float64,
    P::Float64,
    T::Float64,
    η::Float64,
    ΔG0::Float64
    )
    # Catch perverse cases that sometimes arise
    if S <= 0.0
        θs = 1.0
    elseif P <= 0.0
        θs = 0.0
    else
        θs = Q(S, P) / Keq(T, η, ΔG0)
    end
    # θ can be greater than 1, this does not have any impact as q cannot be negative
    return (θs)
end

"""
    Function to find species i's reaction rate for reaction α
    calculate_q(
    S::Float64, Substrate concentration
    P::Float64, Product concentration
    E::Float64, Enzyme copy number
    i::Int64, 
    ps::Microbe, 
    T::Float64, Temperature
    r::Reaction
    )

TBW
"""
function calculate_q(S::Float64,
    P::Float64,
    E::Float64,
    i::Int64,
    ps::Microbe,
    T::Float64,
    r::Reaction)
    # To speed things I don't have a check here to ensure that r.ID matches ps.Reac[i]
    # This is something to check if I start getting errors
    θ = calculate_θ(S, P, T, ps.η[i], r.ΔG0)
    q = ps.kc[i] * E * S * (1 - θs) / (ps.KS[i] + S * (1 + ps.kr[i] * θ))
    # Ensure that negative value cannot be returned
    return (max(q, 0.0))
end

"""
    Function that calculates the Energy Use Efficiency (EUE) for a particular species
    calculating_EUE(
    D :: Float64
    ΔGT::Float64
    q :: Float64
    )

TBW
"""
function calculating_EUE(
    D :: Float64
    ΔGT::Float64
    q :: Float64
    )
    # Initialise an array for vectors of EUEs for each species
    EUEs_over_time = []
    for timepoint in 1:T
        #empty vectors to store EUE values
        species_EUE_list =[]
        #loop over all species
        for species in 1:list_of_species
            numerator_list = []
            denominator_list = []
            for reaction in 1:reactions
                numerator = (1-calculate_D()/calculate_ΔGT()) * calculate_q()
                denominator = calculate_q()
                push!(numerator_list, numerator)
                push!(denominator_list, denominator)
            end
            species_EUE = sum(numerator_list)/sum(denominator_list)
            push!(species_EUE_list, species_EUE)
        end
        push!(EUEs_over_time, species_EUE_list)
    end

    return EUE
end
