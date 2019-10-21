using Syntrophy
using Plots
using DifferentialEquations
using LaTeXStrings
using Measures
import PyPlot

# This is a script to store the functions that plot the figures shown in my Syntrophy SI

# Function to calulate the value that needs to be exceeded for thermodynamic effects
# This is stored here as it is a function that has a mainly illustrative purpose
function Qineq(η::Float64,qm::Float64,m::Float64,kr::Float64,Keq::Float64)
    x = 0.1 # Starting with 10% as a rough guess
    Qi = x*Keq*(η*qm - m)/(η*qm + kr*m)
    return(Qi)
end

# Function to update population and nutrient concentrations
# This is run for a single population utilising a single reaction
function singlepop(du::Array{Float64,1},u::Array{Float64,1},p::Array{Float64,1},nuts::Array{Nut,1},reacs::Array{React,1},
                mic::Microbe,t::Float64)
    # Extract required parameters
    Y = p[1]
    # Extract relevant data from vector of nutrients
    α = nuts.↦:α
    δ = nuts.↦:δ
    con = nuts.↦:cst
    N = length(nuts) # Number of nutrients
    # And relevant data from vector of microbes
    η = mic.η
    m = mic.m # running for single microbe
    M = 1 # Number of microbes
    # Extract reaction stochiometry
    stc = (reacs.↦:stc)[1]
    ΔG0 = (reacs.↦:ΔG0)[1]
    # Now calculate q
    # p[2] = KS, p[3] = qm, p[4] = ΔGATP, p[5] = Temp, p[6] = kr
    q = qrate(u[1:N],p[2],p[3],p[4],ΔG0,p[5],stc,η,p[6])
    # Make vector to store nutrient changes due to consumption
    δX = zeros(N)
    for i = 1:length(stc)
        for j = N+1:N+M
            δX[i] += stc[i]*q*u[j] # Assumes single reaction
        end
    end
    # q has no dependance on population
    # Now update nutrients
    for i = 1:N
        if con[i] == false
            du[i] = α[i]-δ[i]*u[i]+δX[i]
        else
            du[i] = 0
        end
    end
    # Then calculate population changes
    for i = N+1:N+M
        j = i-N
        E = netE(η,q,m)
        if E >= 0.0 # find if growing or decaying
            du[i] = E*Y*u[i] # No dilution rate so can ignore
        else
            du[i] = E*Y*u[i]
        end
    end
    return(du)
end

# A function to make a plot of the maxcosump rate for a microbial population in glucose
function maxcosump()
    # Nutrient variables
    α = 5.55*10^(-6)
    δ = 2.00*10^(-4) #1.00*10^(-4)
    # make nutrients
    # 1 = glucose, 2 = oxegen, 3 = bicarbonate, 4 = hydrogen ion
    nuts = [Nut(1,false,α,δ),Nut(2,true,0,0),Nut(3,false,0,δ),Nut(4,true,0,0)]
    # Now make reactions
    ΔG0 = -2843800.0
    reac = [React(1,[1,2,3,4],[-1,-6,6,6],ΔG0)]
    # make set of microbes
    ηt = [collect(0:0.5:3.5);collect(4:4:32);collect(36:39);collect(39.25:0.25:43)]
    N = length(ηt) # (maximal η)-1
    r = 1 # Only reaction
    m = 2.16*10^(-19) # maintainance
    mics = Array{Microbe,1}(undef,N)
    for i = 1:N
        mics[i] = Microbe(ηt[i],m,r,0.0)
    end
    # Set intial population and nutrient concentrations for all cases
    pops = 100.0
    concs = zeros(length(nuts))
    # define initial concentrations
    concs[1] = 0.0555 # high initial concentration to ensure growth
    concs[2] = 0.21 # High value so oxegen isn't limiting
    concs[3] = 0.0 # No initial concentration
    concs[4] = 1.00*10.0^(-7) # pH 7
    # Define some constants
    ΔGATP = 75000.0 # Gibbs free energy of formation of ATP in a standard cell
    Temp = 312.0 # Temperature that growth is occuring at in Kelvin
    # All other constants require quite a bit of defining
    Y = 2.36*10.0^(13) # yield in cells per mole of ATP
    KS = 2.40*10.0^(-5) # saturation constant (substrate)
    qm = 3.42*10.0^(-18) # maximal rate substrate consumption mol cell s^-1
    kr = 1.0
    p = [Y,KS,qm,ΔGATP,Temp,kr]
    u0 = [concs;pops] # u0 and p same for every microbe
    tspan = (0.0,5000000.0)
    stoc = (reac.↦:stc)[1]

    # Preallocate vector and calculate maxium ATP generation rate
    atpgen = zeros(N)
    ηs = zeros(N)
    θs = zeros(N)
    fP = zeros(N)
    fS = zeros(N)
    fX = zeros(N)
    efP = zeros(N)
    efS = zeros(N)
    efX = zeros(N)
    Qi = zeros(N)
    Qa = zeros(N)

    for i = 1:N
        # Step to ensure steady state is reached for slower cases
        if ηt[i] <= 1.5
            tspan = (0.0,20000000.0)
        else
            tspan = (0.0,5000000.0)
        end
        settle = false
        # New version of function each time for each microbe
        f(du,u,p,t) = singlepop(du,u,p,nuts,reac,mics[i],t)
        # This step will take a while
        prob = ODEProblem(f,u0,tspan,p)
        sol = solve(prob,adaptive=false,dt=100) # turned dt down to make plots look nicer
        # Check if pop has changed more than 0.1% in last 10 time steps
        diff = (sol'[end,5]-sol'[end-10,5])/sol'[end,5]
        if diff < 0.001
            settle = true
        end
        # Print if it hasn't settled down to sufficently steady popultaion
        if settle == false
            println("Problem with eta = $(ηt[i])")
            println(diff)
        end
        # Now calulate max ATP genration rate
        ηc = mics[i].η # ATP usage amount
        qr = qrate(sol'[end,1:4],KS,qm,ΔGATP,ΔG0,Temp,stoc,ηc) # reaction rate
        atpgen[i] = sol'[end,5]*ηc*qr # multiply by population to get total rate
        θs[i] = θT(sol'[end,1:4],stoc,ΔGATP,ΔG0,ηc,Temp) # obtain thermodynamic term
        ηs[i] = ηc # save η for later plotting
        # Save steady state data
        fS[i] = sol'[end,1]
        fP[i] = sol'[end,3]
        fX[i] = sol'[end,5]
        # Save expected values as well
        if ηc == 0
            efS[i] = (α/δ)
            efP[i] = 0.0
            efX[i] = 0.0
            Qi[i] = NaN # Not reasonably defined in this case
        else
            efS[i] = m*KS/((ηc*qm - m)*concs[2]^6)
            if efS[i] > α/δ # Remove cases that require impossibly high substrate concentrations to be viable
                efS[i] = α/δ
            end
            efP[i] = 6*(α/δ - efS[i])
            efX[i] = (ηc/m)*(α - δ*efS[i])
            # Finally predict minimum value for Q
            KeQ = Keq(ΔG0,ηc,ΔGATP,Temp)
            Qi[i] = Qineq(ηc,qm,m,kr,KeQ)
            Qa[i] = QCoef(sol'[end,1:4],stoc)
        end
    end
    pyplot(dpi=200)
    width = 800
    height = 300
    # Plot maximal ATP generation rate
    p1 = plot(ηs,atpgen,xaxis=L"\eta\;\;mol_{ATP}\;(mol_{reaction})^{-1}",yaxis="Maximal ATP production rate mol/s",label="")
    p1 = plot!(p1,title="Production")
    # THIS DOES NOT WORK YET AS I HAVEN'T SET UP MODULE CORRECTLY YET
    px, py = annpos(ηs,atpgen)
    p1 = annotate!(p1,px,py,text("A",17,:black))
    # Plot thermodynamicinhibition term
    p2 = plot(ηs,θs,xaxis=L"\eta\;\;mol_{ATP}\;(mol_{reaction})^{-1}",yaxis=L"\theta\;\;",label="")
    p2 = plot!(p2,title="Inhibition")
    px, py = annpos(ηs,θs)
    p2 = annotate!(p2,px,py,text("B",17,:black))
    plot(p1,p2,layout=(1,2),size = (width,height),left_margin=5mm,right_margin=5mm,top_margin=10mm)
    savefig("Output/SynPlots/VarEff.png")
    # Quick plots of rougher variables
    plot(ηs,fS,title="S vs eta")
    plot!(ηs,efS)
    savefig("Output/SynPlots/Substrate.png")
    plot(ηs,fP,title="P vs eta")
    plot!(ηs,efP)
    savefig("Output/SynPlots/Product.png")
    plot(ηs,fX,title="X vs eta")
    plot!(ηs,efX)
    savefig("Output/SynPlots/Population.png")
    plot(ηs,log10.(Qi),label="Minimal Q",xaxis=L"\eta\;\;mol_{ATP}\;(mol_{reaction})^{-1}")
    plot!(ηs,log10.(Qa),label="Actual Q",yaxis=L"\log{(Q)}")
    savefig("Output/SynPlots/Quotient.png")
    return(nothing)
end

# function to plot a thermodynamically limited population and one without this limitation
function limunlim()
    # Nutrient variables
    α = 5.55*10^(-6)
    δ = 2.00*10^(-4)
    # make nutrients
    # 1 = glucose, 2 = oxegen, 3 = bicarbonate, 4 = hydrogen ion
    nuts = [Nut(1,false,α,δ),Nut(2,true,0,0),Nut(3,false,0,δ),Nut(4,true,0,0)]
    # Now make reactions
    ΔG0 = -2843800.0
    reac = [React(1,[1,2,3,4],[-1,-6,6,6],ΔG0)]
    # microbe variables
    ηs = [32.0,41.5]
    r = 1 # Only reaction
    m = 2.16*10^(-19) # maintainance
    # Considering 1 microbe with no maintaince and no dilution
    mics = [Microbe(ηs[1],m,r,0.0),Microbe(ηs[2],m,r,0.0)]
    # Set intial populations and nutrient concentrations
    pops = 100.0
    concs = zeros(length(nuts))
    # define initial concentrations
    concs[1] = 0.0555 # high initial concentration to ensure growth
    concs[2] = 0.21 # High value so oxegen isn't limiting
    concs[3] = 0.0 # No initial concentration
    concs[4] = 1.00*10.0^(-7) # pH 7
    # Define some constants
    ΔGATP = 75000.0 # Gibbs free energy of formation of ATP in a standard cell
    Temp = 312.0 # Temperature that growth is occuring at in Kelvin
    # All other constants require quite a bit of defining
    Y = 2.36*10.0^(13) # yield in cells per mole of ATP
    KS = 2.40*10.0^(-5) # saturation constant (substrate)
    qm = 3.42*10.0^(-18) # maximal rate substrate consumption mol cell s^-1
    kr = 1.0
    p = [Y,KS,qm,ΔGATP,Temp,kr]
    u0 = [concs;pops]
    tspan = (0.0,5000000.0)
    stoc = (reac.↦:stc)[1]

    # Change python variable name
    pyplot(dpi=200)
    f(du,u,p,t) = singlepop(du,u,p,nuts,reac,mics[1],t)
    prob = ODEProblem(f,u0,tspan,p)
    sol = solve(prob,adaptive=false,dt=500) # turned dt down to make plots look nicer
    # Plot nutrient concentrations on same graph
    Lη = L"\eta"
    p1 = plot(sol.t,[sol'[:,1],sol'[:,3]],ylabel="Concentration M L^-1",label=["substrate" "product"],title="No inhibition ($(Lη) = $(ηs[1]))")
    # Find ideal annotation position and annotate
    px, py = annpos(sol.t,[sol'[:,1];sol'[:,3]])
    p1 = annotate!(p1,px,py,text("A",17,:black))
    # Then plot population on another subplot
    p3 = plot(sol.t,sol'[:,5],label="",ylabel="Cell density L^-1")
    px, py = annpos(sol.t,sol'[:,5])
    p3 = annotate!(p3,px,py,text("C",17,:black))

    # Then do same for limited case
    g(du,u,p,t) = singlepop(du,u,p,nuts,reac,mics[2],t)
    prob = ODEProblem(g,u0,tspan,p)
    sol = solve(prob,adaptive=false,dt=500) # turned dt down to make plots look nicer
    # Plot nutrient concentrations on same graph
    p2 = plot(sol.t,[sol'[:,1],sol'[:,3]],ylabel="Concentration M L^-1",label=["substrate" "product"],title="Inhibition ($(Lη) = $(ηs[2]))")
    px, py = annpos(sol.t,[sol'[:,1];sol'[:,3]])
    p2 = annotate!(p2,px,py,text("B",17,:black))
    # Then plot population on another subplot
    p4 = plot(sol.t,sol'[:,5],label="",ylabel="Cell density L^-1")
    px, py = annpos(sol.t,sol'[:,5])
    p4 = annotate!(p4,px,py,text("D",17,:black))
    # Set overall plot height and width and make combined plot
    width = 800
    height = 600
    plot(p1,p2,p3,p4,layout=(2,2),size = (width, height),xlabel="time s",left_margin=5mm,right_margin=5mm,top_margin=10mm)
    savefig("Output/SynPlots/LimUnlim.png")
    return(nothing)
end


# A function that calulcates maximal ATP production rates for a given supply rate and decay rate
# It then outputs this data to be plotted in another script
function varcosump()
    # Nutrient variables
    α = 5.55*10^(-6) # middle
    δ = 2.00*10^(-4) # middle
    # make nutrients
    # 1 = glucose, 2 = oxegen, 3 = bicarbonate, 4 = hydrogen ion
    nuts = [Nut(1,false,α,δ),Nut(2,true,0,0),Nut(3,false,0,δ),Nut(4,true,0,0)]
    # Now make reactions
    ΔG0 = -2843800.0
    reac = [React(1,[1,2,3,4],[-1,-6,6,6],ΔG0)]
    # make set of microbes
    ηt = [collect(0:0.5:3.5);collect(4:4:32);collect(36:39);collect(39.25:0.25:43)]
    N = length(ηt) # (maximal η)-1
    r = 1 # Only reaction
    m = 2.16*10^(-19) # maintainance
    mics = Array{Microbe,1}(undef,N)
    for i = 1:N
        mics[i] = Microbe(ηt[i],m,r,0.0)
    end
    # Set intial population and nutrient concentrations for all cases
    pops = 100.0
    concs = zeros(length(nuts))
    # define initial concentrations
    concs[1] = 0.0555 # high initial concentration to ensure growth
    concs[2] = 0.21 # High value so oxegen isn't limiting
    concs[3] = 0.0 # No initial concentration
    concs[4] = 1.00*10.0^(-7) # pH 7
    # Define some constants
    ΔGATP = 75000.0 # Gibbs free energy of formation of ATP in a standard cell
    Temp = 312.0 # Temperature that growth is occuring at in Kelvin
    # All other constants require quite a bit of defining
    Y = 2.36*10.0^(13) # yield in cells per mole of ATP
    KS = 2.40*10.0^(-5) # saturation constant (substrate)
    qm = 3.42*10.0^(-18) # maximal rate substrate consumption mol cell s^-1
    p = [Y,KS,qm,ΔGATP,Temp]
    u0 = [concs;pops] # u0 and p same for every microbe
    tspan = (0.0,5000000.0)
    stoc = (reac.↦:stc)[1]

    # Preallocate vector and calculate maxium ATP generation rate
    atpgen = zeros(N)
    ηs = zeros(N)
    θs = zeros(N)
    # Need some method of checking if that results are under control
    for i = 1:N
        settle = false
        # New version of function each time for each microbe
        f(du,u,p,t) = singlepop(du,u,p,nuts,reac,mics[i],t)
        # This step will take a while
        prob = ODEProblem(f,u0,tspan,p)
        sol = solve(prob,adaptive=false,dt=100) # turned dt down to make plots look nicer
        # Check if pop has changed more than 0.1% in last 10 time steps
        diff = (sol'[end,5]-sol'[end-10,5])/sol'[end,5]
        if diff < 0.001
            settle = true
        end
        # Print if it hasn't settled down to sufficently steady popultaion
        if settle == false
            println("Problem with eta = $(ηt[i])")
            println(diff)
        end
        # Now calulate max ATP genration rate
        ηc = mics[i].η # ATP usage amount
        qr = qrate(sol'[end,1:4],KS,qm,ΔGATP,ΔG0,Temp,stoc,ηc) # reaction rate
        atpgen[i] = sol'[end,5]*ηc*qr # multiply by population to get total rate
        θs[i] = θT(sol'[end,1:4],stoc,ΔGATP,ΔG0,ηc,Temp) # obtain thermodynamic term
        ηs[i] = ηc # save η for later plotting
    end
    # Define structure to output
    out = Array{Float64,2}(undef,N+1,3)
    # First line essentially a header
    out[1,1] = NaN
    out[1,2] = α
    out[1,3] = δ
    for i = 1:N
        out[i+1,1] = ηs[i]
        out[i+1,2] = θs[i]
        out[i+1,3] = atpgen[i]
    end
    # Now write out output data to file
    output_file = "Data/delta$(δ)alpha$(α).csv"
    out_file = open(output_file, "w")
    for i = 1:size(out,1)
        line = ""
        for j = 1:size(out,2)
            line *= "$(out[i,j]),"
        end
        # remove surplus ,
        line = line[1:end-1]
        # add a new line
        line *= "\n"
        write(out_file,line)
    end
    return(nothing)
end

# Function to read in sort and plot the maximum atp production rate data
function plotvar()
    # Get all file names
    files = readdir("Data")
    # remove non-csvs
    files = filter!(x->x[end-3:end]==".csv",files)
    # Calculate number of valid files
    L = length(files)
    # setup plots
    pyplot(dpi=200)
    p1 = plot(title="Variation of removal rate",xaxis=L"\eta\;\;mol_{ATP}\;(mol_{reaction})^{-1}",ylabel="Maximal ATP production rate mol/s")
    p2 = plot(title="Variation of supply rate",xaxis=L"\eta\;\;mol_{ATP}\;(mol_{reaction})^{-1}",ylabel="Maximal ATP production rate mol/s")
    p3 = plot(xaxis=L"\eta\;\;mol_{ATP}\;(mol_{reaction})^{-1}",ylabel=L"\theta\;\;")
    p4 = plot(xaxis=L"\eta\;\;mol_{ATP}\;(mol_{reaction})^{-1}",ylabel=L"\theta\;\;")
    # Find length of file
    testfile = "Data/$(files[1])"
    l = countlines(testfile)-1
    # Make data structure to store all data
    datay = zeros(4,3,l)
    datax = zeros(l)
    # And some counters
    i1 = i2 = i3 = i4 = 1
    # Read in files one by one to plot
    for i = 1:L
        infile = "Data/$(files[i])"
        # now read in data
        l = countlines(infile)
        w = 3
        data = zeros(l-1,w)
        # Create storage for supply and removal rates
        α = 0.0
        δ = 0.0
        open(infile, "r") do in_file
            # Use a for loop to process the rows in the input file one-by-one
            k = 0 # first line is header and must be treated seperately
            for line in eachline(in_file)
                # parse line by finding commas
                len = length(line)
                comma = fill(0,w+1)
                j = 1
                for i = 1:len
                    if line[i] == ','
                        j += 1
                        comma[j] = i
                    end
                end
                comma[end] = len+1
                # first line
                if k == 0
                    α = parse(Float64,line[(comma[2]+1):(comma[3]-1)])
                    δ = parse(Float64,line[(comma[3]+1):(comma[4]-1)])
                    # Round both to 3 sf for nicer label names
                    α = round(α; sigdigits=3)
                    δ = round(δ; sigdigits=3)
                else
                    for i = 1:w
                        data[k,i] = parse(Float64,line[(comma[i]+1):(comma[i+1]-1)])
                    end
                end
                k += 1
            end
        end
        # Use α and δ to make nice labels
        La = L"\alpha\;=\;"
        Ld = L"\delta\;=\;"
        lb = "$(La)$(α) $(Ld)$(δ)"
        # Plot the data
        # Filter so that the graphs are better comparisons
        if α ≈ 5.55*10^(-6) && δ ≈ 2.00*10^(-4)
            # Do plotting and then store appropiate data
            p1 = plot!(p1,data[:,1],data[:,3],label=lb,color=i)
            p2 = plot!(p2,data[:,1],data[:,3],label=lb,color=i)
            p3 = plot!(p3,data[:,1],data[:,2],label=lb,color=i)
            p4 = plot!(p4,data[:,1],data[:,2],label=lb,color=i)
            # Store appropiate data
            datay[1,i1,:] = data[:,3]
            datay[2,i2,:] = data[:,3]
            datay[3,i3,:] = data[:,2]
            datay[4,i4,:] = data[:,2]
            # increment appropiate counters
            i1 += 1
            i2 += 1
            i3 += 1
            i4 += 1
        elseif α ≈ 5.55*10^(-6)
            p1 = plot!(p1,data[:,1],data[:,3],label=lb,color=i)
            p3 = plot!(p3,data[:,1],data[:,2],label=lb,color=i)
            datay[1,i1,:] = data[:,3]
            datay[3,i3,:] = data[:,2]
            # increment appropiate counters
            i1 += 1
            i3 += 1
        elseif δ ≈ 2.00*10^(-4)
            p2 = plot!(p2,data[:,1],data[:,3],label=lb,color=i)
            p4 = plot!(p4,data[:,1],data[:,2],label=lb,color=i)
            datay[2,i2,:] = data[:,3]
            datay[4,i4,:] = data[:,2]
            # increment appropiate counters
            i2 += 1
            i4 += 1
        end
        if i == L
            # Store the unchanging x data
            datax = data[:,1]
        end
    end
    # Annotation step
    px, py = annpos(datax,[datay[1,1,:];datay[1,2,:];datay[1,3,:]])
    p1 = annotate!(p1,px,py,text("A",17,:black))
    px, py = annpos(datax,[datay[2,1,:];datay[2,2,:];datay[2,3,:]])
    p2 = annotate!(p2,px,py,text("B",17,:black))
    px, py = annpos(datax,[datay[3,1,:];datay[3,2,:];datay[3,3,:]])
    p3 = annotate!(p3,px,py,text("C",17,:black))
    px, py = annpos(datax,[datay[4,1,:];datay[4,2,:];datay[4,3,:]])
    p4 = annotate!(p4,px,py,text("D",17,:black))
    # set total size of plot and then plot all together
    width = 800
    height = 600
    plot(p1,p2,p3,p4,layout=(2,2),xlabel=L"\eta\;\;mol_{ATP}\;(mol_{reaction})^{-1}",size=(width,height),left_margin=5mm,right_margin=5mm,top_margin=10mm)
    plot!(ann = (:top_left, :auto)) # give figures letters to identify
    savefig("Output/SynPlots/ThermComp.png")
    return(nothing)
end


# This whole script can stand a lot of improvement
# @time maxcosump()
@time limunlim()
# @time plotvar()
