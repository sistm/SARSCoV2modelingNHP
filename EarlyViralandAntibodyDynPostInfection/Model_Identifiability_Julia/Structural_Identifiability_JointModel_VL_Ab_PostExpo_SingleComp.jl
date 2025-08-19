using StructuralIdentifiability
using FileIO


#= Definition of the ODE model using the @ODEmodel macro 
--------------------------------------------------------
Presentation of the model - Single URT compartment:
x1 = T:     Uninfected target cells 
x2 = I1:    Unproductively infected target cells 
x3 = I2:    Productively infected target cells 
x4 = VI:    Infectious viruses 
x5 = VNI:   Non-infectious viruses 
x6 = Vini:  Inoculated viruses 

x7 = S:     SL ASCs
x8 = BAb:   Binding antibodies
x9 = NAb :  Neutralizing antibodies

y1 = gRNA:    genomic viral load
y2 = sgRNA:   subgenomic viral load
y3 = Obs_BAb:   binding antibodies
y4 = Obs_NAb:   neutralizing antibodies
=# 

ODE_model_Single_URT_Comp = @ODEmodel(
    x1'(t) = -beta*(1 - eta*x9(t)/(1+eta*x9(t)))*x1(t)*x4(t) - 1e-3*beta*(1 - eta*x9(t)/(1+eta*x9(t)))*x1(t)*x6(t),
    x2'(t) = beta*(1 - eta*x9(t)/(1+eta*x9(t)))*x1(t)*x4(t) + 1e-3*beta*(1 - eta*x9(t)/(1+eta*x9(t)))*x1(t)*x6(t) - 3*x2(t),
    x3'(t) = 3*x2(t) - delta*x3(t),
    x4'(t) = P*1e-3*x3(t) - 3*x4(t) - beta*(1 - eta*x9(t)/(1+eta*x9(t)))*x1(t)*x4(t),
    x5'(t) = P*(1-1e-3)*x3(t) - 3*x5(t),
    x6'(t) = -20*x6(t) - 1e-3*beta*(1 - eta*x9(t)/(1+eta*x9(t)))*x1(t)*x6(t),

    x7'(t) = 1/100000*(x4(t) + x5(t)) ^ 1/4 + rho*(1-x7(t)/Smax)*x7(t) - 0.3*x7(t),
    x8'(t) = theta*x7(t) - 0.058*x8(t),
    x9'(t) = theta*alpha_NAb*x7(t) - 0.058*x9(t),

    y1(t) = x4(t) + x5(t) + x6(t),
    y2(t) = alpha_VLSG*(x2(t) + x3(t)),
    y3(t) = x8(t),
    y4(t) = x9(t)
)

File_Name = "EarlyViralandAntibodyDynPostInfection/Model_Identifiability_Julia/Structural_Identifiability_JointModel_VL_Ab_PostExpo_SingleComp_Results.txt"

open(File_Name,write=true) do file
    write(file,"-----------------------------------------\n")
end

open(File_Name,write=true,append=true) do file
    write(file,"-- Presentation of the model --\n")
    write(file,"-----------------------------------------\n")
    write(file,"x1 = T:     Uninfected target cells \n")
    write(file,"x2 = I1:    Unproductively infected target cells  \n")
    write(file,"x3 = I2:    Productively infected target cells\n")
    write(file,"x4 = VI:    Infectious viruses  \n")
    write(file,"x5 = VNI:   Non-infectious viruses\n")
    write(file,"x6 = Vini:  Inoculated viruses \n\n")
    write(file,"x7 = S:     SL ASCs \n")
    write(file,"x8 = BAb:   Binding antibodies\n")
    write(file,"x9 = NAb :  Neutralizing antibodies \n")

    write(file,"y1 = gRNA:    genomic viral load \n")
    write(file,"y2 = sgRNA:   subgenomic viral load \n")
    write(file,"y3 = Obs_BAb:   binding antibodies \n")
    write(file,"y4 = Obs_NAb:   neutralizing antibodies \n\n")
end
#= -------------------------------------------------------- =#




#= --- Study of the global indentifiability (without consideration for initial conditions) --- =#
open(File_Name,write=true,append=true) do file
    write(file,"-----------------------------------------\n")
    write(file,"-- Study of the global identifiability --\n")
    write(file,"-----------------------------------------\n")
end

Global_Identifiability_Results = assess_identifiability(ODE_model_Single_URT_Comp)

# Save results
for i in eachindex(Global_Identifiability_Results) 
    res = string(i," -> ",Global_Identifiability_Results[i],"\n")
    open(File_Name,write=true,append=true) do file
        write(file,res)
    end
end 
#= -------------------------------------------------------- =#



#= --- Identification of identifiable functions  --- =#
open(File_Name,write=true,append=true) do file
    write(file,"\n-----------------------------------------\n")
    write(file,"-- Identification of identifiable functions (without states)--\n")
    write(file,"-----------------------------------------\n")
end

Identiable_functions = find_identifiable_functions(ODE_model_Single_URT_Comp,with_states=true) 

# Save results
for i in eachindex(Identiable_functions) 
    res = string(i," -> ",Identiable_functions[i],"\n")
    open(File_Name,write=true,append=true) do file
        write(file,res)
    end
end 
#= -------------------------------------------------------- =#



#= --- Study of the global indentifiability (assuming initial conditions) --- =#
open(File_Name,write=true,append=true) do file
    write(file,"\n-----------------------------------------\n")
    write(file,"-- Study of the global identifiability with IC known on SL ASCs --\n")
    write(file,"-----------------------------------------\n")
end

Global_Identifiability_Results_knownIC = assess_identifiability(ODE_model_Single_URT_Comp,known_ic = [x7])

# Save results
for i in eachindex(Global_Identifiability_Results_knownIC) 
    res = string(i," -> ",Global_Identifiability_Results_knownIC[i],"\n")
    open(File_Name,write=true,append=true) do file
        write(file,res)
    end
end 
#= -------------------------------------------------------- =#