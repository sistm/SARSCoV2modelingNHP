using StructuralIdentifiability
using FileIO


#= Definition of the ODE model using the @ODEmodel macro 
--------------------------------------------------------
Presentation of the model:
x1 = TN:     Uninfected target cells - in nasopharynx
x2 = I1N:    Unproductively infected target cells - in nasopharynx
x3 = I2N:    Productively infected target cells - in nasopharynx
x4 = VIN:    Infectious viruses - in nasopharynx
x5 = VNIN:   Non-infectious viruses - in nasopharynx
x6 = ViniN:  Inoculated viruses - in nasopharynx

x7 = TT:     Uninfected target cells - in trachea
x8 = I1T:    Unproductively infected target cells - in trachea
x9 = I2T:    Productively infected target cells - in trachea
x10 = VIT:    Infectious viruses - in trachea
x11 = VNIT:   Non-infectious viruses - in trachea
x12 = ViniT:  Inoculated viruses - in trachea

x13 = S:     SL ASCs
x14 = BAb:   Binding antibodies
x15 = NAb :  Neutralizing antibodies

y1 = gRNAN:    genomic viral load - in nasopharynx
y2 = sgRNAN:   subgenomic viral load - in nasopharynx
y3 = gRNAT:    genomic viral load - in trachea
y4 = sgRNAT:   subgenomic viral load - in trachea
y5 = Obs_BAb:   binding antibodies
y6 = Obs_NAb:   neutralizing antibodies
=# 

ODE_model_Two_URT_Comp = @ODEmodel(
    x1'(t) = -beta*(1 - eta*x15(t)/(1+eta*x15(t)))*x1(t)*x4(t) - 1e-3*beta*(1 - eta*x15(t)/(1+eta*x15(t)))*x1(t)*x6(t),
    x2'(t) = beta*(1 - eta*x15(t)/(1+eta*x15(t)))*x1(t)*x4(t) + 1e-3*beta*(1 - eta*x15(t)/(1+eta*x15(t)))*x1(t)*x6(t) - 3*x2(t),
    x3'(t) = 3*x2(t) - delta*x3(t),
    x4'(t) = PN*1e-3*x3(t) - 3*x4(t) - beta*(1 - eta*x15(t)/(1+eta*x15(t)))*x1(t)*x4(t),
    x5'(t) = PN*(1-1e-3)*x3(t) - 3*x5(t),
    x6'(t) = -20*x6(t) - 1e-3*beta*(1 - eta*x15(t)/(1+eta*x15(t)))*x1(t)*x6(t),

    x7'(t) = -beta*(1 - eta*x15(t)/(1+eta*x15(t)))*x7(t)*x10(t) - 1e-3*beta*(1 - eta*x15(t)/(1+eta*x15(t)))*x7(t)*x12(t),
    x8'(t) = beta*(1 - eta*x15(t)/(1+eta*x15(t)))*x7(t)*x10(t) + 1e-3*beta*(1 - eta*x15(t)/(1+eta*x15(t)))*x7(t)*x12(t) - 3*x8(t),
    x9'(t) = 3*x8(t) - delta*x9(t),
    x10'(t) = PT*1e-3*x9(t) - 3*x10(t) - beta*(1 - eta*x15(t)/(1+eta*x15(t)))*x7(t)*x10(t),
    x11'(t) = PT*(1-1e-3)*x9(t) - 3*x11(t),
    x12'(t) = -20*x12(t) - 1e-3*beta*(1 - eta*x15(t)/(1+eta*x15(t)))*x7(t)*x12(t),

    x13'(t) = 1/100000*(x4(t) + x5(t) + x10(t) + x11(t))^1/4 + rho*(1-x13(t)/Smax)*x13(t) - 0.3*x13(t),
    x14'(t) = theta*x13(t) - 0.058*x14(t),
    x15'(t) = theta*alpha_NAb*x13(t) - 0.058*x15(t),

    y1(t) = x4(t) + x5(t) + x6(t),
    y2(t) = alpha_VLSG*(x2(t) + x3(t)),
    y3(t) = x10(t) + x11(t) + x12(t),
    y4(t) = alpha_VLSG*(x8(t) + x9(t)),
    y5(t) = x14(t),
    y6(t) = x15(t)
)

File_Name = "EarlyViralandAntibodyDynPostInfection/Model_Identifiability_Julia/Structural_Identifiability_JointModel_VL_Ab_PostExpo_Fullmodel_Results.txt"

open(File_Name,write=true) do file
    write(file,"-----------------------------------------\n")
end
open(File_Name,write=true,append=true) do file
    write(file,"-- Presentation of the model --\n")
    write(file,"-----------------------------------------\n")
    write(file,"x1 = TN:     Uninfected target cells - in nasopharynx\n")
    write(file,"x2 = I1N:    Unproductively infected target cells - in nasopharynx\n")
    write(file,"x3 = I2N:    Productively infected target cells - in nasopharynx\n")
    write(file,"x4 = VIN:    Infectious viruses  - in nasopharynx\n")
    write(file,"x5 = VNIN:   Non-infectious viruses - in nasopharynx\n")
    write(file,"x6 = ViniN:  Inoculated viruses - in nasopharynx\n\n")
    write(file,"x7 = TT:     Uninfected target cells - in trachea\n")
    write(file,"x8 = I1T:    Unproductively infected target cells - in trachea \n")
    write(file,"x9 = I2T:    Productively infected target cells - in trachea\n")
    write(file,"x10 = VIT:    Infectious viruses - in trachea \n")
    write(file,"x11 = VNIT:   Non-infectious viruses - in trachea\n")
    write(file,"x12 = ViniT:  Inoculated viruses - in trachea\n\n")

    write(file,"x13 = S:     SL ASCs \n")
    write(file,"x14 = BAb:   Binding antibodies\n")
    write(file,"x15 = NAb :  Neutralizing antibodies \n")

    write(file,"y1 = gRNAN:    genomic viral load - in nasopharynx\n")
    write(file,"y2 = sgRNAN:   subgenomic viral load - in nasopharynx\n")
    write(file,"y3 = gRNAT:    genomic viral load - in trachea\n")
    write(file,"y4 = sgRNAT:   subgenomic viral load - in trachea\n")
    write(file,"y3 = Obs_BAb:   binding antibodies \n")
    write(file,"y4 = Obs_NAb:   neutralizing antibodies \n\n")
end
#= -------------------------------------------------------- =#


#= --- Study of the local indentifiability --- =#
open(File_Name,write=true,append=true) do file
    write(file,"-----------------------------------------\n")
    write(file,"-- Study of the local identifiability --\n")
    write(file,"-----------------------------------------\n")
end


Local_Identifiability_Results = assess_local_identifiability(ODE_model_Two_URT_Comp)

# Save results
for i in eachindex(Local_Identifiability_Results) 
    res = string(i," -> ",Local_Identifiability_Results[i],"\n")
    open(File_Name,write=true,append=true) do file
        write(file,res)
    end
end 
#= -------------------------------------------------------- =#