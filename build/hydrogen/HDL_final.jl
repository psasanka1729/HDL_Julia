using LinearAlgebra
using SparseArrays
using DifferentialEquations
using PyCall

#= Planck's constant.=#
hbar = 1.054*10^(-34);
#= Mass of an electron.=#
m_e = 9.11*10^(-31);

#= The following is the potential in which the hydrogen experiences from
the bulk and the surrounding silicon atoms. For now, we will approximate
it as a harmonic potential. =#

x_0 = 1.5*10^(-10);
k_potential = parse(Float64,ARGS[1]);
V(x) = (1/2)*k_potential*(x-x_0/2)^2;
dVdx(x) = k_potential*(x-x_0/2);

xn = 4*10^(-10)

#= Silicon hydrogen NN interaction.=#
C_1 = 6*1.6*10^(-31);
U_Si_H(x)= C_1/x
d_U_Si_H(x) = -C_1/x^2; #=Derivative.=#
U_Si_H(xn)/hbar

#=Hopping between silicon atom and the hydrogen atom.=#
C_2 = 2.3*1.6*10^(-19);
d = 10^(-10);
t_Si_H(x) = C_2*exp(-(x-x_0)/d);
d_t_Si_H(x) = -(C_2/d)*exp(-(x-x_0)/d); #=Derivative.=#
t_Si_H(xn)/hbar

#=Hopping between the STM tip and the hydrogen atom.=#
C_3 = 2.3*1.6*10^(-19);
Xi = 10^(-10);
t_STM_H(x) = C_3*exp(-(x_0-x)/Xi);
d_t_STM_H(x) = (C_3/Xi)*exp(-(x_0-x)/Xi);#=Derivative.=#
t_STM_H(xn)/hbar

Nx = 1
Ny_max = 2

U11 = 0.16*1.6*10^(-19); # nearest neighbor potential in first layer
U22 = 0.13*1.6*10^(-19) ;# nearest neighbor potential term in 2nd layer

t11 = 0.8*1.6*10^(-19); # nearest neighbor hopping in first layer
t22 = 0.8*1.6*10^(-19) # nearest nighbor hopping in 2nd layer

t_b1 = 1.5*1.6*10^(-19)  # loss/gain at boundary 1
t_b2 = 1.5*1.6*10^(-19)  # loss/gain at boundary 2
t_b3 = 1.5*1.6*10^(-19)  # loss/gain at boundary 3
t_b4 = 1.5*1.6*10^(-19) ; # loss/gain at boundary 4

#= The following is the potential in which the hydrogen experiences from
the bulk and the surrounding silicon atoms. For now, we will approximate
it as a harmonic potential. =#
#k = 1;
# x_0 is defined above.
#V(x) = 0.5*k*(x-x_0/2)^2;

#=
In the following lines, specify the position
of the silicon atom and the hydrogen atom in the lattice that we are interested in.
=#

Si_position = 6;
H_position = 13;


it=Iterators.product(ntuple(_ -> 0:1, 1+Nx*Ny_max*6)...);
p=collect(it); # possible states with non conserving particle number

# Generating a matrix with zeros.
H = zeros(2^(1+Nx*Ny_max*6),2^(1+Nx*Ny_max*6)); # 6 nearest neighbours.
size(H);

#= The derivative of the Hamiltonian with respect to x. =#
d_H = zeros(2^(1+Nx*Ny_max*6),2^(1+Nx*Ny_max*6));

#= List of bonds for Nx = 1 and Ny = 2 for reference.=#
bonds = [
        # first layer.
        [1,2],[2,3],[3,4],[7,8],[8,9],[9,10],
        # second layer.
        [5,6],[11,12],
        # first to second layer.
        [2,5],[4,6],[8,11],[10,12],
        # second to third layer.
        [8,5],[10,6]
        ];

layer_1 = [];
for Ny=1:Ny_max
    for k=6*(Nx)*(Ny-1)+1:6*(Nx)*(Ny-1)+4*Nx-1
        push!(layer_1,[k,k+1])
    end
end
#layer_1;

layer_2 = [];
for Ny=1:Ny_max
    for k=6*Nx*(Ny-1)+4*Nx+1:2:6*Nx*(Ny-1)+4*Nx+2*Nx-1
        push!(layer_2,[k,k+1])
    end
end
#layer_2

layer_1_2 = [];
for Ny=1:Ny_max
    for n=2:2:4*Nx
        push!(layer_1_2,[Nx*6*(Ny-1)+n,Nx*6*(Ny-1)+4*Nx+Int(n/2)])
    end
end
#layer_1_2

layer_2_3 = [];
for Ny=2:Ny_max
    for n1=2:2:4*Nx
        push!(layer_2_3,[Nx*6*(Ny-1)+n1,Nx*6*(Ny-2)+4*Nx+Int(n1/2)])
    end
end
layer_2_3;

#= The following concatenates all the ordered pairs. =#
bonds = vcat(layer_1,layer_2,layer_1_2,layer_2_3);

first_bd = [];
second_bd = [];
for Ny=1:Ny_max
    if Ny == 1 # First boundary
        
        for k=6*(Nx)*(Ny-1)+1:6*(Nx)*(Ny-1)+4*Nx
            push!(first_bd,k)
        end
        
    elseif Ny == Ny_max # second boundary.
        
        for k=6*(Nx)*(Ny-1)+4*Nx+1:2:6*(Nx)*(Ny-1)+4*Nx+2*Nx-1
            push!(second_bd,k)
            push!(second_bd,k+1)
        end   
    end
end

third_bd = [];
fourth_bd = [];
for Ny=1:Ny_max
    k=1+Nx*6*(Ny-1)
    push!(third_bd,k)
    k=4*Nx+Nx*6*(Ny-1)
    push!(fourth_bd,k)
    #println(k)  
end

function Hamiltonian_constant()
    local H_c = zeros(2^(1+Nx*Ny_max*6),2^(1+Nx*Ny_max*6)); # 6 nearest neighbours.
    #= Gain and loss of fermions at the four boundaries. =#
    for n=1:2^(Nx*Ny_max*6+1)
        # first boundary.
        for j=1:length(first_bd)
            site=first_bd[j]
            # gain
            if p[n][site]==1  # electron is there. # 1 if gain
                q=collect(p[n])  
                q[site]=0 # electron is lost. # 0 if gain 1 if loss
            # loss    
            elseif p[n][site] == 0 # electron is not there. # 0 if loss
                q=collect(p[n])
                q[site]=1 # electron is gained. # 0 if gain 1 if loss
            end                
            #= The phase is calculated below. =#
            #phase1=(-1)^(sum(p[n][k:Nx*Ny_max*6]))
            phase2=1
            if site<Nx*Ny_max*6
                phase2=(-1)^(sum(p[n][site:6*Nx*Ny_max+1]))
            end
            for n1=1:2^(1+Nx*Ny_max*6)
                if collect(p[n1])==q                
                    H_c[n,n1]+=t_b1*phase2
                end
            end                
        end
        # second boundary.
        for j=1:length(second_bd) 
            site = second_bd[j]
            # loss
            if p[n][site]==1 # electron is there.
                q=collect(p[n])
                q[site]=0 # electron is lost.
            # gain   
            elseif p[n][site]==0 # electron is not there.
                q=collect(p[n])
                q[site]=1 # electron is gained.
            end    
            phase2=1
            if site<6*Nx*Ny_max+1
                phase2=(-1)^(sum(p[n][site:6*Nx*Ny_max+1]))
            end
            for n1=1:2^(6*Nx*Ny_max+1)
                if collect(p[n1])==q
                    H_c[n,n1]+=t_b2*phase2
                end
            end
        end
        # third boundary.
        for j=1:length(third_bd) 
            site = third_bd[j]
            # loss
            if p[n][site]==1 # electron is there.
                q=collect(p[n])
                q[site]=0 # electron is lost.
            # gain   
            elseif p[n][site]==0 # electron is not there.
                q=collect(p[n])
                q[site]=1 # electron is gained.
            end  
            phase2=1
            if site<6*Nx*Ny_max+1
                phase2=(-1)^(sum(p[n][site:6*Nx*Ny_max+1]))
            end
            for n1=1:2^(6*Nx*Ny_max+1)
                if collect(p[n1])==q
                    H_c[n,n1]+=t_b3*phase2
                end
            end
        end
        # fourth boundary.
        for j=1:length(fourth_bd) 
            site = fourth_bd[j]
            # loss
            if p[n][site]==1 # electron is there.
                q=collect(p[n])
                q[site]=0 # electron is lost.
            # gain   
            elseif p[n][site]==0 # electron is not there.
                q=collect(p[n])
                q[site]=1 # electron is gained.
            end  
            phase2=1
            if site<6*Nx*Ny_max+1
                phase2=(-1)^(sum(p[n][site:6*Nx*Ny_max+1]))
            end
            for n1=1:2^(6*Nx*Ny_max+1)
                if collect(p[n1])==q
                    H_c[n,n1]+=t_b4*phase2
                end
            end
        end    
    end  
    
    #= NN interaction and NN hopping between Silicon layers. =#
    # Imagine I'm given a list bonds = [[1,2],[2,3],etc]

    for n=1:2^(1+Nx*Ny_max*6) # Loop over basis states
        for j=1:length(bonds)
            site1=bonds[j][1] # First number in the pair.
            site2=bonds[j][2] # Second number in the pair.
            #=
            Assuming NN interaction between any two layers are same (U11).
            NN interaction term has no Hermitian conjugate.
            If both p[n][site] and p[n][site2] are occupied, then
            the potential U11 will be added to the Hamiltonian. 
            =#
            H_c[n,n] += U11 * p[n][site1] * p[n][site2] # NN Interaction term
            #= NN hopping has Hermitian conjugate =#
            if p[n][site1] == 1 && p[n][site2]==0 # NN hopping
                q=collect(p[n]) # collect makes p[n] a vector.
                q[site1]=0
                q[site2]=1
                phase1=(-1)^(sum(p[n][site1:6*Nx*Ny_max+1]))
                phase2=1
                if site2 < 6*Nx*Ny_max+1 
                    phase2=(-1)^(sum(p[n][site2:6*Nx*Ny_max+1]))
                end
                for n1=1:2^(6*Nx*Ny_max+1)
                    if collect(p[n1]) == q
                        H_c[n,n1]+=t22*phase1*phase2
                        H_c[n1,n]+=t22*phase1*phase2 # hermitian conjugate
                    end
                end
            end
        end
    end    
    return H_c
#X0 = X0*10^10;
end;    


#=
The position of the silicon and the hydrogen attached to it is
specified at the top.
=#
#= Array holding the coordinates of the non-zero terms. =#
t_Si_H_Positions = []
#= Array holding the corresponding phases in the position of the elements in the array above. =#
t_Si_H_Phases = []


for n = 1:2^(1+Nx*Ny_max*6) # Iterating over all basis states.
    #=  
    p[n][Si_position] = 1 means there is an electron in the Si atom.
    p[n][H_position] = 0 means there is no electron in the H atom.
    =#
    if p[n][H_position]==1 && p[n][Si_position]==0
        
        q=p[n]
        q=collect(p[n])
        
        q[H_position]=0 # The electron in the Si is lost.
        q[Si_position]=1 # An electron in the H is gained.
        
        phase1=(-1)^(sum(p[n][H_position:1+Nx*Ny_max*6]))
        phase2=1
        
        if H_position < Nx*Ny_max*6      
            phase2=(-1)^(sum(p[n][Si_position:6*Nx*Ny_max+1]))
        end
    
        for n1=1:2^(Nx*Ny_max*6)
            if collect(p[n1])==q
                
                push!(t_Si_H_Positions,[n,n1])
                push!(t_Si_H_Phases,[phase1*phase2])

            end        
        end
    end    
end

#=
The position of the silicon and the hydrogen attached to it is
specified at the top.
=#

U_Si_H_Positions = []
for n = 1:2^(1+Nx*Ny_max*6) # Iterating over all basis states.
    #=  
    p[n][Si_position] = 1 means there is an electron in the Si atom.
    p[n][H_position] = 0 means there is no electron in the H atom.
    =#
    if p[n][H_position]==1 && p[n][Si_position]==1
        push!(U_Si_H_Positions,[n,n])
    end    
end

#=
Exchange of electrons between the STM tip and the hydrogen atom.
This is a non-Hermitian term in the Hamiltonian. The position of
the hydrogen is specified at the top.
=#

t_STM_H_Positions = []
for n = 1:2^(1+Nx*Ny_max*6) # Iterating over all basis states.
    
    # gain
    if p[n][H_position]==1  # electron is there.
        
        q=collect(p[n])  
        q[H_position]=0 # electron is lost.

    # loss    
    elseif p[n][H_position] == 0 # electron is not there.
        
        q=collect(p[n])
        q[H_position]=1 # electron is gained.
        
    end  
    # Iterating over all basis states for appropriate configurations.
    for n1=1:2^(Nx*Ny_max*6)
        if collect(p[n1])==q
            push!(t_STM_H_Positions,[n,n1])
        end
    end
end    

HC = Hamiltonian_constant();
function Hamiltonian_variable(x1)
    x1 = real(x1)
    H_x = zeros(2^(1+Nx*Ny_max*6),2^(1+Nx*Ny_max*6));
    
    #= Nearest neighbour hopping between Si and H atom. =# 
    N1 = length(t_Si_H_Positions)
    for i in 1:N1
        H_x[t_Si_H_Positions[i][1],t_Si_H_Positions[i][2]] += t_Si_H(x1)*t_Si_H_Phases[i][1]
        H_x[t_Si_H_Positions[i][2],t_Si_H_Positions[i][1]] += t_Si_H(x1)*t_Si_H_Phases[i][1] # Hermitian conjugate.
    end
    
    #= Nearest neighbour interaction between the Si and H atom. =#
    N2 = length(U_Si_H_Positions)
    for i in 1:N2
        H_x[U_Si_H_Positions[i][1],U_Si_H_Positions[i][2]] += U_Si_H(x1) # No Hermitian conjugate.
    end
    
    #= STM and the H atom hopping.=#
    N3 = length(t_STM_H_Positions)
    for i in 1:N3
        H_x[t_STM_H_Positions[i][1],t_STM_H_Positions[i][2]] += t_STM_H(x1) # No Hermitian conjugate.
    end
    
    return HC+H_x
end;

function Hamiltonian_variable_without_STM(x1)
    x1 = real(x1)
    H_x = zeros(2^(1+Nx*Ny_max*6),2^(1+Nx*Ny_max*6));
    
    #= Nearest neighbour hopping between Si and H atom. =# 
    N1 = length(t_Si_H_Positions)
    for i in 1:N1
        H_x[t_Si_H_Positions[i][1],t_Si_H_Positions[i][2]] += t_Si_H(x1)*t_Si_H_Phases[i][1]
        H_x[t_Si_H_Positions[i][2],t_Si_H_Positions[i][1]] += t_Si_H(x1)*t_Si_H_Phases[i][1] # Hermitian conjugate.
    end
    
    #= Nearest neighbour interaction between the Si and H atom. =#
    N2 = length(U_Si_H_Positions)
    for i in 1:N2
        H_x[U_Si_H_Positions[i][1],U_Si_H_Positions[i][2]] += U_Si_H(x1) # No Hermitian conjugate.
    end
    
    return HC+H_x
end;

function dHamiltonian(x)
    x = real(x)
    local dHx = zeros(2^(1+6*Nx*Ny_max),2^(1+6*Nx*Ny_max));
    
    #= Nearest neighbour hopping between Si and H atom. =#
    N1 = length(t_Si_H_Positions)
    for i in 1:N1
        dHx[t_Si_H_Positions[i][1],t_Si_H_Positions[i][2]] += d_t_Si_H(x)*t_Si_H_Phases[i][1]
        dHx[t_Si_H_Positions[i][2],t_Si_H_Positions[i][1]] += d_t_Si_H(x)*t_Si_H_Phases[i][1] # Hermitian conjugate.
    end 
    
    #= Nearest neighbour interaction between the Si and H atom. =#
    N2 = length(U_Si_H_Positions)
    for i in 1:N2
        dHx[U_Si_H_Positions[i][1],U_Si_H_Positions[i][2]] += d_U_Si_H(x)
    end  
    
    #= STM and the H atom interaction. =#
    N3 = length(t_STM_H_Positions)
    for i in 1:N3
        dHx[t_STM_H_Positions[i][1],t_STM_H_Positions[i][2]] += d_t_STM_H(x) 
    end    
    return dHx
end;


py"""
f = open('force_data'+'.txt', 'w')
def Write_file_force(x, force):
    f = open('force_data'+'.txt', 'a')
    f.write(str(x)+ '\t' + str(force) +'\n')
"""



#X0 = [0.6,0.7,0.8,1,2,3,4,5]*10^(-10)
X0 = LinRange(1,4,5)*10^(-10)
Force = []
F(x1,Psi1) = -(Psi1'*dHamiltonian(x1)*Psi1)[1]-dVdx(x1)
for xs in X0
    ED = eigen(Hamiltonian_variable_without_STM(xs));
    Eigenvalues = ED.values;
    Eigenvectors = ED.vectors;
    Max_eigenvalue_index = findall(x->imag(x)==maximum(imag(Eigenvalues)), Eigenvalues);
    Max_eigenvalue = Eigenvalues[Max_eigenvalue_index[1]];
    Max_eigenvector = Eigenvectors[1:2^(1+Nx*Ny_max*6),Max_eigenvalue_index[1]:Max_eigenvalue_index[1]];
    ff = F(xs,Max_eigenvector)
    py"Write_file_force"(xs,ff)
    push!(Force,real(ff))
    #println(ff)
end
