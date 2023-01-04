using LinearAlgebra
using SparseArrays
using DifferentialEquations
using PyCall
using Random

#= Planck's constant.=#
hbar = 1.054*10^(-34);
#= Mass of an electron.=#
m_e = 9.11*10^(-31);


#= Silicon hydrogen NN interaction.=#
C_1 = 0.1*1.6*10^(-19);
U_Si_H(x) = C_1/x;
d_U_Si_H(x) = -C_1/x^2; #=Derivative.=#

#=Hopping between silicon atom and the hydrogen atom.=#
C_2 = 2.3*1.6*10^(-19);
x_0 = 1.5*10^(-10);
d = 1;
t_Si_H(x) = C_2*exp(-(x-x_0)/d);
d_t_Si_H(x) = -(C_2/d)*exp(-(x-x_0)/d); #=Derivative.=#

#=Hopping between the STM tip and the hydrogen atom.=#
C_3 = 2.3*1.6*10^(-19);
x_0 = 1.5*10^(-10);
Xi = 1;
t_STM_H(x) = C_3*exp(-(x_0-x)/Xi);
d_t_STM_H(x) = (C_3/Xi)*exp(-(x_0-x)/Xi);#=Derivative.=#

Nx = 1
Ny_max = 1

U11 = 0.16*1.6*10^(-19); # nearest neighbor potential in first layer
U22 = 0.16*1.6*10^(-19) ;# nearest neighbor potential term in 2nd layer
H_U_Si_H = 0.1*10^(-19); # nearest neighbor potential betweeen hydrogen and silicon atom

t11 = 0.8*1.6*10^(-19); # nearest neighbor hopping in first layer
t22 = 0.8*1.6*10^(-19) # nearest nighbor hopping in 2nd layer

t_si = 1.6*10^(-19) # nearest neighor hopping between silicon and hydrogen atom
t_b1 = 1.6*10^(-19)  # loss/gain at boundary 1
t_b2 = 1.6*10^(-19)  # loss/gain at boundary 2
t_b3 = 1.6*10^(-19)  # loss/gain at boundary 3
t_b4 = 1.6*10^(-19) ; # loss/gain at boundary 4
t_stm = 1.6*10^(-19) ; # hopping between the STM tip and the hydrogen atom.

#= The following is the potential in which the hydrogen experiences from
the bulk and the surrounding silicon atoms. For now, we will approximate
it as a harmonic potential. =#
k = 1.6*10^(-19);
# x_0 is defined above.
V(x) = 0.5*k*(x-x_0/2)^2;

#=
In the following lines, specify the position
of the silicon atom and the hydrogen atom in the lattice that we are interested in.
=#

Si_position = 6;
H_position = 7;

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
    local H1 = zeros(2^(1+Nx*Ny_max*6),2^(1+Nx*Ny_max*6)); # 6 nearest neighbours.
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
                    H1[n,n1]+=t_b1*phase2
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
                    H1[n,n1]+=t_b2*phase2
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
                    H1[n,n1]+=t_b3*phase2
                end
            end
        end_
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
                    H1[n,n1]+=t_b4*phase2
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
            H1[n,n] += U11 * p[n][site1] * p[n][site2] # NN Interaction term
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
                        H1[n,n1]+=t22*phase1*phase2
                        H1[n1,n]+=t22*phase1*phase2 # hermitian conjugate
                    end
                end
            end
        end
    end    
    return H1
end    

#=
The position of the silicon and the hydrogen attached to it is
specified at the top.
=#
#= Array holding the coordinates of the non-zero terms. =#
t_Si_H_Positions = []
#= Array holding the corresponding phases in the position of the elements in the array above. =#
t_Si_H_Phases = []

#H_t_Si_H = zeros(2^(1+Nx*Ny_max*6),2^(1+Nx*Ny_max*6))
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
                #H_t_Si_H[n,n1] += phase1*phase2
                #H_t_Si_H[n1,n] += phase1*phase2
                push!(t_Si_H_Positions,[n,n1])
                #push!(t_Si_H_Positions,[n1,n])
                push!(t_Si_H_Phases,[phase1*phase2])
                # The Hamiltonian.
                #H_array[n,n1] = x -> t_Si_H(x)*phase1*phase2
                #H_array[n1,n] = x -> t_Si_H(x)*phase1*phase2
            end        
        end
    end    
end

#=
The position of the silicon and the hydrogen attached to it is
specified at the top.
=#
#H_U_Si_H = zeros(2^(1+Nx*Ny_max*6),2^(1+Nx*Ny_max*6))
U_Si_H_Positions = []
for n = 1:2^(1+Nx*Ny_max*6) # Iterating over all basis states.
    
    #=  
    p[n][Si_position] = 1 means there is an electron in the Si atom.
    p[n][H_position] = 0 means there is no electron in the H atom.
    =#
    if p[n][H_position]==1 && p[n][Si_position]==1
        #H_U_Si_H[n,n] += 1
        push!(U_Si_H_Positions,[n,n])
        #H[n,n] += H_U_Si_H
        #dH[n,n] *= "a"
    end    
end

#=
Exchange of electrons between the STM tip and the hydrogen atom.
This is a non-Hermitian term in the Hamiltonian. The position of
the hydrogen is specified at the top.
=#
#H_t_STM_H = zeros(2^(1+Nx*Ny_max*6),2^(1+Nx*Ny_max*6))
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
        
    for n1=1:2^(Nx*Ny_max*6)
        if collect(p[n1])==q
            #H_t_STM_H[n,n1] += 1
            push!(t_STM_H_Positions,[n,n1])
            #H[n,n1] += t_stm
            #dH_dx[n,n1] *= "c"
        end
    
    end
    
end    

HC = Hamiltonian_constant();
function Hamiltonian_variable(x)
    x = real(x)
    local H_x = zeros(2^(1+Nx*Ny_max*6),2^(1+Nx*Ny_max*6));
    #= Nearest neighbour interaction between Si and H atom. =#
    N1 = length(t_Si_H_Positions)
    for i in 1:N1
        H_x[t_Si_H_Positions[i][1],t_Si_H_Positions[i][2]] += U_Si_H(x)*t_Si_H_Phases[i][1]
        H_x[t_Si_H_Positions[i][2],t_Si_H_Positions[i][1]] += U_Si_H(x)*t_Si_H_Phases[i][1] # Hermitian conjugate.
    end
    #= NN hoppoing between the Si and H atom. =#
    N2 = length(U_Si_H_Positions)
    for i in 1:N2
        H_x[U_Si_H_Positions[i][1],U_Si_H_Positions[i][2]] += t_Si_H(x)
    end
    #= STM and the H atom interaction.=#
    N3 = length(t_STM_H_Positions)
    for i in 1:N3
        H_x[t_STM_H_Positions[i][1],t_STM_H_Positions[i][2]] += t_STM_H(x) 
    end
    return HC+H_x
end

function dHamiltonian(x)
    x = real(x)
    local dHx = zeros(2^(1+6*Nx*Ny_max),2^(1+6*Nx*Ny_max));
    #= Nearest neighbour interaction between Si and H atom. =#
    N1 = length(t_Si_H_Positions)
    for i in 1:N1
        dHx[t_Si_H_Positions[i][1],t_Si_H_Positions[i][2]] += d_U_Si_H(x)*t_Si_H_Phases[i][1]
        dHx[t_Si_H_Positions[i][2],t_Si_H_Positions[i][1]] += d_U_Si_H(x)*t_Si_H_Phases[i][1] # Hermitian conjugate.
    end    
    #= NN hopping between the Si and H atom. =#
    N2 = length(U_Si_H_Positions)
    for i in 1:N2
        dHx[U_Si_H_Positions[i][1],U_Si_H_Positions[i][2]] += d_t_Si_H(x)
    end    
    #= STM and the H atom interaction. =#
    N3 = length(t_STM_H_Positions)
    for i in 1:N3
        dHx[t_STM_H_Positions[i][1],t_STM_H_Positions[i][2]] += d_t_STM_H(x) 
    end    
    return dHx
end

Random.seed!(3000)
psi = rand(Float64,(1,2^(1+Nx*Ny_max*6)));

py"""
f = open('dynamics_data'+'.txt', 'w')
def Write_file(t, x, p):
    f = open('dynamics_data'+'.txt', 'a')
    f.write(str(t) +'\t'+ str(x)+ '\t' + str(p) +'\n')
"""

#= Initial conditions. =#
TT = parse(Float64,ARGS[1]);
t_i = 0.0
x_i = (1.5*10^(-10))/2
p_i = hbar/x_i
#= List to store the t,y and z values. =#
ts = [t_i]
xs = [x_i]
ps = [p_i]
# Time steps. =#
dt = 10^(-16)
# Final time. #
t_end = TT*10^(-13);

#= Initializing the parameters. =#
t = t_i
x = x_i
p = p_i
#= The wavefunction is started with a matrix of random numbers.=#
#psi = rand(Float64,(1,2^(1+Nx*Ny_max*6)));
#= The wavefunction is normalized. =#
Psi_i = psi/norm(psi);

# The coupled ODEs. =#
dxdt(t,x,p) = p/m_e;
dpdt(t,x,p) = -k*(x-x_0/2)-(Psi'*dHamiltonian(x)*Psi)[1];

Psi = Psi_i' #= Psi is a column matrix. =#
while t<t_end

    #= Runge Kutta algorithm of order four. =#
    k1 = dt*dxdt(t,x,p)
    h1 = dt*dpdt(t,x,p)
    k2 = dt*dxdt(t+dt/2, x+k1/2, p+h1/2)
    h2 = dt*dpdt(t+dt/2, x+k1/2, p+h1/2)
    k3 = dt*dxdt(t+dt/2, x+k2/2, p+h2/2)
    h3 = dt*dpdt(t+dt/2, x+k2/2, p+h2/2)
    k4 = dt*dxdt(t+dt, x+k3, p+h3)
    h4 = dt*dpdt(t+dt, x+k3, p+h3)

    global x = x + real(k1+2*k2+2*k3+k4)/6
    global p = p + real(h1+2*h2+2*h3+h4)/6
    global t = t + dt
    push!(ts,t)
    push!(xs,x)
    push!(ps,p)
    py"Write_file"(t,x,p)
    #= The wavefunction at time t+dt. =#
    global Psi = exp(-1im*Hamiltonian_variable(x)*dt/hbar)*Psi
    #println(t)
    #println(x)
    #println(p)
    #println(Psi);
end
#using Plots
#plot(ts,xs)
#plot!(ts,ps)

ts;

ts = 10^(15).*ts; # femtosecond
#10^(10).*xs; # angstrom
ps = 10^(15).*ps;# kg m/fs

#=
using Plots
plot(ts,xs)
plot!(xlabel="t (femtosecond)")
plot!(ylabel="x(t) (metres)")

plot(ts,ps)
plot!(xlabel="t (femtosecond)")
plot!(ylabel="p(t) (kg metres/femtosecond)")
=#

#=
for i=1:length(x)
        py"Write_file"(ts[i],xs[i],ps[i])
end
=#


#Psi = Psi_i';

#dt = 10^(-17)
#@time HH = Hamiltonian_variable(0.1)*(dt)/hbar;

#HE = exp(-1im*HH);


