Nx=1
Ny_max=2

U11=1 # nearest neighbor potential in first layer
U22=1 # nearest neighbor potential term in 2nd layer
Usi=1# nearest neighbor potential betweeen hydrogen and silicon atom

t22=1 # nearest nighbor hopping in 2nd layer
t11=1 # nearest neighbor hopping in first layer
t_si= # nearest neighor hopping between silicon and hydrogen atom
t_b1=1 # loss/gain at boundary 1
t_b2=1 # loss/gain at boundary 2
t_b3=1 # loss/gain at boundary 3
t_b4=1; # loss/gain at boundary 4

it=Iterators.product(ntuple(_ -> 0:1, 1+Nx*Ny_max*6)...);
p=collect(it); # possible states with non conserving particle number

# Generating a matrix with zeros.
H=zeros(2^(1+Nx*Ny_max*6),2^(1+Nx*Ny_max*6)); # 6 nearest neighbours.
size(H)

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
layer_1

layer_2 = [];
for Ny=1:Ny_max
    for k=6*Nx*(Ny-1)+4*Nx+1:2:6*Nx*(Ny-1)+4*Nx+2*Nx-1
        push!(layer_2,[k,k+1])
    end
end
layer_2

layer_1_2 = [];
for Ny=1:Ny_max
    for n=2:2:4*Nx
        push!(layer_1_2,[Nx*6*(Ny-1)+n,Nx*6*(Ny-1)+4*Nx+Int(n/2)])
    end
end
layer_1_2

layer_2_3 = [];
for Ny=2:Ny_max
    for n1=2:2:4*Nx
        push!(layer_2_3,[Nx*6*(Ny-1)+n1,Nx*6*(Ny-2)+4*Nx+Int(n1/2)])
    end
end
layer_2_3

#= The following concatenates all the ordered pairs. =#
bonds = vcat(layer_1,layer_2,layer_1_2,layer_2_3);

length(p[7])

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
        H[n,n]+=U11*p[n][site1]*p[n][site2] # NN Interaction term
        
        
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
                if collect(p[n1])==q
                    H[n,n1]+=t22*phase1*phase2
                    H[n1,n]+=t22*phase1*phase2 # hermitian conjugate
                end
            end
        end
    end
end

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
first_bd

second_bd

third_bd = [];
fourth_bd = [];
for Ny=1:Ny_max
    k=1+Nx*6*(Ny-1)
    push!(third_bd,k)
    k=4*Nx+Nx*6*(Ny-1)
    push!(fourth_bd,k)
    #println(k)  
end
third_bd

fourth_bd

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
                H[n,n1]+=t_b1*phase2
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
                H[n,n1]+=t_b2*phase2
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
                H[n,n1]+=t_b3*phase2
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
                H[n,n1]+=t_b4*phase2
            end
        end
    end    
end     

H;
