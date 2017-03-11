module HXXZ

using Combinatorics

export basis, h_ising, h_xxz, parity_eigenvalues, xrotation_eigenvalues, averageipr_with_delta, five_lowest_energy_eigenstates, magnetization, fidelity, inverse_participation_ratio

function basis(L::Int64, upspins = div(L,2))
    downspins = L - upspins
    #dim = factorial(L)/(factorial(upspins)*factorial(downspins))
    
    arrayup = ones(Int64, upspins)
    arraydown = zeros(Int64, downspins)
    
    onebasisvector = vcat(arrayup,arraydown)
    
    unique(permutations(onebasisvector))
    
end

function spin_z(s::Int64)
    if s == 1
        return 1/2
    elseif s == 0
        return -1/2
    end
end

function h_ising(base; J = 1., delta = 1., open_chain = true)
    dim = length(base)
    h = zeros(dim,dim)
    L = length(base[1])
    for i in 1:dim
        base_element = base[i]
        energy = 0.0
        for j in 1:L-1
            energy += J*delta*(spin_z(base_element[j])*spin_z(base_element[j+1]))
        end
        if !open_chain
            energy += J*delta*(spin_z(base_element[1])*spin_z(base_element[L]))
        end
        h[i,i] = energy
    end
    h
end

""" Returns the energy associated to operate with the flip-flop term a pair of indices, and it changes the state of the
array"""
function flip_flop(a::Array{Int64,1}, indices::Tuple{Int64,Int64}; J = 1)
    b = copy(a)
    
    if b[indices[1]] != b[indices[2]]
        b[indices[1]], b[indices[2]] = b[indices[2]], b[indices[1]]
        
        return b, J/2
        
    else
        return b,0.0
    end    
    
end    


function flip_flop(a::Array{Int64,1};open_c = true)
    len = length(a)
    states_and_energies = Tuple{Array{Int64,1},Float64}[]
    for i in 1:len-1
        flipped = flip_flop(a,(i,i+1))
        if flipped[2] > 0.0
            push!(states_and_energies, flipped)
        end
    end
    if !open_c
        flipped = flip_flop(a,(1,len))
         if flipped[2] > 0.0
            push!(states_and_energies, flipped)
        end
    end
        
    states_and_energies
end

function h_flipflop(basis::Array{Array{Int64,1},1}; J = 1.0, open_chain = true)
    dim = length(basis)
    h = zeros(dim,dim)
    for k in 1:dim
        states_and_energies = flip_flop(basis[k], open_c = open_chain )
        length_states = length(states_and_energies)
        if length_states > 0.
            for j in 1:length_states
                new_index = find(basis .== [states_and_energies[j][1]])
                energy_overlap = states_and_energies[j][2]
                h[k, new_index[]] = energy_overlap
            end
        end
    end
    h
end

function h_xxz(basis::Array{Array{Int64,1},1}; j = 1.0, d = 0.4, OPEN = true)
    h_ising(basis, J =j, delta =d, open_chain = OPEN) + h_flipflop(basis, J = j, open_chain = OPEN )
end

function mirror(a::Array{Int64,1})
    b = copy(a)
    L = length(b)
    if iseven(L)
        for i in 1:div(L,2)
            b[i],b[L+1-i] = b[L+1-i], b[i]
        end
    else
        for i in 1:div(L-1,2)
            b[i], b[L+1-i] = b[L+1-i], b[i]
        end
    end
    b
end


function parity_eigenvalues(basis::Array{Array{Int64,1},1},states::Matrix)
    L = length(basis)
    r = rand(1:L)
    b_pick = basis[r]
    b_mirror = mirror(b_pick)
    mirror_index = find(basis .== [b_mirror]) ## Index associated to the mirror
    eigvalues = zeros(Int64,L)
    for i in 1:L
        while abs(states[:,i][r]) <= 1.e-15
            r = rand(1:L)
            b_pick = basis[r]
            b_mirror = mirror(b_pick)
            mirror_index = find(basis .== [b_mirror]) ## Index associated to the mirror
        end
        sign(states[:,i][r][]) == sign(states[:,i][mirror_index][])  ? eigvalues[i] = 1 : eigvalues[i] = -1
    end
    eigvalues
end

function x_rotation(a::Array{Int64,1})
    b = copy(a)
    L = length(b)
    for i in 1:L
        b[i] == 1 ? b[i] = 0 : b[i] = 1
    end
    b
end

function xrotation_eigenvalues(basis::Array{Array{Int64,1},1},states::Matrix)
    L = length(basis)
    r = rand(1:L)
    b_pick = basis[r]
    b_rotated = x_rotation(b_pick)
    rot_index = find(basis .== [b_rotated]) ## Index associated to the rotation
    eigvalues = zeros(Int64,L)
    for i in 1:L
         while abs(states[:,i][r]) <= 1.e-15
            r = rand(1:L)
            b_pick = basis[r]
            b_rotated = x_rotation(b_pick)
            rot_index = find(basis .== [b_rotated]) 
        end
        sign(states[:,i][r][]) ==  sign(states[:,i][rot_index][])   ? eigvalues[i] = 1 : eigvalues[i] = -1
    end
    eigvalues
end

function participation_ratio(vec::Array{Float64,1})
    len = length(vec)
    pr = 0.0
    for i in 1:len
        pr += vec[i]^4.
    end
    pr
end

function inverse_participation_ratio(vec::Array{Float64,1})
    len = length(vec)
    pr = 0.0
    for i in 1:len
        pr += vec[i]^4.
    end
    ipr = 1./pr
end


"""Returns the average inverse participation ratio in the provided base"""
function averageipr_with_delta(deltainterval::Float64,deltamax::Float64, basis::Array{Array{Int64,1},1})
    deltarange = collect(0.0:deltainterval:deltamax)
    averages = zeros(length(deltarange))
    i = 1
    for delta in deltarange
        ham_matrix = h_xxz(basis, d = delta)
        states = eigvecs(ham_matrix)
        average_ipr = mean([inverse_participation_ratio(states[:,i]) for i in 1:length(states[1,:])])
        averages[i] = average_ipr
        i += 1
    end
    deltarange,averages
end

"""Returns the average inverse participation ratio in a basis diferent from the standard one """
function averageipr_with_delta(deltainterval::Float64,deltamax::Float64, basis::Array{Array{Int64,1},1}, newbasis::Array{Float64,2})
    deltarange = collect(0.0:deltainterval:deltamax)
    averages = zeros(length(deltarange))
    i = 1
    for delta in deltarange
        ham_matrix = h_xxz(basis, d = delta)
        states = transpose(eigvecs(ham_matrix)'*newbasis) #To put the vector as columns
        average_ipr = mean([inverse_participation_ratio(states[:,i]) for i in 1:length(states[1,:])])
        averages[i] = average_ipr
        i += 1
    end
    deltarange,averages
end

function five_lowest_energy_eigenstates(deltamin::Float64,deltainterval::Float64,deltamax::Float64, basis::Array{Array{Int64,1},1})
    deltarange = collect(deltamin:deltainterval:deltamax)
    five_lowest = zeros(5)
    for delta in deltarange
        ham_matrix = h_xxz(basis, d = delta, OPEN = false)
        temp_lowest = sort(eigvals(ham_matrix))[1:5]
        five_lowest = hcat(five_lowest, temp_lowest)
    end
    deltarange, five_lowest[:,2:end]
end


function magnetization(state::Array{Int64,1})
    spins_up = countnz(state)
    spins_down = length(state) - spins_up
    m = spins_up - spins_down
end

function magnetization(basis::Array{Array{Int64,1},1})
    mz = [magnetization(state) for state in basis] 
end

function magnetization(deltamin::Float64,deltainterval::Float64,deltamax::Float64, basis::Array{Array{Int64,1},1}; state = 0)
    deltarange = collect(deltamin:deltainterval:deltamax)
    mbasis = magnetization(basis)
    mgstate = zeros(length(deltarange))
    i = 1
    for delta in deltarange
        ham_matrix = h_xxz(basis, d = delta, OPEN = false)
        ground_state = eigvecs(ham_matrix)[:,state+1]
        mag = 1/2.*dot(ground_state.^2,mbasis)
        mgstate[i] = mag
        i += 1
    end
    deltarange, mgstate
end

function fidelity(deltamin::Float64,deltainterval::Float64,deltamax::Float64, basis::Array{Array{Int64,1},1}; state = 0)
#state = 0 indicates the groundstate
    deltarange = collect(deltamin:deltainterval:deltamax)
    fidelity = zeros(length(deltarange)-1)
    ham_matrix = h_xxz(basis, d = deltarange[1], OPEN = false)
    ground = eigvecs(ham_matrix)[:,1]
    i = 1
    for delta in deltarange[2:end]
        ham_matrix = h_xxz(basis, d = delta, OPEN = false)
        new_ground = eigvecs(ham_matrix)[:,state+1]
        fidelity[i] = dot(ground,new_ground)^2.
        ground = new_ground
        i += 1
    end
    deltarange[1:end-1], fidelity
end


end
