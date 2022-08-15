# This version of the script plots the entropy as a function of epsilon
# it also finds epsilon for which entropy is maximum

# TODO make it find gamma for 100 samples

import DynamicalSystems as DS
import DynamicalSystems: RecurrenceAnalysis as RA
import PyPlot as plt
import Random
import Statistics


println("Selecting best epsilon...")
sample_size = 10000 # 1% of 1e6
series_size = 1000

epsilon_list = collect(0:0.01:1)
entropy_list = []

#system = DS.Systems.betatransformationmap(0.25, β=2.0)
system = DS.Systems.logistic(0.4, r=4.0) # Selecting system
trajectory = DS.trajectory(system, series_size-1) # evolving trajectory
time = [0:1:series_size-1;] # "time", using a sequence from 0 to 999


for epsilon in epsilon_list    
    
    # Getting recurrence matrix and transforming it to a simple numeric matrix
    # populated with 1s and 0s
    matrix = RA.RecurrenceMatrix(trajectory, epsilon, parallel=true)
    grayscale = RA.grayscale(matrix, (1, 0))
    grayscale = trunc.(Int, grayscale)
    

    # Getting indexes for random data    
    random_index = zeros(sample_size)
    Random.rand!(random_index, collect(1:(series_size)^2))
    random_index = trunc.(Int, random_index)

    x_index = zeros(sample_size)
    y_index = zeros(sample_size)
    
    for (i, index) in enumerate(random_index)
        x_index[i] = index % series_size
        y_index[i] = round(index / series_size)
    end
    
    x_index = trunc.(Int, x_index)
    y_index = trunc.(Int, y_index)
    
    # Verifying border problems
    for (i, j) in enumerate(x_index)
        if j >= series_size
            x_index[i] = series_size - 1
        elseif j == 0
            x_index[i] = 1
        end
    end
    
    for (i, j) in enumerate(y_index)
        if j >= series_size
            y_index[i] = series_size - 1
        elseif j == 0
            y_index[i] = 1
        end
    end
    
    
    # Getting matrices from recurrence matrix
    matrices = []
    for (x, y) in zip(x_index, y_index)
        push!(matrices, grayscale[x:x+1, y:y+1])
    end
    
    # Counting each type
    possible_matrices = [[0 0; 0 0], [1 0; 0 0], [0 1; 0 0], [0 0; 1 0],
                         [0 0; 0 1], [1 1; 0 0], [1 0; 1 0], [0 0; 1 1],
                         [0 1; 0 1], [1 0; 0 1], [0 1; 1 0], [1 1; 1 0],
                         [0 1; 1 1], [1 0; 1 1], [1 1; 0 1], [1 1; 1 1]]
    
    current_found = zeros(16)
    for (i, current_matrix) in enumerate(possible_matrices)
        for j in matrices
            if current_matrix == j
                current_found[i] += 1
            end
        end
    end
    
    # normalizing them into probabilities
    current_found = current_found / sample_size

    current_entropy = 0.0
    for i in current_found
        # expression defined in eq. 4
        current_entropy -= i * log(i)
    end

    push!(entropy_list, current_entropy)
end

# finding max entropy
max_entropy_epsilon = [0, 0]
for (i, j) in enumerate(entropy_list)
    if j != NaN
        if max_entropy_epsilon[2] < j
            global max_entropy_epsilon = [i, j]
        end
    end
end

max_entropy_epsilon = epsilon_list[convert(Int, max_entropy_epsilon[1])]
print("\nmax entropy is at epsilon = ", max_entropy_epsilon)

plt.scatter(epsilon_list, entropy_list)
plt.savefig("entropy(epsilon)", dpi=300)
plt.close()




# with the now known best epsilon, compute gamma
# Getting recurrence matrix
println("\nCalculating system gamma...")
matrix = RA.RecurrenceMatrix(trajectory, max_entropy_epsilon, parallel=true)
grayscale = RA.grayscale(matrix, (1, 0))
grayscale = trunc.(Int, grayscale) # not sure if it is still needed


# Getting indexes for random data    
random_index = zeros(sample_size)
Random.rand!(random_index, collect(1:(series_size)^2))
random_index = trunc.(Int, random_index)

x_index = zeros(sample_size)
y_index = zeros(sample_size)
    
for (i, index) in enumerate(random_index)
    x_index[i] = index % series_size
    y_index[i] = round(index / series_size)
end
    
x_index = trunc.(Int, x_index)
y_index = trunc.(Int, y_index)
    
# Verifying border problems
for (i, j) in enumerate(x_index)
    if j >= series_size
        x_index[i] = series_size - 1
    elseif j == 0
        x_index[i] = 1
    end
end
    
for (i, j) in enumerate(y_index)
    if j >= series_size
        y_index[i] = series_size - 1
    elseif j == 0
        y_index[i] = 1
    end
end


# Getting matrices from recurrence matrix
matrices = []
for (x, y) in zip(x_index, y_index)
    push!(matrices, grayscale[x:x+1, y:y+1])
end

# Counting each type
possible_matrices = [[0 0; 0 0], [1 0; 0 0], [0 1; 0 0], [0 0; 1 0],
                     [0 0; 0 1], [1 1; 0 0], [1 0; 1 0], [0 0; 1 1],
                     [0 1; 0 1], [1 0; 0 1], [0 1; 1 0], [1 1; 1 0],
                     [0 1; 1 1], [1 0; 1 1], [1 1; 0 1], [1 1; 1 1]]

current_found1 = zeros(16)
for (i, current_matrix) in enumerate(possible_matrices)
    for j in matrices
        if current_matrix == j
            current_found1[i] += 1
        end
    end
end


# Current system
# Calculating gamma by the paper definition
gamma1 = abs(current_found1[2] - current_found1[3]) +
         abs(current_found1[2] - current_found1[4]) +
         abs(current_found1[2] - current_found1[5]) +
         abs(current_found1[3] - current_found1[4]) +
         abs(current_found1[3] - current_found1[5]) +
         abs(current_found1[4] - current_found1[5])


gamma2 = abs(current_found1[6] - current_found1[7]) +
         abs(current_found1[6] - current_found1[8]) +
         abs(current_found1[6] - current_found1[9]) +
         abs(current_found1[7] - current_found1[8]) +
         abs(current_found1[7] - current_found1[9]) +
         abs(current_found1[8] - current_found1[9])

gamma3 = abs(current_found1[10] - current_found1[11])

gamma4 = abs(current_found1[12] - current_found1[13]) +
         abs(current_found1[12] - current_found1[14]) +
         abs(current_found1[12] - current_found1[15]) +
         abs(current_found1[13] - current_found1[14]) +
         abs(current_found1[13] - current_found1[15]) +
         abs(current_found1[14] - current_found1[15])

system_gamma = (gamma1 + gamma2 + gamma3 + gamma4) / sample_size

print("\nSystem gamma is: ")
print(system_gamma)




# Getting data for shuffled dataset
trajectory = Random.shuffle(trajectory)

matrix = RA.RecurrenceMatrix(trajectory, max_entropy_epsilon, parallel=true)
grayscale = RA.grayscale(matrix, (1, 0))
grayscale = trunc.(Int, grayscale) # not sure if it is still needed

# Getting indexes for random data    
random_index = zeros(sample_size)
Random.rand!(random_index, collect(1:(series_size)^2))
random_index = trunc.(Int, random_index)

x_index = zeros(sample_size)
y_index = zeros(sample_size)
    
for (i, index) in enumerate(random_index)
    x_index[i] = index % series_size
    y_index[i] = round(index / series_size)
end
    
x_index = trunc.(Int, x_index)
y_index = trunc.(Int, y_index)
    
# Verifying border problems
for (i, j) in enumerate(x_index)
    if j >= series_size
        x_index[i] = series_size - 1
    elseif j == 0
        x_index[i] = 1
    end
end
    
for (i, j) in enumerate(y_index)
    if j >= series_size
        y_index[i] = series_size - 1
    elseif j == 0
        y_index[i] = 1
    end
end
    

# Getting matrices from recurrence matrix
matrices_sf = []
for (x, y) in zip(x_index, y_index)
    push!(matrices_sf, grayscale[x:x+1, y:y+1])
end

# for each 100 random matrices in matrices_sf, calculate a gamma

gammas = []
for i in 100:div(sample_size, 100):sample_size
    # Selects a couple matrices, calculate gamma, then go to the next couple
    current_matrices = matrices_sf[i-99:i]

    current_found2 = zeros(16)
    for (i, current_matrix) in enumerate(possible_matrices)
        for j in current_matrices
            if current_matrix == j
                current_found2[i] += 1
            end
        end
    end


    # Calculating gamma by the paper definition
    local gamma1 = abs(current_found2[2] - current_found2[3]) +
                    abs(current_found2[2] - current_found2[4]) +
                    abs(current_found2[2] - current_found2[5]) +
                    abs(current_found2[3] - current_found2[4]) +
                    abs(current_found2[3] - current_found2[5]) +
                    abs(current_found2[4] - current_found2[5])

    local gamma2 = abs(current_found2[6] - current_found2[7]) +
                    abs(current_found2[6] - current_found2[8]) +
                    abs(current_found2[6] - current_found2[9]) +
                    abs(current_found2[7] - current_found2[8]) +
                    abs(current_found2[7] - current_found2[9]) +
                    abs(current_found2[8] - current_found2[9])

    local gamma3 = abs(current_found2[10] - current_found2[11])

    local gamma4 = abs(current_found2[12] - current_found2[13]) +
                    abs(current_found2[12] - current_found2[14]) +
                    abs(current_found2[12] - current_found2[15]) +
                    abs(current_found2[13] - current_found2[14]) +
                    abs(current_found2[13] - current_found2[15]) +
                    abs(current_found2[14] - current_found2[15])

    total_gamma = (gamma1 + gamma2 + gamma3 + gamma4) / div(sample_size, 100)
    push!(gammas, total_gamma)
end

# calculating mean and determining if data is stochastic
mean = Statistics.std(gammas)
sigma = Statistics.std(gammas)

print("\nMean gamma + 3 sigma is: ", mean + 3 * sigma)

if system_gamma < mean + 3 * sigma
    print("\nThe system is stochastic\n")
else
    print("\nThe system is deterministic\n")
end
