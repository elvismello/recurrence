# This version of the script plots the entropy as a function of epsilon

import DynamicalSystems as DS
import DynamicalSystems: RecurrenceAnalysis as RA
import PyPlot as plt
import Random


sample_size = 10000 # 1% of 1e6
series_size = 1000

epsilon_list = collect(0:0.01:1)
entropy_list = []

for epsilon in epsilon_list
    print("\e[2K") # clear whole line
    print("\e[1G") # move cursor to column 1
    print("epsilon is ", epsilon)
    
    system = DS.Systems.logistic(0.4, r=4.0) # Selecting system
    #system = DS.Systems.betatransformationmap(0.25, β=2.0)
    
    trajectory = DS.trajectory(system, series_size-1) # evolving trajectory
    
    time = [0:1:series_size-1;] # "time", using a sequence from 0 to 999
    
    # Getting recurrence matrix and transforming it to a simple numeric matrix
    # populated with 1s and 0s
    matrix = RA.RecurrenceMatrix(trajectory, epsilon, parallel=true)
    grayscale = RA.grayscale(matrix)
    grayscale = trunc.(Int, grayscale)
    
    #PyPlot.imshow(grayscale, cmap="Greys_r", interpolation="none")
    #PyPlot.savefig("juliaTest.png", dpi=600)
    
    
    random_index = zeros(sample_size)
    Random.rand!(random_index)
    
    random_index = random_index * series_size^2
    round.(random_index)
    
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
    
    # TODO Divide current_found by sample_size in order to normalize them into
    # probabilities

    current_found = current_found / sample_size

    current_entropy = 0.0
    for i in current_found
        # expression defined in eq. 4
        current_entropy -= i * log(i)
    end

    push!(entropy_list, current_entropy)
end
print("\n\n")
display(entropy_list)
print("\n\n")

plt.scatter(epsilon_list, entropy_list)
plt.show()