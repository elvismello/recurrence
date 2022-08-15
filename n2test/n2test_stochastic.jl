# TODO generate time series --------------Done

# TODO generate recurrence matrix --------Done

# TODO get results? ----------------------Done?



import DynamicalSystems as DS
import DynamicalSystems: RecurrenceAnalysis as RA
import PyPlot
import Random


epsilon = 0.28 # found empiricaly for Random.rand! filling matrices
sample_size = 100000
series_size = 1000



trajectory = zeros(series_size)
Random.rand!(trajectory)

time = [0:1:series_size-1;] # "time", using a sequence from 0 to 999

# Getting recurrence matrix and transforming it to a simple numeric matrix
# populated with 1s and 0s
matrix = RA.RecurrenceMatrix(trajectory, epsilon, parallel=true)
grayscale = RA.grayscale(matrix)
grayscale = trunc.(Int, grayscale)

PyPlot.scatter(time, trajectory)
PyPlot.savefig("juliaData.png", dpi=600)
PyPlot.close()

PyPlot.imshow(grayscale, cmap="Greys_r", interpolation="none")
PyPlot.savefig("juliaTest.png", dpi=600)
PyPlot.close()


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

print("\n\n")
display(current_found)
print("\n")
display(possible_matrices)
print("\n\n")

PyPlot.scatter(collect(1:size(current_found)[1]), current_found)
PyPlot.show()
PyPlot.close()

# Calculating gamma by the paper definition
gamma1 = abs(current_found[2] - current_found[3]) +
         abs(current_found[2] - current_found[4]) +
         abs(current_found[2] - current_found[5]) +
         abs(current_found[3] - current_found[4]) +
         abs(current_found[3] - current_found[5]) +
         abs(current_found[4] - current_found[5])


gamma2 = abs(current_found[6] - current_found[7]) +
         abs(current_found[6] - current_found[8]) +
         abs(current_found[6] - current_found[9]) +
         abs(current_found[7] - current_found[8]) +
         abs(current_found[7] - current_found[9]) +
         abs(current_found[8] - current_found[9])

gamma3 = abs(current_found[10] - current_found[11])

gamma4 = abs(current_found[12] - current_found[13]) +
         abs(current_found[12] - current_found[14]) +
         abs(current_found[12] - current_found[15]) +
         abs(current_found[13] - current_found[14]) +
         abs(current_found[13] - current_found[15]) +
         abs(current_found[14] - current_found[15])


total_gamma = (gamma1 + gamma2 + gamma3 + gamma4) / sample_size

# multiplying total_gamma by 100, since it should be a percentage
#print("\n\n", "Total gamma is: ", total_gamma * 100, "%\n\n")
print("\n\n", "Total gamma is: ", total_gamma, "\n\n")



