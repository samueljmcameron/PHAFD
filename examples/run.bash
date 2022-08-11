mkdir -p xvtkfiles

mpiexec -np 4 ../a.out -in input.in -nuc nucleationsites.nuc -var label spread -var seed 2948012 -var solution_seed 84991 -var poly_seed1 418029 -var poly_seed2 668591 -var poly_seed3 74910
