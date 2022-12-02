
#mpiexec -np 8 ~/liqliq_with_polymers/build/src/ll_w_p -in input.in -var solution_seed $RANDOM -var poly_seed1 $RANDOM -var volfrac 0.3 -var nucstrength 6.0 -var run 1 -var chi 3.5

mpiexec -np 8 ~/liqliq_with_polymers/build/src/ll_w_p -in restart.in -var solution_seed $RANDOM -var poly_seed1 $RANDOM -var volfrac 0.3 -var nucstrength 6.0 -var run 1 -var chi 3.0

#mpiexec -np 8  ~/liqliq_with_polymers/build/src/ll_w_p -in read.in -var solution_seed $RANDOM -var poly_seed1 $RANDOM -var volfrac 0.3 -var nucstrength 6.0 -var run 1 -var chi 3.0
