#!/bin/bash


# first no tether polymer
echo "first no tether"

/home/mv19141/poly/bin/beadrodpmer-bin -in initpolymers.in -var beads 12 -var x0 -3.0 -var y0 -3.0 -var z0 -3.0 -var xN 3.0 -var yN 3.0 -var zN 3.0 -var seed $RANDOM -var fname first_nt

# first double tether polymer

/home/mv19141/poly/bin/beadrodpmer-bin -in initpolymers.in -var beads 30 -var x0 -10.0 -var y0 -3.0 -var z0 -3.0 -var xN 3.0 -var yN 10.0 -var zN 3.0 -var seed $RANDOM -var fname first_dt

# second no tether polymer


/home/mv19141/poly/bin/beadrodpmer-bin -in initpolymers.in -var beads 20 -var x0 -10.0 -var y0 -20.0 -var z0 -6.0 -var xN -10.0 -var yN -10.0 -var zN -6.0 -var seed $RANDOM -var fname second_nt


# first single tether polymer


/home/mv19141/poly/bin/beadrodpmer-bin -in initpolymers.in -var beads 18 -var x0 15.0 -var y0 15.0 -var z0 15.0 -var xN 25.0 -var yN 25.0 -var zN 25.0 -var seed $RANDOM -var fname first_st



# second single tether polymer

/home/mv19141/poly/bin/beadrodpmer-bin -in initpolymers.in -var beads 30 -var x0 10.0 -var y0 -10.0 -var z0 -10.0 -var xN -10.0 -var yN 10.0 -var zN 10.0 -var seed $RANDOM -var fname second_st
