#! /bin/bash

# Script for building parameter input for ks_spectrum code
# to produce output in the style of nl_spectrum 

paramfile=$1

if [ $# -lt 1 ]
then
    echo "Usage $0 <paramfile>"
    exit 1
fi

source $paramfile

vol3=$[${nx}*${ny}*${nz}]

reload_gauge_cmd="reload_serial ${inlat}"

for ((i=0; i<${nmasses}; i++)); do
case $action in
hisq)
  naik_cmd[$i]="naik_term_epsilon ${naik_term_epsilon[i]}"
;;
asqtad)
  naik_cmd[$i]=""
;;
esac
done

cat <<EOF
prompt 0
nx ${nx}
ny ${ny}
nz ${nz}
nt ${nt}
iseed ${iseed}
job_id ${jobid}
EOF

# Iterate over source time slices
for ((i=0; i<${n_sources}; i++)); do
t0=$[${source_start}+${i}*${source_inc}]
corrfilet=${corrfile}_t${t0}.test-out

cat  <<EOF 

######################################################################
# source time ${t0}
######################################################################

# Gauge field description

${reload_gauge_cmd}
u0 ${u0}
coulomb_gauge_fix
forget
staple_weight 0
ape_iter 0
coordinate_origin 0 0 0 0

# Chiral condensate and related measurements

number_of_pbp_masses 0

# Description of base sources

number_of_base_sources 3

# base source 0

even_wall
subset full
t0 ${t0}
source_label E
forget_source

# base source 1

evenandodd_wall
subset full
t0 ${t0}
source_label q
forget_source

# base source 2

evenminusodd_wall
subset full
t0 ${t0}
source_label o
forget_source

# Description of modified sources

number_of_modified_sources 0

EOF

######################################################################
# Definition of propagators

cat  <<EOF

# Description of propagators

number_of_sets 3

# Parameters for set 0

max_cg_iterations ${max_cg_iterations}
max_cg_restarts 5
check yes
momentum_twist 0 0 0
time_bc antiperiodic
precision ${precision}

source 0

number_of_propagators ${nmasses}
EOF

# Propagators for even wall source

for ((m=0; m<${nmasses}; m++)); do

cat  <<EOF

# propagator ${m}

mass ${mass[$m]}
${naik_cmd[$m]}
error_for_propagator ${error_for_propagator[$m]}
rel_error_for_propagator 0

fresh_ksprop
forget_ksprop

EOF

done

cat  <<EOF

# Parameters for set 1

max_cg_iterations ${max_cg_iterations}
max_cg_restarts 5
check yes
momentum_twist 0 0 0
time_bc antiperiodic
precision ${precision}

source 1

number_of_propagators ${nmasses}
EOF

# Propagators for evenandodd_wall source

for ((m=0; m<${nmasses}; m++)); do

cat  <<EOF

# propagator $[${nmasses}+${m}]

mass ${mass[$m]}
${naik_cmd[$m]}
error_for_propagator ${error_for_propagator[$m]}
rel_error_for_propagator 0

fresh_ksprop
forget_ksprop

EOF

done

# Parameters for set 2

cat  <<EOF

max_cg_iterations ${max_cg_iterations}
max_cg_restarts 5
check yes
momentum_twist 0 0 0
time_bc antiperiodic
precision ${precision}

source 2

number_of_propagators ${nmasses}
EOF

# Propagators for evenminusodd_wall source

for ((m=0; m<${nmasses}; m++)); do

cat  <<EOF

# propagator $[2*${nmasses}+${m}]

mass ${mass[$m]}
${naik_cmd[$m]}
error_for_propagator ${error_for_propagator[$m]}
rel_error_for_propagator 0

fresh_ksprop
forget_ksprop

EOF

done

######################################################################
# Definition of quarks

cat  <<EOF

number_of_quarks $[3*${nmasses}]

EOF

for ((m=0; m<${nmasses}; m++)); do

cat  <<EOF

# mass ${m}

propagator ${m}

identity
op_label d
forget_ksprop

EOF

done

for ((m=0; m<${nmasses}; m++)); do

cat  <<EOF

# mass ${m}

propagator $[${nmasses}+${m}]

identity
op_label d
forget_ksprop

EOF

done

for ((m=0; m<${nmasses}; m++)); do

cat  <<EOF

# mass ${m}

propagator $[2*${nmasses}+${m}]

identity
op_label d
forget_ksprop

EOF

done

######################################################################
# Specification of Mesons

cat  <<EOF
# Description of mesons

number_of_mesons $[4*${nmasses}]

EOF

k=0

for ((m=0; m<${nmasses}; m++)); do

cat  <<EOF

# pair ${k} mass ${m} even wall / even wall

pair ${m} ${m}
spectrum_request meson

save_corr_fnal ${corrfilet}
r_offset 0 0 0 ${t0}

number_of_correlators 2

correlator PION_PS p000  1 / 16 pion5  0 0 0 E E E
correlator PION_SC p000  1 / 16 pion05 0 0 0 E E E           

EOF

k=$[${k}+1]

done

for ((m=0; m<${nmasses}; m++)); do

cat  <<EOF

# pair ${k} mass ${m} evenandodd wall / evenandodd wall

pair $[${nmasses}+${m}] $[${nmasses}+${m}]
spectrum_request meson

save_corr_fnal ${corrfilet}
r_offset 0 0 0 ${t0}

number_of_correlators 1

correlator PION_PS_a p000  1 / 16 pion5  0 0 0 E E E

EOF

k=$[${k}+1]

done

for ((m=0; m<${nmasses}; m++)); do

cat  <<EOF

# pair ${k} mass ${m} evenminusodd wall / evenminusodd wall

pair $[2*${nmasses}+${m}] $[2*${nmasses}+${m}]
spectrum_request meson

save_corr_fnal ${corrfilet}
r_offset 0 0 0 ${t0}

number_of_correlators 1

correlator PION_PS_b p000  1 / 16 pion5  0 0 0 E E E

EOF

k=$[${k}+1]

done

for ((m=0; m<${nmasses}; m++)); do

cat  <<EOF

# pair ${k} mass ${m} evenandodd wall / evenminusodd wall

pair $[${nmasses}+${m}] $[2*${nmasses}+${m}]
spectrum_request meson

save_corr_fnal ${corrfilet}
r_offset 0 0 0 ${t0}

number_of_correlators 1

correlator PION_SC p000  1 / 16 pion05  0 0 0 E E E

EOF

k=$[${k}+1]

done


######################################################################
# Specification of baryons

cat  <<EOF
# Description of baryons

number_of_baryons $[2*${nmasses}]

EOF

k=0

for ((m=0; m<${nmasses}; m++)); do

cat  <<EOF

# triplet ${k} mass ${m} even wall

triplet ${m} ${m} ${m}
spectrum_request baryon

save_corr_fnal ${corrfilet}
r_offset 0 0 0 ${t0}

number_of_correlators 1

correlator NUCLEON  1 / 64 nucleon

EOF

k=$[${k}+1]

done

for ((m=0; m<${nmasses}; m++)); do

cat  <<EOF

# triplet ${k} mass ${m} evenandodd wall

triplet $[${nmasses}+${m}] $[${nmasses}+${m}] $[${nmasses}+${m}]
spectrum_request baryon

save_corr_fnal ${corrfilet}
r_offset 0 0 0 ${t0}

number_of_correlators 2

correlator NUCLEON  1 / 64 nucleon
correlator DELTA    1 / 64 delta

EOF

k=$[${k}+1]

done

reload_gauge_cmd="continue"

done # sources
