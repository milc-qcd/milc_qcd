#! /bin/bash

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

number_of_base_sources 1

# base source 0

corner_wall
subset full
t0 ${t0}
source_label c
forget_source

# Description of modified sources

number_of_modified_sources 0

EOF

######################################################################
# Definition of propagators for one set

cat  <<EOF

# Description of propagators

number_of_sets 1

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

# Propagators for random wall source

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

######################################################################
# Definition of quarks

cat  <<EOF

number_of_quarks ${nmasses}

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

######################################################################
# Specification of Mesons

cat  <<EOF
# Description of mesons

number_of_mesons $[$[${nmasses}+1]*${nmasses}/2]

EOF

k=0

for ((m=0; m<${nmasses}; m++)); do
for ((n=${m}; n<${nmasses}; n++)); do

cat  <<EOF

# pair ${k} (masses ${m} ${n})

pair ${m} ${n}
spectrum_request meson

save_corr_fnal ${corrfilet}
r_offset 0 0 0 ${t0}

number_of_correlators 8

correlator PION_5  p000  1 * 1 pion5  0 0 0 E E E
correlator PION_05 p000  1 * 1 pion05 0 0 0 E E E           
correlator RHO_i   p000  1 * 3 rhox   0 0 0 E E E           
correlator RHO_i   p000  1 * 3 rhoy   0 0 0 E E E           
correlator RHO_i   p000  1 * 3 rhoz   0 0 0 E E E           
correlator RHO_i0  p000  1 * 3 rhox0  0 0 0 E E E           
correlator RHO_i0  p000  1 * 3 rhoy0  0 0 0 E E E           
correlator RHO_i0  p000  1 * 3 rhoz0  0 0 0 E E E           

EOF

k=$[${k}+1]

done
done

######################################################################
# Specification of baryons

cat  <<EOF
# Description of baryons

number_of_baryons $[$[${nmasses}+2]*$[${nmasses}+1]*${nmasses}/6]

EOF

k=0

for ((m0=0; m0<${nmasses}; m0++)); do
for ((m1=${m0}; m1<${nmasses}; m1++)); do
for ((m2=${m1}; m2<${nmasses}; m2++)); do

cat  <<EOF

# triplet ${k} (masses ${m0} ${m1} ${m2})

triplet ${m0} ${m1} ${m2}
spectrum_request baryon

save_corr_fnal ${corrfilet}
r_offset 0 0 0 ${t0}

number_of_correlators 1

correlator NUCLEON  1 * 1 nucleon

EOF

k=$[${k}+1]

done
done
done

reload_gauge_cmd="continue"

done # sources
