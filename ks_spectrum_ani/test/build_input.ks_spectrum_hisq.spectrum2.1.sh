#! /bin/bash

# Script for building parameter input for ks_spectrum code
# to produce output in the style of spectrum2

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

cat <<EOF

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

number_of_pbp_masses ${nmasses}

max_cg_iterations 300
max_cg_restarts 5
npbp_reps 1
prec_pbp 1

EOF

for ((m=0; m<${nmasses}; m++)); do

cat <<EOF
mass ${mass[m]}
${naik_cmd[m]}
error_for_propagator ${error_for_propagator[m]}
rel_error_for_propagator 0
EOF

done

cat <<EOF
# Description of base sources

number_of_base_sources 1

# base source 0

corner_wall
subset full
t0 0
source_label C
forget_source

# Description of completed sources

number_of_modified_sources 1

# source 0

source 0

identity
op_label I
forget_source

# Description of propagators

number_of_sets 1

# Parameters common to all members of this set

max_cg_iterations ${max_cg_iterations}
max_cg_restarts 5
check yes
momentum_twist 0 0 0
time_bc antiperiodic
precision ${precision}

source 0

number_of_propagators ${nmasses}

EOF

for ((m=0; m<${nmasses}; m++)); do

cat <<EOF
# Propagators for mass ${m}

# propagator ${m}

mass ${mass[m]}
${naik_cmd[m]}
error_for_propagator ${error_for_propagator[m]}
rel_error_for_propagator 0

fresh_ksprop
forget_ksprop

EOF

done

cat <<EOF
# Definition of quarks

number_of_quarks ${nmasses}

EOF

for ((m=0; m<${nmasses}; m++)); do

cat <<EOF
# mass ${m}

propagator ${m}

identity
op_label d
forget_ksprop
EOF

done

cat <<EOF
# Description of mesons

number_of_mesons ${nmasses}

EOF

for ((m=0; m<${nmasses}; m++)); do

cat <<EOF

# pair ${m}

pair ${m} ${m}
spectrum_request meson

save_corr_fnal ${corrfilet}
r_offset 0 0 0 ${t0}

number_of_correlators 8

correlator PION_PS p000  1 * 1 pion5  0 0 0 E E E
correlator PION_SC p000  1 * 1 pion05 0 0 0 E E E           
correlator RHO_VT  p000  1 * 3 rhox   0 0 0 E E E           
correlator RHO_VT  p000  1 * 3 rhoy   0 0 0 E E E           
correlator RHO_VT  p000  1 * 3 rhoz   0 0 0 E E E           
correlator RHO_PV  p000  1 * 3 rhox0  0 0 0 E E E           
correlator RHO_PV  p000  1 * 3 rhoy0  0 0 0 E E E           
correlator RHO_PV  p000  1 * 3 rhoz0  0 0 0 E E E           

EOF

done

cat  <<EOF

# Description of baryons

number_of_baryons ${nmasses}

EOF

for ((m=0; m<${nmasses}; m++)); do

cat <<EOF
 
# mass ${m}

triplet ${m} ${m} ${m}
spectrum_request baryon

save_corr_fnal ${corrfilet}
r_offset 0 0 0 ${t0}

number_of_correlators 1

correlator NUCLEON  1 * 1 nucleon

EOF

done


done
