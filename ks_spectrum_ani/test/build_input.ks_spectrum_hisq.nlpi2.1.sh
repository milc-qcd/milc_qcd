#! /bin/bash

paramfile=$1

if [ $# -lt 1 ]
then
    echo "Usage $0 <paramfile>"
    exit 1
fi

source $paramfile

vol3=$[${nx}*${ny}*${nz}]

norm=`echo $vol3 | awk '{print 1./$1}'`

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

evenandodd_wall
subset full
t0 ${t0}
source_label q
forget_source

# Description of modified sources

number_of_modified_sources 2

# source 1

source 0

funnywall1
op_label f1
forget_source

# source 2

source 0

funnywall2
op_label f2
forget_source

EOF

######################################################################
# Definition of propagators for two sets

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

# Propagators for set 0

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

# Propagators for set 1

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

cat  <<EOF

# Parameters for set 2

max_cg_iterations ${max_cg_iterations}
max_cg_restarts 5
check yes
momentum_twist 0 0 0
time_bc antiperiodic

source 2

number_of_propagators ${nmasses}
EOF

# Propagators for set 2

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

# Quarks for even and odd wall source

for ((m=0; m<${nmasses}; m++)); do

cat  <<EOF

# mass ${m}

propagator ${m}

identity
op_label d
forget_ksprop

EOF

done

# Quarks with funnywall1 source

for ((m=0; m<${nmasses}; m++)); do

cat  <<EOF

# mass ${m} WP

propagator $[${nmasses}+${m}]

identity
op_label d
forget_ksprop

EOF

done

# Quarks with funnywall2 source

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

number_of_mesons $[2*${nmasses}]

EOF

for ((m=0; m<${nmasses}; m++)); do

n1=$[${nmasses}+${m}]
n2=$[2*${nmasses}+${m}]

cat  <<EOF

# Even and odd wall with funnywall1
# pair 0 (mass ${m} )

pair ${m} ${n1}
spectrum_request meson

save_corr_fnal ${corrfilet}
r_offset 0 0 0 ${t0}

number_of_correlators 6

# Normalization is 1/vol3
correlator PION_5  p000  1 * ${norm} pion5  0 0 0 E E E
correlator PION_i5 p000  1 * ${norm} pioni5 0 0 0 E E E           
correlator PION_i  p000  1 * ${norm} pioni  0 0 0 E E E
correlator PION_s  p000  1 * ${norm} pions  0 0 0 E E E           
correlator RHO_i   p000  1 * ${norm} rhoi   0 0 0 E E E
correlator RHO_s   p000  1 * ${norm} rhois  0 0 0 E E E           

# pair 1 mass ${m}

pair ${m} ${n2}
spectrum_request meson

save_corr_fnal ${corrfilet}
r_offset 0 0 0 ${t0}

number_of_correlators 6

# Normalization is 1/vol3
correlator PION_05 p000  1 * ${norm} pion05 0 0 0 E E E
correlator PION_ij p000  1 * ${norm} pionij 0 0 0 E E E           
correlator PION_i0 p000  1 * ${norm} pioni0 0 0 0 E E E
correlator PION_0  p000  1 * ${norm} pion0  0 0 0 E E E           
correlator RHO_i0  p000  1 * ${norm} rhoi0  0 0 0 E E E
correlator RHO_0   p000  1 * ${norm} rho0   0 0 0 E E E           

EOF

done

cat  <<EOF

# Description of baryons

number_of_baryons 0

EOF
done
