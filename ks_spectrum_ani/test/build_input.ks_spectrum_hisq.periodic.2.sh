#! /bin/bash

paramfile=$1

if [ $# -lt 1 ]
then
    echo "Usage $0 <paramfile>"
    exit 1
fi

source $paramfile

vol3=$[${nx}*${ny}*${nz}]

ppnorm=`echo $vol3 | awk '{print 1./($1**2)}'`
pwnorm=`echo $vol3 | awk '{print 1./($1**3)}'`

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

point
subset full
origin 0 0 0 0
source_label d
forget_source

# Description of completed sources

number_of_modified_sources 0
EOF

######################################################################
# Definition of propagators for two sets

cat  <<EOF

# Description of propagators

number_of_sets 1

# Parameters for set 0

max_cg_iterations ${max_cg_iterations}
max_cg_restarts 5
check yes
momentum_twist 0 0 0
time_bc periodic
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
save_serial_fm_ksprop ks.test.l4448_m${mass[m]}_t0_x0_y0_z0

EOF

done

######################################################################
# Definition of quarks

cat  <<EOF

number_of_quarks $[2*${nmasses}]

EOF

# Quarks with point source and point sink

for ((m=0; m<${nmasses}; m++)); do

cat  <<EOF

# mass ${m} PP

propagator ${m}

identity
op_label d
forget_ksprop

EOF

done

# Quarks with point source and wall sink

for ((m=0; m<${nmasses}; m++)); do

cat  <<EOF

# mass ${m} PW

propagator ${m}

evenandodd_wall
op_label EO
forget_ksprop

EOF

done

######################################################################
# Specification of Mesons

cat  <<EOF
# Description of mesons

number_of_mesons $[$[${nmasses}+1]*${nmasses}]

EOF

k=0

for ((m=0; m<${nmasses}; m++)); do
for ((n=${m}; n<${nmasses}; n++)); do

mpp=${m}
mpw=$[${nmasses}+${m}]

npp=${n}
npw=$[${nmasses}+${n}]

cat  <<EOF

# pair ${k} (masses ${m} ${n} PP)

pair ${mpp} ${npp}
spectrum_request meson

save_corr_fnal ${corrfilet}
r_offset 0 0 0 ${t0}

number_of_correlators 1

# Normalization is 1/(vol3)^2

correlator POINT_KAON_5 p000  1 * ${ppnorm} pion5  0 0 0 E E E

# pair $[${k}+1] (masses ${m} ${n} PW)

pair ${mpw} ${npw}
spectrum_request meson

save_corr_fnal ${corrfilet}
r_offset 0 0 0 ${t0}

number_of_correlators 1

# Normalization is 1/(vol3)^3

correlator WALL_KAON_5 p000  1 * ${pwnorm} pion5  0 0 0 E E E

EOF

k=$[${k}+2]

done
done

cat  <<EOF

# Description of baryons

number_of_baryons 0

EOF

reload_gauge_cmd="continue"

done # sources
