#! /bin/bash

paramfile=$1

if [ $# -lt 1 ]
then
    echo "Usage $0 <paramfile>"
    exit 1
fi

source $paramfile

vol3=$[${nx}*${ny}*${nz}]

ppnorm=`echo $vol3 | awk '{print 1./(3*$1**2)}'`
pwnorm=`echo $vol3 | awk '{print 1./(3*$1**3)}'`
wpnorm=`echo $vol3 | awk '{print 4./(3*$1**2)}'`
wwnorm=`echo $vol3 | awk '{print 4./(3*$1**3)}'`

reload_gauge_cmd="reload_serial ${inlat}"

for ((i=0; i<${fpi_nmasses}; i++)); do
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

number_of_base_sources 2

# base source 0

random_color_wall
subset full
t0 ${t0}
ncolor ${nrand_source}
momentum 0 0 0
source_label r
forget_source

# base source 1

evenandodd_wall
subset full
t0 ${t0}
source_label q
forget_source

# Description of completed sources

number_of_modified_sources 0

EOF

######################################################################
# Definition of propagators for two sets

cat  <<EOF

# Description of propagators

number_of_sets 2

# Parameters for set 0

max_cg_iterations ${max_cg_iterations}
max_cg_restarts 5
check yes
momentum_twist 0 0 0
time_bc antiperiodic
precision ${precision}

source 0

number_of_propagators ${fpi_nmasses}
EOF

# Propagators for random wall source

for ((m=0; m<${fpi_nmasses}; m++)); do

cat  <<EOF

# propagator ${m}

mass ${fpi_mass[$m]}
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

number_of_propagators ${fpi_nmasses}
EOF

# Propagators for evenodd wall source

for ((m=0; m<${fpi_nmasses}; m++)); do

cat  <<EOF

# propagator $[${fpi_nmasses}+${m}]

mass ${fpi_mass[$m]}
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

number_of_quarks $[4*${fpi_nmasses}]

EOF

# Quarks with point source and point sink

for ((m=0; m<${fpi_nmasses}; m++)); do

cat  <<EOF

# mass ${m} PP

propagator ${m}

identity
op_label d
forget_ksprop

EOF

done

# Quarks with wall source and point sink

for ((m=0; m<${fpi_nmasses}; m++)); do

cat  <<EOF

# mass ${m} WP

propagator $[${fpi_nmasses}+${m}]

identity
op_label d
forget_ksprop

EOF

done

# Quarks with point source and wall sink

for ((m=0; m<${fpi_nmasses}; m++)); do

cat  <<EOF

# mass ${m} PW

propagator ${m}

evenandodd_wall
op_label EO
forget_ksprop

EOF

done

# Quarks with wall source and wall sink

for ((m=0; m<${fpi_nmasses}; m++)); do

cat  <<EOF

# mass ${m} WW

propagator $[${fpi_nmasses}+${m}]

evenandodd_wall
op_label EO
forget_ksprop

EOF

done

######################################################################
# Specification of Mesons

cat  <<EOF
# Description of mesons

number_of_mesons $[$[${fpi_nmasses}+1]*${fpi_nmasses}*2]

EOF

k=0

for ((m=0; m<${fpi_nmasses}; m++)); do
for ((n=${m}; n<${fpi_nmasses}; n++)); do

mpp=${m}
mwp=$[${fpi_nmasses}+${m}]
mpw=$[2*${fpi_nmasses}+${m}]
mww=$[3*${fpi_nmasses}+${m}]

npp=${n}
nwp=$[${fpi_nmasses}+${n}]
npw=$[2*${fpi_nmasses}+${n}]
nww=$[3*${fpi_nmasses}+${n}]

cat  <<EOF

# pair ${k} (masses ${m} ${n} PP)

pair ${mpp} ${npp}
spectrum_request meson

save_corr_fnal ${corrfilet}
r_offset 0 0 0 ${t0}

number_of_correlators 1

# Normalization is 1/[3*(vol3)^2]

correlator POINT_KAON_5 p000  1 * ${ppnorm} pion5  0 0 0 E E E

# pair $[${k}+1] (masses ${m} ${n} PW)

pair ${mpw} ${npw}
spectrum_request meson

save_corr_fnal ${corrfilet}
r_offset 0 0 0 ${t0}

number_of_correlators 1

# Normalization is 1/[3*(vol3)^3]

correlator WALL_KAON_5 p000  1 * ${pwnorm} pion5  0 0 0 E E E

# pair $[${k}+2] (masses ${m} ${n} WP)

pair ${mwp} ${nwp}
spectrum_request meson

save_corr_fnal ${corrfilet}
r_offset 0 0 0 ${t0}

number_of_correlators 1

# Normalization is 4/[3*(vol3)^2]

correlator POINT_KAON_5 p000  1 * ${wpnorm} pion5  0 0 0 E E E

# pair $[${k}+3] (masses ${m} ${n} WW)

pair ${mww} ${nww}
spectrum_request meson

save_corr_fnal ${corrfilet}
r_offset 0 0 0 ${t0}

number_of_correlators 1

# Normalization is 4/[3*(vol3)^3]

correlator WALL_KAON_5 p000  1 * ${wwnorm} pion5  0 0 0 E E E

EOF

k=$[${k}+4]

done
done

cat  <<EOF

# Description of baryons

number_of_baryons 0

EOF

reload_gauge_cmd="continue"

done # sources
