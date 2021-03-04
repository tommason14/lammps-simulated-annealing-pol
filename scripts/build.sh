#!/usr/bin/env bash
cwd="$(pwd)"
job="$(ls *job)"
for xyz in $(ls *xyz)
do
newdir="$(echo $xyz | sed 's/.xyz$//')"
echo "Processing $newdir" | tee -a debug.txt
[[ -d $newdir ]] || mkdir $newdir
data="$(echo $xyz | sed 's/.xyz$/.data/')"
poldata="$(echo $xyz | sed 's/.xyz$/-pol.data/')"
cp $xyz $newdir
cp $job $newdir
cd $newdir
python3 ../lmp_opls.py $xyz ../il.ff >> ../debug.txt 2>&1

#########################################################################################################
#  Remove pair coefficients from datafile, and prepare temporary input file for input_polariser to read #
#########################################################################################################

sed -n '/Pair Coeffs/,/^\s*[A-Z]/p' $data | grep '^[0-9]' | awk '{print "pair_coeff",$1,$1,"lj/cut/coul/long",$2,$3,$4,$5}' > tmp.in
# input_polariser looks for extra lines past the pair coeffs
cat tmp.in <(echo) > tmptmp.in && mv tmptmp.in tmp.in
# remove pair coeffs
cat <(sed '/Pair Coeffs/Q' $data) <(sed -ne '/Bond Coeffs/,$ p' $data) > tmp.data

polarizer -f ../cdrude.dff tmp.data $poldata >> ../debug.txt 2>&1
input_polariser tmp.in -d $poldata >> ../debug.txt 2>&1

####################################################
#  Expand boxsize to ± 30 Å, with walls at ± 20 Å  #
####################################################

sed -i 's/.*xlo.*/-30.000  30.000  xlo xhi/;s/.*ylo.*/-30.000  30.000  ylo yhi/;s/.*zlo.*/-30.000  30.000 zlo zhi/' $poldata >> ../debug.txt 2>&1

#############################
#  Now write an input file  #
#############################

atoms=$(sed -n '/Atoms/,/Bonds/p' $poldata |
grep -v 'DP$' |
grep '^\s*[0-9]' |
awk '{print $3}' |
sort -nu |
tr '\n' ' ' |
sed 's/ $//')

cores=$(sed -n '/Atoms/,/Bonds/p' $poldata |
grep 'DC$' |
awk '{print $3}' |
sort -nu |
tr '\n' ' ' |
sed 's/ $//')

drudes=$(sed -n '/Atoms/,/Bonds/p' $poldata |
grep 'DP$' |
awk '{print $3}' |
sort -nu |
tr '\n' ' ' |
sed 's/ $//')

shake=$(sed -n '/Bond Coeffs/,/Angle Coeffs/p' $poldata |
grep -v "DC-DP" |
grep "# H\|-H" |
awk '{print $1}' |
tr '\n' ' ' |
sed 's/ $//')

drudefix=$(sed -n '/Masses/,/^[A-Z]/p' $poldata |
grep "^\s*[0-9]" |
awk '{if ($NF=="DC"){print "C"} else if ($NF=="DP"){print "D"} else {print "N"}}' | 
tr '\n' ' ' |
sed 's/ $//')

elements=$(cat $poldata |
python3 -c "
import sys
import re
atomic_wt = ['H',  'Li', 'B',  'C',
             'N',  'O',  'F',  'Ne',
             'Na', 'Mg', 'Al', 'Si',
             'P',  'S',  'Cl', 'Ar',
             'K',  'Ca', 'Ti', 'Fe',
             'Zn', 'Se', 'Br', 'Kr',
             'Mo', 'Ru', 'Sn', 'Te',
             'I',  'Xe', 'D',  'M']
def atomic_symbol(name):
    if name[:2] in atomic_wt:
        return name[:2]
    elif name[0] in atomic_wt:
        return name[0]
    else:
        return name
types = []
found = False
for line in sys.stdin:
    if 'Masses' in line:
        found = True
        continue
    if len(types) > 0 and re.search('^\s*[A-Z]', line):
        break
    if found and len(line.split()) > 1:
        if line.strip().endswith('DC'):
            types.append(line.split()[-2])
        elif line.strip().endswith('DP'):
            types.append('D')
        else:
            types.append(line.split()[-1])

types = [atomic_symbol(name) for name in types]
print(' '.join(types))
")

cat << EOF > anneal.in
units           real
boundary        p p p
neighbor        2.0 bin
neigh_modify    every 1 delay 0 check yes

atom_style      full
bond_style      harmonic
angle_style     harmonic
dihedral_style  opls
improper_style  cvff
pair_style      hybrid/overlay lj/cut/coul/long 12.0 12.0 coul/long/cs 12.0 thole 2.6 12.0
kspace_style    pppm 1e-5
special_bonds   lj/coul 0.0 0.0 0.5

read_data $poldata extra/special/per/atom 3 
include pair-drude.lmp 

group ATOMS type $atoms
group CORES type $cores
group DRUDES type $drudes

fix DRUDE all drude $drudefix
fix SHAKE ATOMS shake 0.0001 20 0 b $shake
fix MOM all momentum 100 linear 1 1 1
comm_modify vel yes 

timestep        1

fix wall1 all   wall/lj93 xlo  -20.0 5.0 3.7 2.5 pbc yes
fix wall2 all   wall/lj93 xhi  20.0 5.0 3.7 2.5 pbc yes
fix wall3 all   wall/lj93 ylo  -20.0 5.0 3.7 2.5 pbc yes
fix wall4 all   wall/lj93 yhi  20.0 5.0 3.7 2.5 pbc yes
fix wall5 all   wall/lj93 zlo  -20.0 5.0 3.7 2.5 pbc yes
fix wall6 all   wall/lj93 zhi  20.0 5.0 3.7 2.5 pbc yes

restart         1000 run.restart1 run.restart2

variable TA equal 400
variable TD equal 1

##########################################
#  equilibrate for 1 ns, cool for 10 ns  #
##########################################

velocity ATOMS create \${TA} \${TA} dist gaussian
velocity DRUDES create \${TD} \${TD} dist gaussian

fix nvt all tgnvt/drude temp \${TA} \${TA} 100 \${TD} 20

thermo 100
thermo_style custom step f_nvt[1] f_nvt[2] f_nvt[3] press vol density pe ke etotal enthalpy evdwl ecoul epair 

dump d1 ATOMS  custom 10000 equilibrate.lmp element xu yu zu
dump_modify d1 element $elements
dump_modify d1 sort id
run             1000000
undump d1
unfix nvt

fix nvt all tgnvt/drude temp \${TA} 10 100 \${TD} 20

dump d1 ATOMS   custom 10000 cooling.lmp element xu yu zu
dump_modify d1  element $elements
dump_modify d1  sort id
run             10000000
EOF
echo "----" >> ../debug.txt
rm tmp.in tmp.data $xyz
cd "$cwd"
done
