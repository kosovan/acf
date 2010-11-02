# Test script for simulating a hydrogel network

# polymer parameters
set MPC 20;
set bond_length 0.95;
set type_node 0;
set type_mon 1;

# diffusing labels
set n_labels 50;
set type_label 2;
set q_label 0;

# lattice size = box length
set lat_a [expr 4.0/sqrt(3.0)*($MPC+1.0)*$bond_length];
set bl [expr 1*$lat_a];
setmd box_l $bl $bl $bl;
puts "box_length: [setmd box_l]";

# integration parameters
set time_step 0.0125;
setmd time_step $time_step;
setmd skin 0.4;
set temperature 1.0;
set gamma 2.0;
thermostat langevin $temperature $gamma;
# overall simulation time 
set tot_time 20000; # total simulation time in reduced units
set round_time 1; # time elapsed between two integration rounds
set int_times [expr int($tot_time/$round_time)];
set int_steps [expr int($round_time/$time_step)];;
# warmup parameters
set ljcap 10;
set warm_times 10;
set warm_steps 10;
set mindist 0.95;
set act_mindist 0.0;

# Interaction parameters

#FENE
set fene_k 7.0;
set fene_cut 2.0;

# WCA
set wca_epsilon 1.0;
set wca_sigma 1.0;
set wca_cut [expr pow(2,1.0/6.0)];
set wca_shift [expr pow([expr $wca_sigma/$wca_cut],6)-pow([expr $wca_sigma/$wca_cut],12)];
set wca_off 0.0;


# directory setup
set name "freehg$MPC\_L$n_labels\_g$gamma\_test";
set wd [pwd];
set data_d "/home/kosovan/data/$name";
file mkdir $data_d;
cd $data_d;
set vtf 1;
set vtf_store 1;
set vtf_file_name "$name.vtf";
set vtf_file [open "$vtf_file_name" "w"];

# observable output
set store_time 0;
set cfg_index 0;
# Create vmd script
set vmd_file [open "vmd_movie.script" "w"]
puts $vmd_file "source /home/kosovan/Espresso/espresso-2.1.5/scripts/vmd_plg.tcl"
puts $vmd_file "#loadseries $name.vmd"
puts $vmd_file "mol load vtf $vtf_file_name"
puts $vmd_file "mol modstyle 0 0 CPK 2.500000 0.000000 8.000000 5.000000"
puts $vmd_file "mol modcolor 0 0 Name"
puts $vmd_file "nearclip set 0.01"
puts $vmd_file "logfile vmd.log"
puts $vmd_file "#scale by 1.7"
puts $vmd_file "axes location off"
puts $vmd_file "logfile off"
puts $vmd_file "pbc set [setmd box_l]"
puts $vmd_file "pbc box_draw"
#puts $vmd_file "animate forward"
close $vmd_file

# set the interactions
inter 0 FENE $fene_k $fene_cut;

# set the wca interaction 
inter $type_label $type_label lennard-jones $wca_epsilon $wca_sigma $wca_cut $wca_shift $wca_off;

# setup the polymer lattice
set n_part_types 0;
set id_node_min 0;
set ic_node_max 7;
set id_node_min [expr [setmd n_part]]; 
set id_pol_min [expr [setmd n_part]+8]; 
# no polymer for free diffusion
#diamond $lat_a $bond_length $MPC;
set id_pol_max [expr [setmd n_part]-1];

# add the lables
set id_ci_min [setmd n_part];
counterions $n_labels charge $q_label type $type_label;
set id_ci_max [expr [setmd n_part]-1];


# write down the configuration to check that everything is correctly set up
set config 0;

# file to write down particle configurations
set test_file_name "msd.test";
set test_file [open $test_file_name "w"];

# warmup integration 
for {set i 0} {$i<$warm_times} {incr i} {
  if {![expr $i%10]} { puts -nonewline "Warmup run $i, act_mindist: $act_mindist mindist: $mindist\n"; flush stdout; }
  inter ljforcecap $ljcap;
  integrate $warm_steps;
  set act_mindist [analyze mindist];
  incr ljcap 10;
}
inter ljforcecap 0;
puts "\nWarmup finished"; flush stdout;

puts "\nTune cells:\n[tune_cells]"; flush stdout;

# Write down the state just before the main integration
polyBlockWriteAll "$name.warm";

# Main integration 
puts "Start main integration";
puts "interactions: [inter]"; flush stdout;
set start_time [clock seconds];
puts "start at $start_time (HW clock)";
set time [setmd time];
for {set i 0} {$i<$int_times} {incr i} {
  #puts "integration round $i"; flush stdout;
  integrate $int_steps;
  set time [setmd time];
  if {![expr $i%100]} { puts -nonewline "Main run $i, time=$time\n"; flush stdout; }
  # we set some store_time by default but it shoudl be checked afterwards!
  #if {$time>$store_time} { analyze append; } 
  if {$time>$store_time} { 
	  analyze append; 
	  # write down the configuration for reading 
	  for {set pid $id_ci_min} {$pid<=$id_ci_max} {incr pid} { puts -nonewline $test_file "[part $pid print pos] "; }
	  puts $test_file " ";
  } 
  # write down the configuration for later analysis
  set fn [open "$name.[format %04d $cfg_index]" "w"]; 
  blockfile $fn write variable  {time};
  blockfile $fn write particles {id type pos v};
  close $fn;
  incr cfg_index;
  # write down a checkpoint
  polyBlockWriteAll "$name.end" all -;
  #if {$vtf && ![expr $i%$vtf_store]} {writevcf $vtf_file short folded pids all; incr config 1;}
}
close $test_file;
# write down the last configuration and the last checkpoint
polyBlockWriteAll "$name.end" all -;
set fn [open "$name.[format %04d $cfg_index]" "w"]; 
blockfile $fn write particles {id type pos v};
close $fn;
close $vtf_file;
puts "\nIntegration finished"; 
set end_sim_time [clock seconds];

puts "Analyzing van Hove ACF"; flush stdout;
set vh_type $type_label;
set vh_rmin 0;
set vh_rmax 25;
set vh_rbins 50;
set dt [expr $time_step*$int_steps];
set vh [analyze vanhove $vh_type $vh_rmin $vh_rmax $vh_rbins ];
set vh_file [open "$name.vh" "w"];
puts $vh_file $vh;
close $vh_file;

# write down vh data in the form readable for Gnuplot
set msdfile [open "$name.msd" "w"];
set msdt [lindex $vh 0];
set msd [lindex $msdt 1];
set maxt [llength $msd];
for {set t 0} {$t<$maxt} {incr t} {	
	puts $msdfile "[expr ($t+1)*$dt] [lindex $msd $t]" 
}
close $msdfile;

set vh_gp_file [open "$name.vh.gp" "w"];
set rfac [expr 1.0*$vh_rmax/$vh_rbins];
set Grt [lindex [lindex $vh 1] 1];
set maxt [llength $Grt];
for {set t 0} {$t<$maxt} {incr t} {
	set Gr [lindex $Grt $t];
	for {set r 0} {$r<$vh_rbins} {incr r} {
		puts $vh_gp_file "[format %6e [expr $r*$rfac]] [format %6e [expr $t*$dt]] [format %6e [lindex $Gr $r]]"
	}
	puts $vh_gp_file " ";
}
close $vh_gp_file;

puts "dt: $dt"; flush stdout;
set end_time [clock seconds];
set duration [expr $end_time-$start_time];
set sim_duration [expr $end_sim_time-$start_time];
set vh_duration [expr $end_time-$end_sim_time];
puts "end at $start_time, total duration: $duration (s) = [expr $duration/60] (min);";
puts "van Hove duration: $vh_duration (s) = [expr $vh_duration/60] (min);\nsimulation duration: $sim_duration (s) = [expr $sim_duration/60] (min)";
puts "MPC: $MPC, total particles: [setmd n_part]";

puts "Happily finished :-)"; exit;

# TODO

# integration at constant pressure
# electrostatics and counterions
