#!/usr/bin/tcl

proc Fasta_Count {argv} {

	set f_in1 [open [lindex $argv 0] "r"]
	set q_len [lindex $argv 1]
	set f_out [open [lindex $argv 2] "w"]
	set line_limit [lindex $argv 3]

	global q_array_all
	global q_array_a
	global q_array_t
	global q_array_g
	global q_array_c
	global q_array_n
	global q_array_x
	global q_array_z

	### HEADER OUT FILE ###

	puts $f_out "\#\tA_abs\tT_abs\tG_abs\tC_abs\tN_abs\tX_abs\tZ_abs\tALL_abs\t\*\*\*\tALL_sum\t\*\*\*\tA_fr\tT_fr\tG_fr\tC_fr\tATGC\tAT\tGC\tNXZ\tALL_fr"

	####### SET ZERO FOR DATAPOINTS ######

	set q 0
	while { $q < $q_len } {
	set q_array_all($q) 0
	set q_array_a($q)   0
	set q_array_t($q)   0
	set q_array_g($q)   0
	set q_array_c($q)   0
	set q_array_n($q)   0
	set q_array_x($q)   0
	set q_array_z($q)   0
	incr q
    }

	####### READ TABLE INTO MEMORY #######
	set l 1
	while { [gets $f_in1 current_line] >= 0 && $l <= $line_limit} {
	set current_line [string toupper $current_line]
	set current_data [split   $current_line ""]
	set data_length  [llength $current_data]
	set k 0
	if { [lindex $current_data 0] != ">" && $data_length >= 24 } {
		while { $k < $data_length } {
			set current_q   [lindex  $current_data $k]
			incr q_array_all($k)
			if { $current_q == "A" } {
				incr q_array_a($k)
			}
			if { $current_q == "T" } {
				incr q_array_t($k)
			}
			if { $current_q == "G" } {
				incr q_array_g($k)
			}
			if { $current_q == "C" } {
				incr q_array_c($k)
			}
			if { $current_q == "N" } {
				incr q_array_n($k)
			}
			if { $current_q == "X" } {
				incr q_array_x($k)
			}
			if { $current_q == "*" } {
				incr q_array_z($k)
			}
			incr k 
			}
			incr l
			set l_mod [expr fmod($l,1000)]
			if { $l_mod == 0 } {
				puts " ...  $l  ... "
				puts $current_data
			}
		}
	}
	close $f_in1

	set n 0
	while { $n < $q_len } {

	set a_fract [expr int(round($q_array_a($n)*1000.00/$q_array_all($n)))]
	set a_fract [expr $a_fract/10.0]
	set t_fract [expr int(round($q_array_t($n)*1000.00/$q_array_all($n)))]
	set t_fract [expr $t_fract/10.0]
	set g_fract [expr int(round($q_array_g($n)*1000.00/$q_array_all($n)))]
	set g_fract [expr $g_fract/10.0]
	set c_fract [expr int(round($q_array_c($n)*1000.00/$q_array_all($n)))]
	set c_fract [expr $c_fract/10.0]
	set n_fract [expr int(round($q_array_n($n)*1000.00/$q_array_all($n)))]
	set n_fract [expr $n_fract/10.0]
	set x_fract [expr int(round($q_array_x($n)*1000.00/$q_array_all($n)))]
	set x_fract [expr $x_fract/10.0]
	set z_fract [expr int(round($q_array_z($n)*1000.00/$q_array_all($n)))]
	set z_fract [expr $z_fract/10.0]

	set q_sum [expr $q_array_a($n) + $q_array_t($n) + $q_array_g($n) + $q_array_c($n) + $q_array_n($n) + $q_array_x($n) + $q_array_z($n)]

	set atgc_sum [expr $q_array_a($n) + $q_array_t($n) + $q_array_g($n) + $q_array_c($n)]
	set nxz_sum  [expr $q_array_n($n) + $q_array_x($n) + $q_array_z($n)]

	set at_sum [expr $q_array_a($n) + $q_array_t($n)]
	set gc_sum [expr $q_array_g($n) + $q_array_c($n)]

	set all_fract [expr int(round($q_sum*1000.00/$q_array_all($n)))]
	set all_fract [expr $all_fract/10.0]

	set at_fract [expr int(round($at_sum*1000.00/$q_array_all($n)))]
	set at_fract [expr $at_fract/10.0]
	set gc_fract [expr int(round($gc_sum*1000.00/$q_array_all($n)))]
	set gc_fract [expr $gc_fract/10.0]

	set atgc_fract [expr int(round($atgc_sum*1000.00/$q_array_all($n)))]
	set atgc_fract [expr $atgc_fract/10.0]
	set nxz_fract [expr int(round($nxz_sum*1000.00/$q_array_all($n)))]
	set nxz_fract [expr $nxz_fract/10.0]

	puts $f_out "$n\t$q_array_a($n)\t$q_array_t($n)\t$q_array_g($n)\t$q_array_c($n)\t$q_array_n($n)\t$q_array_x($n)\t$q_array_z($n)\t$q_array_all($n)\t\*\*\*\t$q_sum\t\*\*\*\t$a_fract\t$t_fract\t$g_fract\t$c_fract\t$atgc_fract\t$at_fract\t$gc_fract\t$nxz_fract\t$all_fract"
	puts "$n\t$q_array_a($n)\t$q_array_t($n)\t$q_array_g($n)\t$q_array_c($n)\t$q_array_n($n)\t$q_array_x($n)\t$q_array_z($n)\t$q_array_all($n)\t\*\*\*\t$q_sum\t\*\*\*\t$a_fract\t$t_fract\t$g_fract\t$c_fract\t$atgc_fract\t$at_fract\t$gc_fract\t$nxz_fract\t$all_fract"
	incr n
	}

	close $f_out
	puts ""
	puts "DONE"
}

if {$argc != 4} {
	puts "Program usage:"
	puts "Sequence_File, Total_Length, output_file, line_limit"
} else {
	Fasta_Count $argv
}

