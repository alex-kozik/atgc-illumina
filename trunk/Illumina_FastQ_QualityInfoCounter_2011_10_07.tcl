#!/usr/bin/tcl

proc FastQ_Count {argv} {

	set phred_diff 64
	set f_in1 [open [lindex $argv 0] "r"]
	set q_len [lindex $argv 1]
	set f_out [open [lindex $argv 2] "w"]
	set first_Q_line [lindex $argv 3]
	set mod_line_number [lindex $argv 4]
	set line_limit [lindex $argv 5]
	set q_line_spacing [expr $mod_line_number - 4]

	global q_array_all
	global q_array_s
	global q_array_a
	global q_array_b
	global q_array_c
	global q_array_d
	global q_array_f
	global q_array_x

	set cutoff_A 33
	set cutoff_B 25
	set cutoff_C 20
	set cutoff_D 15
	set cutoff_F 10
	set cutoff_X 2
	
	### HEADER OUT FILE ###

	puts $f_out "\#\tS_abs\tA_abs\tB_abs\tC_abs\tD_abs\tF_abs\tX_abs\tALL_abs\t\*\*\*\tSABCD\tALL_sum\t\*\*\*\tS_33\tA_33\tB_25\tC_20\tD_15\tF_10\tX_02\tALL_fr\t\*\*\*\tS_D\tF_X\tALL"

	####### SET ZERO FOR DATAPOINTS ######

	set q 0
	while { $q < $q_len } {
	set q_array_all($q) 0
	set q_array_s($q)   0
	set q_array_a($q)   0
	set q_array_b($q)   0
	set q_array_c($q)   0
	set q_array_d($q)   0
	set q_array_f($q)   0
	set q_array_x($q)   0
	incr q
	}

	####### READ TABLE INTO MEMORY #######
	set l 1
	set q_lines 0
	while { [gets $f_in1 current_line] >= 0 && $l <= $line_limit } {
		### set current_line [string toupper $current_line]
		set current_data [split   $current_line ""]
		set data_length  [llength $current_data]
		set k 0
		set q_line_counter [expr $l + $q_line_spacing]
		set q_line_mod [expr fmod($q_line_counter,$mod_line_number)]
		### puts "$l\t$current_data"
		if { $q_line_mod == 0 && $data_length >= 24 } {
			### puts "$l\t$current_data"
			while { $k < $data_length } {
				set current_q_ascii [lindex $current_data $k]
				scan [string index $current_q_ascii 0] %c current_q
				set current_q [expr $current_q - $phred_diff]
				### puts "$current_q_ascii\t$current_q"
				incr q_array_all($k)
				if { $current_q <= $cutoff_X } {
					incr q_array_x($k)
				}
				if { $current_q >  $cutoff_X && $current_q < $cutoff_F } {
					incr q_array_f($k)
				}
				if { $current_q >= $cutoff_F && $current_q < $cutoff_D } {
					incr q_array_d($k)
				}
				if { $current_q >= $cutoff_D && $current_q < $cutoff_C } {
					incr q_array_c($k)
				}
				if { $current_q >= $cutoff_C && $current_q < $cutoff_B } {
					incr q_array_b($k)
				}
				if { $current_q >= $cutoff_B && $current_q < $cutoff_A } {
					incr q_array_a($k)
				}
				if { $current_q >= $cutoff_A } {
					incr q_array_s($k)
				}
				incr k 
			}
			incr q_lines
			set l_mod [expr fmod($q_lines,1000)]
			if { $l_mod == 0 } {
				puts " ...  $l  ... "
				puts " TEST: $current_line "
			}
		}
		incr l
	}
	
	close $f_in1

	set n 0
	while { $n < $q_len } {

	set sabcd_sum [expr $q_array_s($n) + $q_array_a($n) + $q_array_b($n) + $q_array_c($n) + $q_array_d($n)]
	set fx_sum   [expr $q_array_f($n) + $q_array_x($n)]

	### FRACTION CALCULATION ###
	### puts "$q_array_a($n)\t$q_array_all($n)"
	set s_fract [expr int(round($q_array_s($n)*1000.00/$q_array_all($n)))]
	set s_fract [expr $s_fract/10.0]
	set a_fract [expr int(round($q_array_a($n)*1000.00/$q_array_all($n)))]
	set a_fract [expr $a_fract/10.0]
	set b_fract [expr int(round($q_array_b($n)*1000.00/$q_array_all($n)))]
	set b_fract [expr $b_fract/10.0]
	set c_fract [expr int(round($q_array_c($n)*1000.00/$q_array_all($n)))]
	set c_fract [expr $c_fract/10.0]
	set d_fract [expr int(round($q_array_d($n)*1000.00/$q_array_all($n)))]
	set d_fract [expr $d_fract/10.0]
	set f_fract [expr int(round($q_array_f($n)*1000.00/$q_array_all($n)))]
	set f_fract [expr $f_fract/10.0]
	set x_fract [expr int(round($q_array_x($n)*1000.00/$q_array_all($n)))]
	set x_fract [expr $x_fract/10.0]

	set sabcd_fract [expr int(round($sabcd_sum*1000.00/$q_array_all($n)))]
	set sabcd_fract [expr $sabcd_fract/10.0]
	set fx_fract [expr int(round($fx_sum*1000.00/$q_array_all($n)))]
	set fx_fract [expr $fx_fract/10.0]

	set q_sum [expr $q_array_s($n) + $q_array_a($n) + $q_array_b($n) + $q_array_c($n) + $q_array_d($n) + $q_array_f($n) + $q_array_x($n)]

	set all_fract [expr int(round($q_sum*1000.00/$q_array_all($n)))] 
	set all_fract [expr $all_fract/10.0]

	puts $f_out "$n\t$q_array_s($n)\t$q_array_a($n)\t$q_array_b($n)\t$q_array_c($n)\t$q_array_d($n)\t$q_array_f($n)\t$q_array_x($n)\t$q_array_all($n)\t\*\*\*\t$sabcd_sum\t$q_sum\t\*\*\*\t$s_fract\t$a_fract\t$b_fract\t$c_fract\t$d_fract\t$f_fract\t$x_fract\t$all_fract\t\*\*\*\t$sabcd_fract\t$fx_fract\t$all_fract"
	puts "$n\t$q_array_s($n)\t$q_array_a($n)\t$q_array_b($n)\t$q_array_c($n)\t$q_array_d($n)\t$q_array_f($n)\t$q_array_x($n)\t$q_array_all($n)\t\*\*\*\t$sabcd_sum\t$q_sum\t\*\*\*\t$s_fract\t$a_fract\t$b_fract\t$c_fract\t$d_fract\t$f_fract\t$x_fract\t$all_fract\t\*\*\*\t$sabcd_fract\t$fx_fract\t$all_fract"
	incr n
    }

	close $f_out
	puts ""
	puts "DONE"
}

if {$argc != 6} {
	puts "Program usage:"
	puts "Quality_File, Total_Length, output_file, first_Q_line, mod_line_number, line_limit"
} else {
	FastQ_Count $argv
}

