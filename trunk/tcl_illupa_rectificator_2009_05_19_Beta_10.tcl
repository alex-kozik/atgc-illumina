#!/usr/bin/tcl

proc Run_Illupator { argv } {

	set dir_path [lindex $argv 0]
	set reg_expr [lindex $argv 1]
	set out_name [lindex $argv 2]
	set cycles_n [lindex $argv 3]
	set cutoff_v [lindex $argv 4]

	puts ""
	puts "    -= PARAMETERS =-    "
	puts ""
	puts "DIRECTORY:     $dir_path"
	puts "FILE_PATTERN:  $reg_expr"
	puts "OUTPUT_FILE:   $out_name"
	puts "CYCLES_N:      $cycles_n"
	puts "CUTOFF_VAL:    $cutoff_v"
	puts ""

	### FOR THIS PARTICULAR PROJECT SUFFIX MUST BE _pos.txt ###
	if { $reg_expr != "_pos.txt" } {
		puts ""
		puts "FOR THIS PARTICULAR PROJECT SUFFIX MUST BE _pos.txt"
		puts ""
		exit
	}
	### END OF DUMMY CONDITION ###

	set file_out_log [open $out_name\.Log "w"]
	set file_out_sig [open $out_name\.Sig2.Good "w"]
	set file_out_bad [open $out_name\.Sig2.Bad "w"]
	set file_out_xyz [open $out_name\.XY.Coords "w"]

	puts $file_out_log ""
	puts $file_out_log "    -= PARAMETERS =-    "
	puts $file_out_log ""
	puts $file_out_log "DIRECTORY:     $dir_path"
	puts $file_out_log "FILE_PATTERN:  $reg_expr"
	puts $file_out_log "OUTPUT_FILE:   $out_name"
	puts $file_out_log "CYCLES_N:      $cycles_n"
	puts $file_out_log "CUTOFF_VAL:    $cutoff_v"
	puts $file_out_log ""

	global pos_file_array
	global int_file_array
	global nse_file_array
	global trio_file_array
	global data_coords_array
	global data_int_array
	global data nse_array
	global cycle_array

	set pattern_len [string length $reg_expr]

	set pos_file_list [ glob -directory $dir_path -type f *$reg_expr ]
	set pos_file_list [ lsort $pos_file_list ]
	set file_count 1
	puts ""
	puts "   -= FILES TO WORK WITH =-   "
	puts ""
	puts $file_out_log ""
	puts $file_out_log "   -= FILES TO WORK WITH =-   "
	puts $file_out_log ""
	foreach pos_file_name $pos_file_list {
		set pos_file_len [ string length $pos_file_name ]
		set trio_file_array($file_count) "POSITION"
		puts "\t$file_count\t$pos_file_name\t$trio_file_array($file_count)"
		puts $file_out_log "\t$file_count\t$pos_file_name\t$trio_file_array($file_count)"
		set pos_file_array($file_count) $pos_file_name
		set prefix_file_name [ string range $pos_file_name 0 [ expr $pos_file_len - $pattern_len - 1 ] ]
		set int_file_name "$prefix_file_name\_int.txt.p"
		set nse_file_name "$prefix_file_name\_nse.txt.p"
		set int_file_array($file_count) $int_file_name
		set nse_file_array($file_count) $nse_file_name
		puts "\t\t$int_file_name"
		puts "\t\t$nse_file_name"
		puts ""
		puts $file_out_log "\t\t$int_file_name"
		puts $file_out_log "\t\t$nse_file_name"
		puts $file_out_log ""
		incr file_count
	}

	set total_count 1
	set array_count 1
	puts ""
	puts "  -= GATHERING POSITION INFO PER TILE =-  "
	puts ""
	puts $file_out_log ""
	puts $file_out_log "  -= GATHERING POSITION INFO PER TILE =-  "
	puts $file_out_log ""
	while { $array_count <= [ expr $file_count - 1 ] } {
		set temp_pos [open $pos_file_array($array_count) "r"]
		set pos_file_len [ string length $pos_file_array($array_count) ]
		### DEBUG ### puts $pos_file_len
		set tile_address [string range $pos_file_array($array_count) [expr $pos_file_len - $pattern_len-6] [expr $pos_file_len - $pattern_len-1]]
		puts ""
		puts "\t$array_count\t WORKING WITH ...  $pos_file_array($array_count)"
		puts "\t\t TILE_ADDRESS:   $tile_address"
		puts $file_out_log ""
		puts $file_out_log "\t$array_count\t WORKING WITH ...  $pos_file_array($array_count)"
		puts $file_out_log "\t\t TILE_ADDRESS:   $tile_address"
		set tile_list [ split $tile_address "_" ]
		set lane_id [ lindex $tile_list 0 ]
		set tile_id [ lindex $tile_list 1 ]
		set line_count 1
		set max_temp_line_number 1
		while { [ gets $temp_pos current_line ] >= 0 } {
			set xy_coords $current_line
			set   xy_coords  [string trim $xy_coords]
			regsub -all {  } $xy_coords " " xy_coords
			regsub -all {  } $xy_coords " " xy_coords
			regsub -all {  } $xy_coords " " xy_coords
			regsub -all {  } $xy_coords " " xy_coords
			regsub -all { } $xy_coords "\t" xy_coords
			set x_coord [ lindex $xy_coords 0 ]
			set y_coord [ lindex $xy_coords 1 ]
			set data_coords_array($line_count) "$x_coord\t$y_coord"
			puts $file_out_xyz "$total_count\t$line_count\t$lane_id\t$tile_id\t$x_coord\t$y_coord"
			incr line_count
			incr total_count
		}
		close $temp_pos
		set max_temp_line_number [ expr $line_count - 1 ]
		puts "\t\t MAX_LINE_NUMBER:  $max_temp_line_number"
		puts ""
		puts $file_out_log "\t\t MAX_LINE_NUMBER:  $max_temp_line_number"
		puts $file_out_log ""

		set temp_int [open $int_file_array($array_count) "r"]
		set temp_nse [open $nse_file_array($array_count) "r"]

		set int_line_count 1
		set int_array_count 1
		set cycle_count 1
		set cutoff_count 0
		while { [ gets $temp_int current_line ] >= 0 } {
			if { $cycle_count == 1 } {
				set cycle_array($int_array_count) 0
			}
			set comment_test [ string range $current_line 0 0 ]
			set   current_line  [string trim $current_line]
			regsub -all {  } $current_line " " current_line
			regsub -all {  } $current_line " " current_line
			regsub -all {  } $current_line " " current_line
			regsub -all {  } $current_line " " current_line
			set end_cycle_test [ string range $current_line 0 9 ]
			if { $end_cycle_test == "#END CYCLE" } {
				set int_array_count 1
				puts "INT_FILE: $current_line\tLINE: $int_line_count - END OF CYCLE:  $cycle_count"
				puts $file_out_log "INT_FILE: $current_line\tLINE: $int_line_count - END OF CYCLE:  $cycle_count"
				incr cycle_count
				### set cycle_array($int_array_count) $cutoff_count
			}
			if { $comment_test != "#" } {
				set acgt_list [split $current_line " "]
				set num_sorted [lsort -real $acgt_list]
				set max_acgt   [lindex  $num_sorted 3]
				if { $max_acgt >= $cutoff_v } {
					incr cycle_array($int_array_count)
				}
				set value_a [lindex $acgt_list 0]
				set value_c [lindex $acgt_list 1]
				set value_g [lindex $acgt_list 2]
				set value_t [lindex $acgt_list 3]
				regsub -all {\..*} $value_a "" value_a
				regsub -all {\..*} $value_c "" value_c
				regsub -all {\..*} $value_g "" value_g
				regsub -all {\..*} $value_t "" value_t
				set acgt_string "$value_a $value_c $value_g $value_t"
				set data_int_array($cycle_count,$int_array_count) $acgt_string
				incr int_line_count
				incr int_array_count
			}
		}

		puts ""
		puts $file_out_log ""

		set nse_line_count 1
		set nse_array_count 1
		set cycle_count 1
		set max_cycle 1
		set max_clust 1
		while { [ gets $temp_nse current_line ] >= 0 } {
			set comment_test [ string range $current_line 0 0 ]
			set   current_line  [string trim $current_line]
			regsub -all {  } $current_line "" current_line
			regsub -all {  } $current_line "" current_line
			regsub -all {  } $current_line "" current_line
			regsub -all {  } $current_line "" current_line
			set end_cycle_test [ string range $current_line 0 9 ]
			if { $end_cycle_test == "#END CYCLE" } {
				set nse_array_count 1
				puts "NSE_FILE: $current_line\tLINE: $nse_line_count - END OF CYCLE:  $cycle_count"
				puts $file_out_log "NSE_FILE: $current_line\tLINE: $nse_line_count - END OF CYCLE:  $cycle_count"
				set max_cycle $cycle_count
				incr cycle_count
			}
			if { $comment_test != "#" } {
				set acgt_list [split $current_line " "]
				set value_a [lindex $acgt_list 0]
				set value_c [lindex $acgt_list 1]
				set value_g [lindex $acgt_list 2]
				set value_t [lindex $acgt_list 3]
				regsub -all {\..*} $value_a " " value_a
				regsub -all {\..*} $value_c " " value_c
				regsub -all {\..*} $value_g " " value_g
				regsub -all {\..*} $value_t " " value_t
				set acgt_string "$value_a $value_c $value_g $value_t"
				set data_nse_array($cycle_count,$nse_array_count) $acgt_string
				set max_clust $nse_array_count
				incr nse_line_count
				incr nse_array_count
			}
		}

		puts ""
		puts "INTENSITY FILE LINE COUNT:  $int_line_count"
		puts "NOISE FILE LINE COUNT:      $nse_line_count"
		puts ""

		if { $nse_array_count == $int_array_count } {
			puts ""
			puts "SO FAR SO GOOD:  IDENTICAL NUMBER OF LINES IN BOTH FILES"
			puts ""
		}
		if { $nse_array_count != $int_array_count } {
			puts ""
			puts "INT NSE DATA MIS-MATCH"
			puts "      EXIT ...        "
			puts ""
			exit
		}

		set cluster_count_out 1
		while { $cluster_count_out <= $max_clust } {
			set cycle_count_out 1
			set current_coords $data_coords_array($cluster_count_out)
			set current_acgt_cycles  $cycle_array($cluster_count_out)
			if { $current_acgt_cycles >= $cycles_n } {
				puts -nonewline $file_out_sig "$lane_id\t$tile_id\t$current_coords\t"
				while { $cycle_count_out <= $max_cycle } {
					set current_int $data_int_array($cycle_count_out,$cluster_count_out)
					puts -nonewline $file_out_sig $current_int
					if { $cycle_count_out < $max_cycle } {
						puts -nonewline $file_out_sig "\t"
					}
					if { $cycle_count_out == $max_cycle } {
						puts -nonewline $file_out_sig "\n"
					}
					incr cycle_count_out
				}
			}
			if { $current_acgt_cycles < $cycles_n } {
				puts -nonewline $file_out_bad "$lane_id\t$tile_id\t$current_coords\t"
				while { $cycle_count_out <= $max_cycle } {
					set current_int $data_int_array($cycle_count_out,$cluster_count_out)
					puts -nonewline $file_out_bad $current_int
					if { $cycle_count_out < $max_cycle } {
						puts -nonewline $file_out_bad "\t"
					}
					if { $cycle_count_out == $max_cycle } {
						puts -nonewline $file_out_bad "\n"
					}
					incr cycle_count_out
				}
			}
			set l_mod [expr fmod($cluster_count_out,10000)]
			if { $l_mod == 0 } {
				puts " WORKING WITH CLUSTER:  $cluster_count_out  OUT OF:  $max_clust "
				puts " LANE: $lane_id  TILE: $tile_id "
			}
			incr cluster_count_out
		}

		puts $file_out_log ""
		puts $file_out_log "INTENSITY FILE LINE COUNT:  $int_line_count"
		puts $file_out_log "NOISE FILE LINE COUNT:      $nse_line_count"
		puts $file_out_log ""

		close $temp_int
		close $temp_nse

		incr array_count

	}

	close $file_out_log
	close $file_out_sig
	close $file_out_bad
	close $file_out_xyz

# THE END #
}

if {$argc != 5} {
	puts ""
	puts "Program usage:"
	puts "Directory_1  Pattern_2  OutputFile_3  NumberOfCycles_4  CutOffValue_5"
	puts ""
	puts "Program does analysis of Illumina  INT  NSE  POS  files"
	puts "and generates SIG2 style output for downstream basecall"
	puts ""
} else {
	Run_Illupator $argv
}

