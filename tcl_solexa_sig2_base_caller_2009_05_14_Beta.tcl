#!/usr/bin/tcl

proc Sig2_Table {argv} {

    #### PARAMETERS ####
    set upper_cut 200
    set upper_dif 60
    set good_diff 100
    set quality_upper_limit 1000
    set good_length 24
    set gc_upper 0.8
    set gc_lower 0.2
    ##  END OF PARAM  ##

    set print_bad_file "TRUE"
    ### set print_bad_file "FALSE"

    set f_in1 [open [lindex $argv 0] "r"]
    set f_out [open [lindex $argv 1] "w"]
    set f_out_qual  [open [lindex $argv 1].qual "w"]
    set f_out_clean [open [lindex $argv 1].trim.fasta "w"]
    set f_out_clean_qual  [open [lindex $argv 1].trim.qual  "w"]
    set f_out_fastq [open [lindex $argv 1].trim.fastq  "w"]
    set f_bad [open [lindex $argv 2] "w"]
    set f_bcl [open [lindex $argv 3] "w"]
    set exp_id [lindex $argv 4]
    set cycle_first [lindex $argv 5]
    set cycle_last  [lindex $argv 6]

    ### TRIMMING WITH ADAPTER OPTIONS ###
    set number_of_runs $cycle_last
    set trimming_with_adapter "FALSE"
    ### set trimming_with_adapter "TRUE"
    set adapter_seqs "TCGTAT"
    set adapter_status "_NOT_FOUND_"
    ### END OF TRIMMING WITH ADAPTER ###


    ####### READ ID LIST TO EXTRACT #######
    set k 1
    while {[gets $f_in1 current_line] >= 0} {
	set data_str [split $current_line "\t"]
	set data_len [llength $data_str]
	set line_id  [lindex $data_str 0]
	set cell_id  [lindex $data_str 1]
	set coord_x  [lindex $data_str 2]
	set coord_y  [lindex $data_str 3]

	puts "line: $k\tlength: $data_len"

	if { $k <= 10000 } {
		puts $f_bcl "\>$exp_id$line_id\_$cell_id\_$coord_x\_$coord_y  COUNT: $k"
	}

	set dna_list ""
	set qual_list ""
	set fastq_list ""

	### set n 4
	### 5-TH COLUMN CONTAIN DATA FOR FIRST CYCLE ###
	set data_list_start [expr  4 + $cycle_first - 1]
	set data_list_end   [expr  4 + $cycle_last -  1]
	set n $data_list_start

	# set first_x 0
	set first_x $number_of_runs
	set x_found "FALSE"

	while { $n <= [expr $data_list_end] } {
		set num_str  [lindex $data_str $n]
		set num_str  [string trim $num_str]
		# four times
		regsub -all {  } $num_str " " num_str
		regsub -all {  } $num_str " " num_str
		regsub -all {  } $num_str " " num_str
		regsub -all {  } $num_str " " num_str
		regsub -all { } $num_str "\t" num_str
		set val_pos  [expr $n - 3]
		# ACGT values
		set val_A [lindex $num_str 0]
		set val_C [lindex $num_str 1]
		set val_G [lindex $num_str 2]
		set val_T [lindex $num_str 3]
		###
		# puts $num_str
		set num_sorted [lsort -real $num_str]
		# puts $num_sorted
		set value_max1 [lindex $num_sorted 3]
		set value_max2 [lindex $num_sorted 2]
		set value_max_diff [expr $value_max1 - $value_max2]
		### QUALITY SCORES ### 2008 March 25
		set current_quality "X"
		set fastq_quality "!"
		set value_ratio     "X.XXXXXXXX"
		if { $value_max1 != 0 } {
			set value_ratio [expr $value_max2/$value_max1]
			if { $value_ratio <= 0.2 } {
				set current_quality "A"
				set fastq_quality "I"
			}
			if { $value_ratio <= 0.4 && $value_ratio > 0.2 } {
				set current_quality "B"
				set fastq_quality "9"
			}
			if { $value_ratio <= 0.6 && $value_ratio > 0.4 } {
				set current_quality "C"
				set fastq_quality "5"
			}
			if { $value_ratio <= 0.8 && $value_ratio > 0.6 } {
				set current_quality "D"
				set fastq_quality "0"
			}
			if { $value_ratio > 0.8 } {
				set current_quality "F"
				set fastq_quality "+"
			}
			if { $value_ratio > 0.9 } {
				set fastq_quality "#"
			}
		}
		if { $value_max1 >= $quality_upper_limit } {
			set current_quality [string toupper $current_quality]
		}
		if { $value_max1 <  $quality_upper_limit } {
			set current_quality [string tolower $current_quality]
		}
		if { $value_max1 <  $upper_cut } {
			set current_quality "x"
		}
		###
		set value_status "NONE"
		set base_call "N"
		# set first_x 0
		# set x_found "FALSE"
		if { $val_A >= $upper_cut || $val_C >= $upper_cut || $val_G >= $upper_cut || $val_T >= $upper_cut } {
			set value_status "HOPE"
			set base_call "X"
			if { $val_A >= [expr $val_C + $upper_dif] && $val_A >= [expr $val_G + $upper_dif] && $val_A >= [expr $val_T + $upper_dif] } {
				set base_call "A"
				set value_status "GOOD"
			}
			if { $val_C >= [expr $val_A + $upper_dif] && $val_C >= [expr $val_G + $upper_dif] && $val_C >= [expr $val_T + $upper_dif] } {
                                set base_call "C"
				set value_status "GOOD"
                        }
			if { $val_G >= [expr $val_A + $upper_dif] && $val_G >= [expr $val_C + $upper_dif] && $val_G >= [expr $val_T + $upper_dif] } {
                                set base_call "G"
				set value_status "GOOD"
                        }
			if { $val_T >= [expr $val_A + $upper_dif] && $val_T >= [expr $val_C + $upper_dif] && $val_T >= [expr $val_G + $upper_dif] } {
                                set base_call "T"
				set value_status "GOOD"
                        }
		}
		if { $value_max_diff < $good_diff } {
			set base_call [string tolower $base_call]
		}
		if { $val_A == 0.0 && $val_C == 0.0 && $val_G == 0.0 && $val_T == 0.0 } {
			set value_status "CRAP"
			set base_call "*"
		}
		if { $base_call == "*" && $x_found == "FALSE" } {
			set first_x [expr $n - 3]
			set x_found "TRUE"
		}

		### PUT ALL TOGETHER ###
		set dna_list [ lappend dna_list $base_call ]
		set qual_list [ lappend qual_list $current_quality ]
		set fastq_list [ lappend fastq_list $fastq_quality ]

		if { $first_x == $number_of_runs } {
			set cut_q "--"
		}
		if { $first_x != $number_of_runs } {
			set cut_q $first_x
		}

		if { $k <= 10000 } {
			puts $f_bcl "$val_pos\:\tA:$val_A\tC:$val_C\tG:$val_G\tT:$val_T\t$value_status\tBASE_CALL:-> $base_call\tFIRST_FAIL: $cut_q  MAX_DIFF: $value_max_diff  RATIO: $value_ratio  QUALITY: $current_quality"
		}
		incr n
	}

	# puts $f_out "$num_str"
	set dna_string [join $dna_list ""]
	set qual_string [join $qual_list ""]
	set fastq_string [join $fastq_list ""]
	set dna_clean $dna_string

	set string_quality [string range $dna_string 0 [expr $good_length-1]]

	regsub -all {n} $string_quality "" string_quality
	regsub -all {N} $string_quality "" string_quality
	regsub -all {x} $string_quality "" string_quality
	regsub -all {X} $string_quality "" string_quality
	regsub -all {\*} $string_quality "" string_quality

	regsub -all {n.*} $dna_clean "" dna_clean
	regsub -all {N.*} $dna_clean "" dna_clean
	regsub -all {x.*} $dna_clean "" dna_clean
	regsub -all {X.*} $dna_clean "" dna_clean
	regsub -all {\*.*} $dna_clean "" dna_clean 

	####### SPECIAL CASE WITH ADAPTER TRIMMING #######
	if { $trimming_with_adapter == "TRUE" } {
		set dna_clean [string toupper $dna_clean]
		set adapter_status "_NOT_FOUND_"
		set adapter_found [string first $adapter_seqs $dna_clean]
		if { $adapter_found >= 0 } {
			set adapter_status "_$adapter_seqs\_"
		}
		regsub -all {$adapter_seqs.*} $dna_clean "" dna_clean
	}
	##################################################

	set clean_length [string length $dna_clean]

	set clean_upper [string toupper $dna_clean]

	regsub -all {A} $clean_upper "" clean_upper
	regsub -all {T} $clean_upper "" clean_upper

	set gc_length [string length $clean_upper]

	set gc_ratio -1
	set gc_status "UNDEFINED"

	if { $clean_length != 0 } {
		set gc_ratio [expr $gc_length*1.00/$clean_length]
	}

	set gc_status "UNKNOWN"

	if { $gc_ratio >= $gc_upper && $clean_length != 0 } {
		set gc_status "UPPER_GC"
	}

	if { $gc_ratio <= $gc_lower && $clean_length != 0 } {
		set gc_status "LOWER_GC"
	}

	if { $gc_ratio > $gc_lower && $gc_ratio < $gc_upper && $clean_length != 0 } {
		set gc_status "___GC___"
	}

	# puts $string_quality

	set string_quality_len [string length $string_quality]

	set string_qual "UNKNOWN"
	if { $string_quality_len >= $good_length } {
		set string_qual "LONG__Q: $string_quality_len"
	}
	if { $string_quality_len <  $good_length } {
		set string_qual "SHORT_Q  $string_quality_len"
	}

	set fail_str $first_x

	if { $fail_str == $number_of_runs  } {
		set fail_str "--"
	}

	if { $clean_length >= $good_length } {

		set clean_quality_string [string range $qual_string 0 [expr $clean_length - 1]]
		set fastq_string_clean [string range $fastq_string 0 [expr $clean_length - 1]]
		set fffff_quality_string $clean_quality_string
		regsub -all {f} $fffff_quality_string "" fffff_quality_string
		regsub -all {F} $fffff_quality_string "" fffff_quality_string
		set abcd_qual_length [string length $fffff_quality_string]
		set ffff_qual_length [expr $clean_length - $abcd_qual_length]
		set abcd_qual_fract  [expr ($abcd_qual_length*1.00/$clean_length)*100]
		set abcd_qual_fract  [expr ceil($abcd_qual_fract)]
		set abcd_qual_fract  [expr  int($abcd_qual_fract)]

		puts $f_out_clean "\>$exp_id$line_id\_$cell_id\_$coord_x\_$coord_y  |  LINE: $k  \[ LENGTH: $clean_length \]  \( $gc_status: $gc_length\/$clean_length \)  \[ ABCD_Q:$abcd_qual_length  F_QUAL:$ffff_qual_length \]  \[ Q_FRACT:$abcd_qual_fract \]  "
		puts $f_out_clean "$dna_clean"
		puts $f_out_clean ""

		puts $f_out_clean_qual "\>$exp_id$line_id\_$cell_id\_$coord_x\_$coord_y  |  \( $gc_status: $gc_length\/$clean_length \)  \[ ABCD_Q:$abcd_qual_length  F_QUAL:$ffff_qual_length \]  \[ Q_FRACT:$abcd_qual_fract \]  "
		puts $f_out_clean_qual $clean_quality_string
		puts $f_out_clean_qual ""

		if { $gc_status == "___GC___" } {
			if { $trimming_with_adapter == "TRUE" } {
				puts $f_out_fastq "\@$exp_id$line_id\_$cell_id\_$coord_x\_$coord_y  $adapter_status \: $adapter_found"
			}
			if { $trimming_with_adapter == "FALSE" } {
				puts $f_out_fastq "\@$exp_id$line_id\_$cell_id\_$coord_x\_$coord_y"
				set dna_clean [string toupper $dna_clean]
			}
			puts $f_out_fastq $dna_clean
			puts $f_out_fastq "\+"
			puts $f_out_fastq $fastq_string_clean
			puts $f_out_fastq ""
		}

	}

	if { $first_x >= [expr $good_length+1] } {

		# COUNT ABCD and F QUALITY #
		set qall_length [string length $qual_string]
		set qfff_string $qual_string
		regsub -all {x} $qfff_string "" qfff_string
		regsub -all {X} $qfff_string "" qfff_string
		regsub -all {f} $qfff_string "" qfff_string
		regsub -all {F} $qfff_string "" qfff_string

		set abcd_length [string length $qfff_string]
		set qfx_length  [expr $qall_length - $abcd_length]
		set abcd_fract  [expr ($abcd_length*1.00/$qall_length)*100]
		set abcd_fract  [expr ceil($abcd_fract)]
		set abcd_fract  [expr  int($abcd_fract)]

		if { $abcd_length >= 20 } {

			puts $f_out "\>$exp_id$line_id\_$cell_id\_$coord_x\_$coord_y  |  LINE: $k  \[ FAIL: $fail_str \]  ::  $string_qual :: CLEAN_LENGTH: $clean_length  \[ ABCD_Q:$abcd_length  F_QUAL:$qfx_length \]  \[ Q_FRACT:$abcd_fract \]  "
			puts $f_out "$dna_string"
			puts $f_out ""

			puts $f_out_qual "\>$exp_id$line_id\_$cell_id\_$coord_x\_$coord_y  \[ ABCD_Q:$abcd_length  F_QUAL:$qfx_length \]  \[ Q_FRACT:$abcd_fract \]  "
			puts $f_out_qual "$qual_string"
			puts $f_out_qual ""

		}

		if { $abcd_length <  20 && $print_bad_file == "TRUE" } {
			puts $f_bad "\>$exp_id$line_id\_$cell_id\_$coord_x\_$coord_y  |  LINE: $k  \[ FAIL: $fail_str \]  ::  $string_qual :: CLEAN_LENGTH: $clean_length "
			puts $f_bad "$dna_string"
			puts $f_bad ""
		}

	}

    if { $first_x  < [expr $good_length+1] && $print_bad_file == "TRUE" } {
        puts $f_bad "\>$exp_id$line_id\_$cell_id\_$coord_x\_$coord_y  |  LINE: $k  \[ FAIL: $fail_str \]  ::  $string_qual :: CLEAN_LENGTH: $clean_length "
		puts $f_bad "$dna_string"
		puts $f_bad ""
    }

	if { $k <= 10000 } {
		puts $f_bcl ""
	}

	incr k
    }

    close $f_in1
    close $f_out
    close $f_out_qual
    close $f_out_clean
    close $f_out_clean_qual
    close $f_out_fastq
    close $f_bad
    close $f_bcl
    puts ""
    puts "DONE"
}

if {$argc != 7} {
    puts "Program usage:"
    puts "Table_to_Process, output_file_good, output_file_bad, output_file_debug, id_prefix, cycle_first, cycle_last"
} else {
    Sig2_Table $argv
}

