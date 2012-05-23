#!/usr/bin/tcl

proc Process_Tables {argv} {
	
	set mod_value 100000
	
	set f_in1 [open [lindex $argv 0] "r"]
	set f_in2 [open [lindex $argv 1] "r"]
	set f_out [open [lindex $argv 2] "w"]
	set condt [lindex $argv 3]
	
	global id_array
	
	####### READ ID TABLE INTO MEMORY #######
	set l 0
	while { [gets $f_in2 current_line] >= 0 } {
		set current_data [split   $current_line "\t"]
		set current_id   [lindex  $current_data 0]
		set id_array($current_id) $current_line
		incr l
		set n_mod [expr fmod($l,$mod_value)]
		if { $n_mod == 0 } {
			puts " ... $l ... "
		}
	}
	close $f_in2
	puts " $l IDs found "
	puts "--------------------"
	
	####### READ DATA TABLE TO EXTRACT #######
	set k 0
	set m 0
	set e 0
	while {[gets $f_in1 current_line] >= 0} {
		set current_data   [split   $current_line "\t"]
		set current_id     [lindex  $current_data 0]
		## set current_id     [lindex  $current_data 1]
	
		set query_id [info exists id_array($current_id)]
		if {$query_id == 0} {
			## puts "  $current_id  WAS NOT EXTRACTED ???  *** $k LINES  "
			if { $condt == "NOT" } {
				puts $f_out $current_line
			}
			incr m
		}
		if {$query_id == 1} {
			## puts "  $current_id  WAS     EXTRACTED !!!  *** $k LINES  "
			if { $condt == "YES" } {
				puts $f_out $current_line
			}
			incr e
		}
		incr k
		set n_mod [expr fmod($k,$mod_value)]
		if { $n_mod == 0 } {
			puts " $k lines processed - $e lines extracted "
		}
	}
	
	puts "$m DATAPOINT OUT OF $k WERE NOT EXTRACTED"
	puts "$e DATAPOINT OUT OF $k WERE     EXTRACTED"
	
	close $f_in1
	close $f_out
	puts ""
	puts "DONE"
}

if {$argc != 4} {
	puts "Program usage:"
	puts "Data_File_Table, List_of_IDs, output_file, condition_YES_NOT"
} else {
    Process_Tables $argv
}
