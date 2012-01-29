#!/usr/bin/tcl

proc FastQ2Tab {argv} {
	
	set f_in1 [open [lindex $argv 0] "r"]
	set f_out [open [lindex $argv 1] "w"]
	set l_mod [lindex $argv 2]
	
	set l 1
	
	while { [gets $f_in1 current_line] >= 0 } {
		set line_mod [expr fmod($l,$l_mod)]
		puts -nonewline $f_out $current_line
		set line_delimeter \t
		if { $line_mod == 0 && $l != 0 } {
			set line_delimeter \n
		}
		puts -nonewline $f_out $line_delimeter
		set progress_mod [expr fmod($l,100000)]
		if { $progress_mod == 0 } {
			puts " ... $l lines processed ... "
		}
		incr l
	}
	close $f_out
	close $f_in1
	puts ""
	puts " Well Done - $l Lines Total "
	puts ""
}

if {$argc != 3} {
	puts "Program usage:"
	puts "input_file, output_file, mod_line_number"
} else {
	FastQ2Tab $argv
}
