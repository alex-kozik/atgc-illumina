#!/usr/bin/tcl

proc Process_Qseq {argv} {

    set length_cut 60

    set f_in1  [open [lindex $argv 0] "r"]
    set f_out  [open [lindex $argv 1] "w"]
    set f_fas  [open [lindex $argv 1].fa "w"]
    set f_trm  [open [lindex $argv 1].trim.fa "w"]
    set tag_l  [lindex $argv 2]
    set r_len  [lindex $argv 3]

    ####### READ AND PROCESS QSEQ DATA #######
    set l 0
    set g 0
    while { [gets $f_in1 current_line] >= 0 } {
	set current_data [split $current_line "\t"]
	set machine [lindex $current_data  0]
        set run_id  [lindex $current_data  1]
	set lane_n  [lindex $current_data  2]
	set tile_n  [lindex $current_data  3]
	set coordx  [lindex $current_data  4]
	set coordy  [lindex $current_data  5]
	set indx_0  [lindex $current_data  6]
	set read_n  [lindex $current_data  7]
	set seq_str [lindex $current_data  8]
	set seq_qlt [lindex $current_data  9]
	set seq_flt [lindex $current_data 10]
	if { $tag_l >= 1 } {
		set seq_tag [string range $seq_str   0 [expr $tag_l - 1]]
		set seq_trm [string range $seq_str   $tag_l       $r_len]
		set qlt_trm [string range $seq_qlt   $tag_l       $r_len]
		set fasta_id_str "$run_id\:$machine\:$lane_n\:$tile_n\:$coordx\:$coordy\#$indx_0\/$read_n  \[ADPT\: _$seq_tag\_\]  "
	}
	if { $tag_l == 0 } {
		set seq_tag "_NONE_"
		set seq_trm $seq_str
		set qlt_trm $seq_qlt
		set fasta_id_str "$run_id\:$machine\:$lane_n\:$tile_n\:$coordx\:$coordy\#$indx_0\/$read_n  "
	}
	set crap_found [string first "\." $seq_str]
	if { $crap_found < 0 && $seq_flt == 1 } {
		### set fasta_id_str "$run_id\:$machine\:$lane_n\:$tile_n\:$coordx\:$coordy\#$indx_0\/$read_n  \[ADPT\: _$seq_tag\_\]  "
		puts $f_out "$l\t\ DATA_LINE_START @$fasta_id_str _SEQS_NEW_LINE_ $seq_trm\ +$fasta_id_str _QUAL_NEW_LINE_ $qlt_trm END_OF_DATA_LINE "
		puts $f_fas "\>$fasta_id_str $l  "
		puts $f_fas $seq_trm
		puts $f_fas ""
		incr g
	}

	regsub -all {\..*} $seq_trm "" trimmed_fasta
	set clean_length [string length $trimmed_fasta]
	if { $clean_length >= $length_cut && $seq_flt == 1 } {
		puts $f_trm "\>$fasta_id_str $l  LEN: $clean_length "
		puts $f_trm $trimmed_fasta
		puts $f_trm ""
	}

        set k_mod [expr fmod($l,1000)]
        if { $k_mod == 0 } {
                puts " - $g  OUT OF  - $l "
        }
	incr l
    }
    close $f_in1
    close $f_out
    close $f_fas
    close $f_trm

    puts ""
    puts "DONE"
}

if {$argc != 4} {
    puts "Program usage:"
    puts "qseq_file_to_process, output_file, tag_length, read_length"
} else {
    Process_Qseq $argv
}

