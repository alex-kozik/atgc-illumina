#!/usr/bin/tcl

proc Process_Qseq {argv} {

    ### GC CONTENT WINDOW ###
    set gc_lower 0.2
    set gc_upper 0.8

	### ADAPTOR TRIMMING STATUS ###
	set adaptor_trim "TRUE"
	# set adaptor_trim "FALSE"
	
	### LOW COMPLEXITY TRIM ###
	set lowcompl_trim "TRUE"
	# set lowcompl_trim "FALSE"
	
    ### INPUT - OUTPUT FILES ###
    set f_in1  [open [lindex $argv 0] "r"]
    set f_tab  [open [lindex $argv 1].trim.tab "w"]
    set f_fsq  [open [lindex $argv 1].trim.fastq "w"]
    set f_fst  [open [lindex $argv 1].trim.fasta "w"]

    ### PARAMETERS ###
    set tag_l  [lindex $argv 2]
    set r_len  [lindex $argv 3]
    set length_min [lindex $argv 4]

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

	### LIMIT SEQS TO SELECTED LENGTH ###
	set seq_lmt [string  range  $seq_str  0  $r_len]
	set seq_qlt [string  range  $seq_qlt  0  $r_len]

	### SELECT ILLUMINA HQ FILTERED READS ONLY ###
	if { $seq_flt == 1 } {
		### FIND B - LOW QUALITY CALLS ###
		regsub -all {B.*} $seq_qlt "" trimmed_quality
		set clean_length_q [string length $trimmed_quality]
		### FIND AND REMOVE BAD BASECALLS ###
		regsub -all {\..*} $seq_lmt "" trimmed_fasta
		
		set clean_length_t [string length $trimmed_fasta]
		
		### FIND SMALLEST ###
		if { $clean_length_q <= $clean_length_t } {
			set clean_length_a $clean_length_q
		}
		if { $clean_length_t <= $clean_length_q } {
			set clean_length_a $clean_length_t
		}
		
		### SELECT SUBSTRING ###
		# set trimmed_fasta [string range $seq_lmt  $tag_l  [expr $clean_length_a-1]]
		set trimmed_fasta [string range $seq_lmt  0  [expr $clean_length_a-1]]
		
		set adaptor_str ""
		### ADAPTOR TRIMMING ###
		if { $adaptor_trim == "TRUE" } {
			set adp_before [string length $trimmed_fasta]
			set trimmed_fasta [Adaptor_Trimming $trimmed_fasta]
			set adp_after  [string length $trimmed_fasta]
			if { $adp_after < $adp_before } {
				set adaptor_str   " _ADP_TRM_ "
			}
		}
		
		set low_compl_str ""
		### LOW COMPLEXITY TRIMMING ###
		if { $lowcompl_trim == "TRUE" } {
			set low_before [string length $trimmed_fasta]
			set trimmed_fasta [LowCompl_Trimming $trimmed_fasta]
			set low_after  [string length $trimmed_fasta]
			if { $low_after < $low_before } {
				set low_compl_str " LCMPL_TRM "
			}
		}
		
		set clean_length_s [string length $trimmed_fasta]
		
		### FIND SMALLEST ###
		if { $clean_length_q <= $clean_length_s } {
			set clean_length $clean_length_q
		}
		if { $clean_length_s <= $clean_length_q } {
			set clean_length $clean_length_s
		}
		
		### SELECT SUBSTRING ###
		set seq_trm [string range $seq_lmt  $tag_l  [expr $clean_length-1]]
		set qlt_trm [string range $seq_qlt  $tag_l  [expr $clean_length-1]]
		set trm_len [string length $seq_trm]
		### GC CONTENT ###
		set clean_upper [string toupper $seq_trm]
		regsub -all {A} $clean_upper "" clean_upper
		regsub -all {T} $clean_upper "" clean_upper
		set gc_length [string length $clean_upper]
		set gc_ratio -1
		if { $trm_len != 0 } {
			set gc_ratio [expr $gc_length*1.00/$trm_len]
		}
		set gc_status "UNKNOWN"
		if { $gc_ratio >= $gc_upper && $trm_len != 0 } {
                	set gc_status "UPPER_GC"
		}
		if { $gc_ratio <= $gc_lower && $trm_len != 0 } {
			set gc_status "LOWER_GC"
        	}
		if { $gc_ratio > $gc_lower && $gc_ratio < $gc_upper && $trm_len != 0 } {
			set gc_status "___GC___"
		}
		### TAG INFO ###
		if { $tag_l == 0 } {
			set seq_tag "_NONE_"
			set fasta_id_str "$run_id\:$machine\:$lane_n\:$tile_n\:$coordx\:$coordy\#$indx_0\/$read_n  "
		}
		if { $tag_l >= 1 } {
			set seq_tag [string range $seq_str   0 [expr $tag_l - 1]]
			set fasta_id_str "$run_id\:$machine\:$lane_n\:$tile_n\:$coordx\:$coordy\#$indx_0\/$read_n  \[ADPT\: _$seq_tag\_\]  TRM_LEN: $trm_len "
		}
		### SELECT ONLY HQ LONG READS ABOVE CUTOFF VALUE ###
		if { $clean_length >= $length_min && $gc_status == "___GC___" } {
			### TAB-DELIMITED FILE ###
			puts $f_tab "$l\t\ DATA_LINE_START @$fasta_id_str _SEQS_NEW_LINE_ $seq_trm\ +$fasta_id_str _QUAL_NEW_LINE_ $qlt_trm END_OF_DATA_LINE "
			### FAST-Q FILE ###
			puts $f_fsq "\@$fasta_id_str $l  TOT_LEN: $clean_length  \( $gc_status: $gc_length\/$trm_len \) $adaptor_str $low_compl_str"
			puts $f_fsq $seq_trm
			puts $f_fsq "\+"
			puts $f_fsq $qlt_trm
			puts $f_fsq ""
			### FAST-A FILE ###
			puts $f_fst "\>$fasta_id_str $l  TOT_LEN: $clean_length  \( $gc_status: $gc_length\/$trm_len \) $adaptor_str $low_compl_str"
			puts $f_fst $seq_trm
			puts $f_fst ""
			incr g
		}
	}

        set k_mod [expr fmod($l,1000)]
        if { $k_mod == 0 } {
                puts " - $g  OUT OF  - $l "
        }
	incr l
    }
    close $f_in1
    close $f_tab
    close $f_fsq
    close $f_fst

    puts " -+-+-+-+-+-+-+-+-+-+-+-+-+- "
    puts " - $g  OUT OF  - $l "
    puts " -+-+-+-+-+-+-+-+-+-+-+-+-+- "

    puts ""
    puts "DONE"
}

proc Adaptor_Trimming { trimmed_fasta } {
	### NGB00361.1:1-92 Illumina PCR Primer
	regsub -all -line {AGATCGGAAGAGCGTC.*} $trimmed_fasta "" trimmed_fasta
	regsub -all -line {AGATCGGAAGAGCGTC$} $trimmed_fasta "" trimmed_fasta
	regsub -all -line {AGATCGGAAGAGCGT$} $trimmed_fasta "" trimmed_fasta
	### NGB00362.1:1-61 Illumina Paired End PCR Primer 2.0
	regsub -all -line {AGATCGGAAGAGCGGT.*} $trimmed_fasta "" trimmed_fasta
	regsub -all -line {AGATCGGAAGAGCGGT$} $trimmed_fasta "" trimmed_fasta
	regsub -all -line {AGATCGGAAGAGCGG$} $trimmed_fasta "" trimmed_fasta
	### COMMON ADAPTOR ###
	regsub -all -line {AGATCGGAAGAGCG$} $trimmed_fasta "" trimmed_fasta
	regsub -all -line {AGATCGGAAGAGC$} $trimmed_fasta "" trimmed_fasta
	regsub -all -line {AGATCGGAAGAG$} $trimmed_fasta "" trimmed_fasta
	regsub -all -line {AGATCGGAAGA$} $trimmed_fasta "" trimmed_fasta
	regsub -all -line {AGATCGGAAG$} $trimmed_fasta "" trimmed_fasta
	regsub -all -line {AGATCGGAA$} $trimmed_fasta "" trimmed_fasta
	regsub -all -line {AGATCGGA$} $trimmed_fasta "" trimmed_fasta
	return $trimmed_fasta
}

proc LowCompl_Trimming { trimmed_fasta } {
	# regsub -all -line {^A{8,}} $trimmed_fasta "" trimmed_fasta
	# regsub -all -line {^T{8,}} $trimmed_fasta "" trimmed_fasta
	# regsub -all -line {^G{8,}} $trimmed_fasta "" trimmed_fasta
	# regsub -all -line {^C{8,}} $trimmed_fasta "" trimmed_fasta
	regsub -all -line {A{8,}$} $trimmed_fasta "" trimmed_fasta
	regsub -all -line {T{8,}$} $trimmed_fasta "" trimmed_fasta
	regsub -all -line {G{8,}$} $trimmed_fasta "" trimmed_fasta
	regsub -all -line {C{8,}$} $trimmed_fasta "" trimmed_fasta
	return $trimmed_fasta
}

if {$argc != 5} {
    puts "Program usage:"
    puts "qseq_file_to_process, output_file, tag_length, read_length, min_length_cutoff"
} else {
    Process_Qseq $argv
}
