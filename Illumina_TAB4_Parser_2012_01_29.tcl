#!/usr/bin/tcl

proc Process_Qseq {argv} {
	
	set mod_val 100000
	
	### GC CONTENT WINDOW ###
	set gc_lower 0.2
	set gc_upper 0.8
	
	### ADAPTOR TRIMMING STATUS ###
	set adaptor_trim "TRUE"
	# set adaptor_trim "FALSE"
	
	### LOW COMPLEXITY TRIM ###
	set lowcompl_trim "TRUE"
	# set lowcompl_trim "FALSE"
	
	### DNA TAG TRIM ###
	# set dna_tag_trim "TRUE"
	set dna_tag_trim "FALSE"
	
	### TAB FILE WITH TRIMMED QSEQ DATA ###
	set tab_data "FALSE"
	# set tab_data "TRUE"
	
	### USE BAR-CODE LIBRARY ###
	# set bar_code_lib_status "TRUE"
	set bar_code_lib_status "FALSE"
	# set add_dummy_T "TRUE"
	set add_dummy_T "FALSE"
	global bar_code_array
	# set rev_comp_bar "TRUE"
	set rev_comp_bar "FALSE"
	
	### INPUT  FILES ###
	### QSEQ FILE ###
	set f_in1  [open [lindex $argv 0] "r"]
	### BAR CODE DICTIONARY ###
	if { $bar_code_lib_status == "TRUE" } {
		puts "                           "
		puts " -= BAR CODE DICTIONARY =- "
		set d 0
		set f_in2  [open [lindex $argv 0]._bar_codes_lib.tab "r"]
		while { [gets $f_in2 current_line] >= 0 } {
			set current_data [split $current_line "\t"]
			set bar_code_seq [lindex $current_data  0]
			if { $rev_comp_bar == "TRUE" } {
				set bar_code_seq [Reverse_Complement_String $bar_code_seq]
			}
			### ADD T TO BAR CODE SEQUENCE ###
			if { $add_dummy_T == "TRUE" } {
				set bar_code_seq "$bar_code_seq\T"
			}
			set bar_code_id  [lindex $current_data  1]
			incr d
			puts "$bar_code_seq\t$bar_code_id\t -=- \t$d"
			set bar_code_array($bar_code_seq) $bar_code_id
		}
		close $f_in2
		puts " -= +++++++++++++++++++ =- "
		puts "                           "
	}
	
	### OUTPUT FILES ###
	set f_seq_1  [open [lindex $argv 1]._all_1.fastq "w"]		; # ALL HIGH QUALITY 1 QSEQ DATA CONVERTED TO FASTQ
	if { $tab_data == "TRUE" } {
		set f_tab  [open [lindex $argv 1].trim.tab "w"]
	}
	set f_fsq  [open [lindex $argv 1].trim.fastq "w"]
	set f_fst  [open [lindex $argv 1].trim.fasta "w"]
	set f_log  [open [lindex $argv 1].trim.xlog "w"]
	set f_gch  [open [lindex $argv 1].trim.xgch.fa "w"]		; # HIGH GC CONTENT
	set f_gcl  [open [lindex $argv 1].trim.xgcl.fa "w"]		; # LOW GC CONTENT
	
	### PARAMETERS ###
	set tag_l  [lindex $argv 2]
	set r_len  [lindex $argv 3]
	set length_min [lindex $argv 4]
	set input_format [lindex $argv 5]
	set bad_q_char [lindex $argv 6]
	
	puts $f_log  "TAG__L: $tag_l"
	puts $f_log  "READ_L: $r_len"
	puts $f_log  "MIN__L: $length_min"
	
	####### READ AND PROCESS QSEQ DATA #######
	set l 0
	set g 0
	set adp_count 0		; # ADAPTOR COUNT
	set hpl_count 0		; # HOMOPOLYMER COUNT
	set dna_count 0		; # DNA TAG COUNT
	set gc_low_count 0
	set gc_high_count 0
	set gc_ok_count 0
	
	while { [gets $f_in1 current_line] >= 0 } {
	set current_data [split $current_line "\t"]
	if { $input_format == "tab-4" } {
		set spot_id  [lindex $current_data  0]
		set seq_str  [lindex $current_data  1]
		set qual_id  [lindex $current_data  2]
		set seq_qlt  [lindex $current_data  3]
	}
	if { $input_format != "tab-4" } {
		puts "Wrong File Format"
		exit
	}
	
	### LIMIT SEQS TO SELECTED LENGTH ###
	set seq_lmt [string  range  $seq_str  0  $r_len]
	set seq_qlt [string  range  $seq_qlt  0  $r_len]
	
	
	set fasta_id_head [string  range  $spot_id  1  end]
	### puts $f_seq_1 "\@$fasta_id_head  $l "
	### regsub -all {\.} $seq_lmt "N" seq_lmt_N
	### puts $f_seq_1 $seq_lmt_N
	### puts $f_seq_1 "\+ "
	### puts $f_seq_1 $seq_qlt
	
	### FIND B - LOW QUALITY CALLS ###
	if { $bad_q_char == "B" } {
		regsub -all {B.*} $seq_qlt "" trimmed_quality
	}
	if { $bad_q_char == "!" } {
		regsub -all {\!.*} $seq_qlt "" trimmed_quality
		regsub -all {\".*} $seq_qlt "" trimmed_quality
		regsub -all {\#.*} $seq_qlt "" trimmed_quality
	}
	
	set clean_length_q [string length $trimmed_quality]
	### FIND AND REMOVE BAD BASECALLS ###
	regsub -all {\..*} $seq_lmt "" trimmed_fasta
	regsub -all {N.*} $seq_lmt "" trimmed_fasta
	
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
			incr adp_count
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
			incr hpl_count
		}
	}
	
	set dna_tag_str ""
	### DNA TAG TRIM ###
	if { $dna_tag_trim == "TRUE" } {
		set dna_before [string length $trimmed_fasta]
		set trimmed_fasta [DNA_Tag_Trimming $trimmed_fasta]
		set dna_after  [string length $trimmed_fasta]
		if { $dna_after < $dna_before } {
			set dna_tag_str "_TCGTATGCC_"
			incr dna_count
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
	if { $trm_len >= $length_min } {
		set gc_ratio [expr $gc_length*1.00/$trm_len]
		set gc_intgr [expr int(round($gc_ratio*100))]
	}
	set gc_status "UNKNOWN"
	if { $gc_ratio > $gc_upper && $trm_len >= $length_min } {
		set gc_status "HIGH_GC"
		incr gc_high_count
	}
	if { $gc_ratio < $gc_lower && $trm_len >= $length_min } {
		set gc_status "LOW__GC"
		incr gc_low_count
		}
	if { $gc_ratio >= $gc_lower && $gc_ratio <= $gc_upper && $trm_len >= $length_min } {
		set gc_status "___GC___"
		incr gc_ok_count
	}
	### TAG INFO ###
	if { $tag_l == 0 } {
		set seq_tag "_NONE_"
		set fasta_id_str $fasta_id_head
	}
	if { $tag_l >= 1 } {
		set seq_tag [string range $seq_str   0 [expr $tag_l - 1]]
		if { $bar_code_lib_status == "FALSE" } {
			set fasta_id_str "$fasta_id_head  \[ADPT\: _$seq_tag\_\]  TRM_LEN: $trm_len "
		}
		if { $bar_code_lib_status == "TRUE" } {
			set bar_code_query [info exists bar_code_array($seq_tag)]
			if { $bar_code_query == 0 } {
				set bar_code_id "XXXX"
			}
			if { $bar_code_query == 1 } {
				set bar_code_id $bar_code_array($seq_tag)
			}
			set fasta_id_str "_$bar_code_id\_\:$fasta_id_head  \[ADPT\: _$seq_tag\_\]  TRM_LEN: $trm_len "
		}
	}
	
	### SELECT ONLY HQ LONG READS ABOVE CUTOFF VALUE ###
	if { $clean_length >= $length_min && $gc_status == "___GC___" } {
		### TAB-DELIMITED FILE ###
		if { $tab_data == "TRUE" } {
			puts $f_tab "$l\t\ DATA_LINE_START @$fasta_id_str _SEQS_NEW_LINE_ $seq_trm\ +$fasta_id_str _QUAL_NEW_LINE_ $qlt_trm END_OF_DATA_LINE "
		}
		### FAST-Q FILE ###
		# puts $f_fsq "\@$fasta_id_str $l  TOT_LEN: $clean_length  \( $gc_status: $gc_length\/$trm_len \) $adaptor_str $low_compl_str"
		puts $f_fsq "\@$fasta_id_str $l  TOT_LEN: $clean_length  \( $gc_status: $gc_intgr \) $adaptor_str $low_compl_str $dna_tag_str"
		puts $f_fsq $seq_trm
		puts $f_fsq "\+ "
		puts $f_fsq $qlt_trm
		puts $f_fsq ""
		### FAST-A FILE ###
		# puts $f_fst "\>$fasta_id_str $l  TOT_LEN: $clean_length  \( $gc_status: $gc_length\/$trm_len \) $adaptor_str $low_compl_str"
		puts $f_fst "\>$fasta_id_str $l  TOT_LEN: $clean_length  \( $gc_status: $gc_intgr \) $adaptor_str $low_compl_str $dna_tag_str"
		puts $f_fst $seq_trm
		puts $f_fst ""
		incr g
	}
	
	if { $clean_length >= $length_min && $gc_status == "HIGH_GC" } {
			puts $f_gch "\>$fasta_id_str $l  TOT_LEN: $clean_length  \( $gc_status: $gc_intgr \) $adaptor_str $low_compl_str $dna_tag_str"
			puts $f_gch $seq_trm
			puts $f_gch ""
		}
	
	if { $clean_length >= $length_min && $gc_status == "LOW__GC" } {
		puts $f_gcl "\>$fasta_id_str $l  TOT_LEN: $clean_length  \( $gc_status: $gc_intgr \) $adaptor_str $low_compl_str $dna_tag_str"
		puts $f_gcl $seq_trm
		puts $f_gcl ""
	}
		set k_mod [expr fmod($l,$mod_val)]
		if { $k_mod == 0 } {
			puts " - $g  OUT OF  - $l    - ADP_TRM: $adp_count    - LCMPL_TRM: $hpl_count    - BAR_COUNT: $dna_count"
			puts "LOW_GC: $gc_low_count      HIGH_GC: $gc_high_count      OK_GC: $gc_ok_count"
			puts " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - "
		}
	incr l
	}
	
	close $f_in1
	if { $tab_data == "TRUE" } {
		close $f_tab
	}
	close $f_fsq
	close $f_fst
	close $f_gch
	close $f_gcl
	close $f_seq_1
	
	puts " -+-+-+-+-+-+-+-+-+-+-+-+-+- "
	puts " - $g  OUT OF  - $l    - ADP_TRM: $adp_count    - LCMPL_TRM: $hpl_count    - BAR_COUNT: $dna_count"
	puts "LOW_GC: $gc_low_count      HIGH_GC: $gc_high_count      OK_GC: $gc_ok_count"
	puts $f_log " - $g  OUT OF  - $l    - ADP_TRM: $adp_count    - LCMPL_TRM: $hpl_count    - BAR_COUNT: $dna_count"
	puts $f_log "LOW_GC: $gc_low_count      HIGH_GC: $gc_high_count      OK_GC: $gc_ok_count"
	puts " -+-+-+-+-+-+-+-+-+-+-+-+-+- "
	
	close $f_log
	
	puts ""
	puts "DONE"
}

proc Adaptor_Trimming { trimmed_fasta } {
	### NGB00361.1:1-92 Illumina PCR Primer
	# regsub -all -line {AGATCGGAAGAGCGTC.*} $trimmed_fasta "" trimmed_fasta
	# regsub -all -line {AGATCGGAAGAGCGTC$} $trimmed_fasta "" trimmed_fasta
	# regsub -all -line {AGATCGGAAGAGCGT$} $trimmed_fasta "" trimmed_fasta
	### NGB00362.1:1-61 Illumina Paired End PCR Primer 2.0
	# regsub -all -line {AGATCGGAAGAGCGGT.*} $trimmed_fasta "" trimmed_fasta
	# regsub -all -line {AGATCGGAAGAGCGGT$} $trimmed_fasta "" trimmed_fasta
	# regsub -all -line {AGATCGGAAGAGCGG$} $trimmed_fasta "" trimmed_fasta
	### COMMON ADAPTOR ###
	regsub -all -line {AGATCGGAAGAGCG.*} $trimmed_fasta "" trimmed_fasta
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

proc DNA_Tag_Trimming { trimmed_fasta } {
	regsub -all -line {TCGTATGCC.*} $trimmed_fasta "" trimmed_fasta
	regsub -all -line {TCGTATGCC$} $trimmed_fasta "" trimmed_fasta
	return $trimmed_fasta
}

proc Reverse_Complement_String { string_frw } {
	
	global sleep_time
	
	### set string_rev [string reverse $string_frw]		; # Tcl version 8.5
	
	set string_rev {}
	set i [string length $string_frw]
	while { $i > 0 } {
		append string_rev [string index $string_frw [incr i -1]]
	}
	
	set string_rev_compl [string map { A T G C C G T A } $string_rev]
	return $string_rev_compl
	
}

if {$argc != 7} {
	puts "Program usage:"
	puts "qseq_file_to_process, output_file, tag_length, read_length, min_length_cutoff, tab-4, bad-B/bad-!"
} else {
	Process_Qseq $argv
}
