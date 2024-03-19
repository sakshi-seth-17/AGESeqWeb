#!/usr/bin/perl

use strict;
use warnings;
use Storable;
use HTTP::Server::Simple::CGI;
use JSON::XS;# qw(decode_json);
use Time::Piece;
use Cwd;

my $port = 8087;
my $hostname = "172.30.18.104";




{
    package MyHTTPServer;
    use base qw(HTTP::Server::Simple::CGI);
    use CGI;
	
	sub ageSeq {
	
	#my ($dat_dir, $file_design, $final_out) = @_;  # Retrieve the arguments
	my ($dat_dir, $file_design, $final_out, $mismatch_cutoff, $min_cutoff, $wt_like_report, $indel_report ) = @_;  # Retrieve the arguments
	
	print "\ninside fun ".$file_design;
	
	my $blat = '';   # working directory or PATH
	#my $blat = '/usr/local/blat/latest/bin/blat';  # Your prefered full location
	
	##########################
	## setting can be changed here

	## setting for reports
	#my $mismatch_cutoff   = 0.1 ;    # mismatch rate to filter low quality alignment, default 0.1 (10 %)
	#my $min_cutoff        = 0   ;    # cutoff to filter reads with low abundance, default  0
	#my $wt_like_report    = 20   ;    # report top xx WT like records, default 20 
	#my $indel_report      = 50   ;    # report top xx records with indel, default  50 
	my $remove_files      = 1  ;  # keep (0) or delete (1) intermediate files, default = 1
	
	###############################################################
	# default setting
	
	my $remove = 'rm';
	my $osname = $^O;
	
	print "You are running ",$osname," system\n";
	if(uc($osname )=~ m/MSWIN/){
		$blat  = 'blat.exe';
		$remove = 'del';
	}
	else{
		if($blat eq ''){
			# use default
			if(-e 'blat'){
				$blat = './blat';
			}
			else{
				$blat = 'blat';
			}
		}
	}
	
	
	my  $design_fas = $file_design.'_DESIGN.fa';
	my @psl_file_data = ();
	
	if(not defined ($file_design )){
		die "Design file is needed\n";
	}
	
	
	if(not defined ($final_out )){
		$final_out = 'AGE_output.txt';
	}
	##########################################################
	###  step 1 load design file
	
	my $link_br = &check_line_break($file_design);
	
	open DESIGN,"<$file_design" or die "File  $file_design  not found error here\n";
	open DESIGNFAS,">$design_fas" or die "File  $file_design  not found error here\n";
	
	my $num = 0;
	my %design_hash = ();
	local $/ = $link_br;
	
	while (my $line = <DESIGN>) {
		$num ++;
		if($num == 1){
			next;
		}
		if($line !~ m/\w/){
			next;
		}
		
		chomp  $line;
		my @temp = split(/\t/,$line);
		my $id  = $temp[0];
		my $seq = $temp[1];
		$seq =~ s/\s+//g;
		$seq =~ m/([^\r\n]+)/;
		$seq = $1;
		$id  =~ s/\s+/\_/g;
		
		$design_hash{$id} = $seq;
		print DESIGNFAS ">$id\n$seq\n";
	}
	close(DESIGNFAS);
	
		
	#######################################################
	##############  step 2 - load read files  #############
	
	#### fastq files
	opendir(ARGV,"$dat_dir");
	my @files=grep(/\.fq$|\.fastq$/,readdir ARGV);
	my %total_read_num = ();
	my $num_c          = 0; 
	my @fasta_files    = ();
	
	foreach my $file(@files){	
			my $file_in   = $dat_dir.'/'.$file;
			if(uc($osname)=~ m/MSWIN/){
				$file_in   = $dat_dir."\\".$file;
			}
			else{
				if((substr $file_in,0,1) ne '/'){
					$file_in = './'.$file_in;
				}
			}
					
			my $fasta_out = $file.'_INPUT.fa';         
			
			push  @fasta_files   , $fasta_out; ## from fastq
	
			my $check_fastq2fasta = &fastq2fasta($file_in,$fasta_out);
			my $total_num  = $check_fastq2fasta;
			$total_read_num{$fasta_out} = $total_num;
			print "$total_num read in this $file\n";
	}
	
	###  fasta files
	opendir(ARGV,"$dat_dir");
	@files=grep(/\.txt$|\.fas$|\.fa$|\.fasta$|\.seq$/,readdir ARGV);
	
	my @single_fas = ();
	foreach my $file(@files){	
			my $file_in   = $dat_dir.'/'.$file;
			
			if(uc($osname )=~ m/MSWIN/){
				$file_in   = $dat_dir."\\".$file;
			}	 
			
			### check read num
			my $seq_num = check_read_num ($file_in);
			
			if($seq_num>1){
				my $fasta_out = $file.'_INPUT.fa';   
				push  @fasta_files   , $fasta_out;        ## from fasta 	             	 	     
				my $total_num  =  &fasta2fasta($file_in,$fasta_out);
				# load number
				$total_read_num{$fasta_out} = $total_num;
			}
			else{
				&single_fasta2array($file_in,$file)  ;
			}
	}
	
	if(@single_fas>0){
		my $fasta_out = 'Sanger_INPUT.fa';   
		push  @fasta_files   , $fasta_out;      ## from merged fasta 	
		open SINGLE_FAS, ">$fasta_out" || die " can not write single fasta \n";
		my $num = @single_fas+0;
		
		for (@single_fas){
			my ($seq_id, $seq_in) = @{$_};
			print SINGLE_FAS  ">",$seq_id,"\n",$seq_in,"\n";
		}
		close(SINGLE_FAS) ;
		$total_read_num{$fasta_out} = $num;
	}
	
	
	########################################################	  
	## Here will put step3 to step 5 into one loop
	
	
	open (REPORT,">$final_out") || die "cannot write\n";
	my @files_fas = @fasta_files;
	
	foreach my $file_blat_in (@files_fas){	
		############  step 3 - blat  ###########################
		my $blat_out        =  $file_blat_in.'_blat_crispr.psl';	 
		my $command_line =  'blat_exe/blat/'.$blat.' '.' -tileSize=7 -oneOff=1  -maxGap=20 -minIdentity=70 -minScore=20 '.$file_blat_in.' '.  $design_fas.' '.$blat_out;
		
		system("$command_line");
	
		print "blat job $file_blat_in is done\n"; 
		##########  step 4 - convert psl -> bed  ##############   
		
		my $bed_hash_address = &psl2bed ($blat_out );	## address of one hash
		
		##########  step 5 - get sequences number from bed ##############  
		
		my $read_num_address = &MAIN ($file_blat_in,$bed_hash_address);## address of one hash
		
	
		## remove files
		
		if($remove_files == 1){
				$command_line = $remove.' '.$file_blat_in;
				system("$command_line");
	
				$command_line = $remove.' '.$blat_out ;
				system("$command_line");
		}
	}
	
	
	if($remove_files == 1){
		my $command_line2 = $remove.' '.$design_fas ;
		system("$command_line2");
	}
	
	
	#########################################################################################
	## functions
	
		
	sub MAIN{
		my $fas_file_in         = shift @_; 
		my $bed_file_address_in = shift @_; 
		
		my %bed_hash_in = %{$bed_file_address_in};
		
	
		#########################################################
		# read convereted reads fasta file
			
			open (FAS,"$fas_file_in ")||die "can not open fasta file: $fas_file_in ";
			
			my %query_with_seq2num  = ();
			my %seq2alignment = ();
			my %seq_part2assignment = ();
			
			local $/ = "\n"; # correct
			
			while (my $lineIn =<FAS>){
				chomp $lineIn;
				if($lineIn =~ m/^\>/){
						my $id = substr $lineIn,1;
						
						my @numbers = split(/\_/, $id);
						my $total_num = $numbers[-1];
						
						$lineIn = <FAS>;
						chomp $lineIn;
						my $seq = $lineIn;
										
						if(exists $bed_hash_in{$id}){
								# $psl_hash{$read_target}{$query}  # format of hash
								## asign read to target reference 
							
							##########################################################
							my @read_asignment = ();
							## loop for best
							for my $read_key (keys %{$bed_hash_in{$id}}){
									### get annotation from array
									## Q: design T: read
									my @get_array = @{$bed_hash_in{$id}{$read_key}};
									my ($strand, $q_name,$q_size,$q_start,$q_end,
													$t_name,$t_size,$t_start,$t_end,
													$block_num,$block_size,$q_start_s,$t_start_s) =  @get_array ;
													
									my $seq_out = substr $seq,$t_start,($t_end-$t_start);
										
									if($strand eq '-'){										
										$seq_out = &rev_com ($seq_out);
									}
												## get orignal
										my $ref_ori_seq = '';
										if(exists $design_hash{$q_name}){
													$ref_ori_seq = $design_hash{$q_name};
										}
										
																					
										my $edit_note = &get_editing_note(\@get_array,$ref_ori_seq, $seq_out,$seq);			
																					# target   # part of read
										my ($align_ref,$align_read,$note) = split(/\t/, $edit_note);
	
										my $len_ori = length($ref_ori_seq);
										my @align_1 = split(//, $align_ref);
										my @align_2 = split(//, $align_read);
										my $site_snp = 0;				
										my @site_snp = ();
										my $iden_num = 0;
										
										for my $a1(@align_1){
											my $a2 = shift @align_2;
											if($a1 ne '-' and $a2 ne '-'){
												if($a1 ne $a2){
													push @site_snp, $site_snp;
												}
											}
											
											if($a1 eq $a2){
												$iden_num++;
											}
											
											$site_snp++;
										}
										$iden_num = $iden_num-2;
										
										## filter on SNP number
										my $snp_num_for_filter = @site_snp+0;
										
										my $dis_max = 0;
										my $dis_current = 0;
										my @continuous_site_temp = ();
										my @continuous_site = ();
										my @site_snp_checking = @site_snp;
										my $start = shift @site_snp_checking ;
										for my $site_checking (@site_snp_checking){
											my $dis_in = $site_checking -$start;
											if($dis_in <4){
												$dis_current++;
												push  @continuous_site_temp,$site_checking ;
												if($dis_current > $dis_max){
													$dis_max = $dis_current;
													@continuous_site = @continuous_site_temp;
												}
											}
											else{
												$dis_current = 0;
												@continuous_site_temp = ();
												
											}
	
											$start = $site_checking;
										}
										
					
										if( ($snp_num_for_filter/$len_ori) > 0.5){ # first filtering
											if($dis_max <5){
												next;  
											}
										}
										
										
										if($dis_max>=5){
											@site_snp = ();
											my @align_x = split(//, $align_read);
											for my $iii(@continuous_site ){
												$align_x[$iii] = 'N';
											}
											$align_read = join('',@align_x);
											$note = 'strange editing';
											my @ttt = ($align_ref,$align_read,$note) ;				                     	  
											
											$edit_note = join("\t", @ttt );
										}
										
										
										
										
										my @ooo = ($iden_num,$ref_ori_seq,$seq_out,$edit_note,\@site_snp,\@get_array);
										push @read_asignment,[@ooo];
							}## for loop for best hit of each fasta reads
							
							## get the best hits
							
							if(@read_asignment<1){
								next;
							}
					
							##  assign the first hit
							@read_asignment = sort {$b->[0] <=> $a->[0]}@read_asignment ;
							my ($iden_num,$ref_ori_seq,$seq_out,$edit_note,$array_snp_addr,$bed_array_addr) = @{$read_asignment[0]};
							
							my @snp_filtering = @{$array_snp_addr};				        
							
							if((@snp_filtering+0)/length($ref_ori_seq) >$mismatch_cutoff){
									next;
							}
							
							
							## after asignment 
							
							## check whether seq_out is assigned earlier
							## $seq_out has been asigned to one record
							## just assign reads to that one althogh it may be get one new assignment
							if( exists $seq_part2assignment{$seq_out}){
								my $first_assignment = $seq_part2assignment{$seq_out};
								$query_with_seq2num{$first_assignment}{$seq_out} = $query_with_seq2num{$first_assignment}{$seq_out}+$total_num;
								next;
							}
							
									
											# number for part of reads
							my ($strand, $q_name,$q_size,$q_start,$q_end,
													$t_name,$t_size,$t_start,$t_end,
													$block_num,$block_size,$q_start_s,$t_start_s) = @{$bed_array_addr};        
		
							# this the first asignment	
							if(not exists $seq2alignment{$seq_out}){
									my @out = ($ref_ori_seq,$edit_note,$array_snp_addr);
									$seq2alignment{$seq_out} = \@out;	
						
							}
										# assigned
							if(not exists $seq_part2assignment{$seq_out})	{
								$seq_part2assignment{$seq_out} = $q_name;
								$query_with_seq2num{$q_name}{$seq_out} = $total_num;
							}	
							
														
									
						}# exists  bed file
				}# lineIN
			}# while file
			close(FAS);
			
			### summary data for output
			my ($hash_addr1, $hash_addr2) = &get_out_buffer($fas_file_in,\%query_with_seq2num, \%seq2alignment);
			
			### write output
			
			&write_buffer_out($hash_addr1, $hash_addr2);
			
	} # Main function
			
		
		
		
	sub get_out_buffer{
		my $fas_file_in = $_[0];
		my %query_with_seq2num_in =  %{$_[1]};
		my %seq2alignment_in       = %{$_[2]};
		
		my %buffer_out = ();
		my %buffer_num = ();
		
		my $total_non_redun = 0;
		## push output in buffer before printing 
		my $test_case_num = $total_read_num{$fas_file_in};
	
		for my $ref_name (sort keys %query_with_seq2num_in){
			my %in_hash = %{$query_with_seq2num_in{$ref_name}};
			
			## each group
			my $report_wt_count = 0;
			my $report_indel_count = 0;
			my $other_num = 0;
			my $sub_hit_num = 0;
			my $indel_hit_num = 0;
			
			
			my @other_record = ();
			### loop for records
			for my $part_of_read (sort {$in_hash{$b} <=> $in_hash{$a}} keys %in_hash){
				
					my $hitnum = $in_hash{$part_of_read};
					
					if($hitnum < $min_cutoff) {
						next;
					}
					
					if($hitnum/$test_case_num < 0.00005) {
						next;
					}
			
					## get Editing Note for output only
							
				if(exists $seq2alignment_in{$part_of_read}){
					my ($ref_ori_seq,$edit_note,$array_snp_addr)= @{$seq2alignment_in{$part_of_read}};
					my ($align_ref,$align_read,$note) = split(/\t/,$edit_note);
						my @snp_out = @{$array_snp_addr};
						my $snp_num = @snp_out +0;
						# deep sequencing SNP
						if($test_case_num >500 and $hitnum <=2 and $snp_num/length($ref_ori_seq)>0.05){
							next;
						}
						
						my $snp_note = '';
						if($snp_num>0 ){
								$snp_note = $snp_num.' SNP('.join(',',@snp_out).')';
						}
						
						$sub_hit_num = $sub_hit_num+$hitnum;
						$total_non_redun = $total_non_redun + $hitnum;
						
						
						# output here
						if($ref_ori_seq ne $design_hash{$ref_name} ){
							print "error in read assingment \n";
						}
		
						
					my @out_line = ($fas_file_in,$ref_name,$ref_ori_seq, $part_of_read,$hitnum, $align_ref,$align_read,$note,$snp_note);
											## indel 			 	 
					if($edit_note !~ m/no.+indel/ ){
							
						if((($test_case_num >500 and $hitnum <=3) or ($hitnum /$test_case_num <0.001 ) ) and $edit_note =~ m/strange/){
							next;   # filter strange case in deep sequencing with few reads
						}
	
						
						$report_indel_count ++;
						# total indel
						$indel_hit_num = $indel_hit_num + $hitnum;
						
						if($report_indel_count <= $indel_report ){
								#print REPORT join("\t",@out_line),"\n";	 
								## in buffer
								if(not exists $buffer_out{$fas_file_in}{$ref_name}{"data"}){
									my @bbb = ();
									push @bbb, [@out_line];
									$buffer_out{$fas_file_in}{$ref_name}{"data"} = \@bbb;
								}
								else{
									push @{$buffer_out{$fas_file_in}{$ref_name}{"data"}},[@out_line];
								}
						}
						else{
								$other_num = $other_num+$hitnum;
								@other_record = ($fas_file_in,$ref_name,$ref_ori_seq, 'others',$hitnum, '-', '-', '-', '-');	
						}
	
					}
					else{
						$report_wt_count ++;
						if($report_wt_count <= $wt_like_report  ){
							#print REPORT join("\t",@out_line),"\n";	 
							# in buffer
							if(not exists $buffer_out{$fas_file_in}{$ref_name}{"data"}){
									my @bbb = ();
									push @bbb, [@out_line];
									$buffer_out{$fas_file_in}{$ref_name}{"data"} = \@bbb;
								}
								else{
									push @{$buffer_out{$fas_file_in}{$ref_name}{"data"}},[@out_line];
								}
						}
						else{
							$other_num = $other_num+$hitnum;
							@other_record = ($fas_file_in,$ref_name,$ref_ori_seq, 'others',$hitnum, '-', '-', '-', '-');	 	
						}
						
					}# wt like report	 
					
				} # if exists 
	
			}# for part of reads
			
			## other report
			if($other_num>0){
						$other_record[-5]	 = $other_num;
						#print REPORT join("\t",@other_record),"\n";
						
						if(not exists $buffer_out{$fas_file_in}{$ref_name}{"data"}){
									my @bbb = ();
									push @bbb, [@other_record];
									$buffer_out{$fas_file_in}{$ref_name}{"data"} = \@bbb;
						}
						else{
									push @{$buffer_out{$fas_file_in}{$ref_name}{"data"}},[@other_record];
						}
			}
			
				#my @summary_record = ($fas_file_in,$ref_name,$ref_ori_seq, 'total_hit:'.$total_non_redun,'sub_hit:'.$sub_hit_num, 'indel_hit:'.$indel_hit_num);
					## total read number
			my $total_num =  $total_read_num{$fas_file_in};
			
								#    0           1          2                           
			# my @summary_record = ($fas_file_in,$ref_name,'Total Reads: '., 'Total Hits: '.$total_non_redun,'Sub Hits: '.$sub_hit_num, 'Indel Hits: '.$indel_hit_num);	 
				my @summary_record = ('Sub Hits: '.$sub_hit_num, 'Indel Hits: '.$indel_hit_num);	 
				
				$buffer_out{$fas_file_in}{$ref_name}{'sum'}   =   \@summary_record ;
				$buffer_num{$fas_file_in}{'total'} =   $total_num ;
				$buffer_num{$fas_file_in}{'sub'}   =   $total_non_redun ;
				#print REPORT join("\t",@summary_record),"\n";
				
		}# for ref_seq
		
		
		return (\%buffer_out,\%buffer_num);
	} # sub 
	
	
	### write the output
	
	sub write_buffer_out{
		my %hash_out = %{$_[0]};
		my %hash_out_num = %{$_[1]};
		my $print_tracking = 0;
		for my $fas_file_in (sort keys %hash_out){
			
			my  $total_num       = $hash_out_num{$fas_file_in}{'total'} ;
			my  $total_non_redun = $hash_out_num{$fas_file_in}{'sub'} ;	
				
			for my $ref_name (sort keys %{$hash_out{$fas_file_in}}){
				if(not exists $hash_out{$fas_file_in}{$ref_name}{"data"}){
					# print "No data for $fas_file_in $ref_name \n";
					next;
				}
				my @data = @{$hash_out{$fas_file_in}{$ref_name}{"data"}};
				
				for (@data){
					$print_tracking++;
					if($print_tracking==1){
						print REPORT "INPUT","\t","Target\t","TargetSequence\t","ReadSequence\t","Read#","\t",
										"AlignedTarget\t","AlignedRead\t", "Indels\t","SNPs\n";
					}
					print REPORT join("\t", @{$_}),"\n"; 
				}
				
				# sum
				my  @sum_p2 = @{$hash_out{$fas_file_in}{$ref_name}{'sum'}} ;
				
				
				my  @sum_p1 = ($fas_file_in,$ref_name,'Total Reads: '.$total_num, 'Total Hits: '.$total_non_redun);
				
				print REPORT join("\t", @sum_p1),"\t", join("\t", @sum_p2),"\n";
				
			}
			
			
			##### print summary region
		print REPORT "Summary \t -----------------------------------------------------\n";
		print REPORT "Sum:INPUT\t","Target\t","AlignedTarget\t","AlignedRead\t","Total Hits\t","Sub Hits\t","Indel or WT Hits\t","Indel or WT rate %\t","Pattern\n";
			
			
			
			for my $ref_name (sort keys %{$hash_out{$fas_file_in }}){
				if(not exists $hash_out{$fas_file_in}{$ref_name}{"data"}){
					# print "No data for $fas_file_in $ref_name \n";
					next;
				}
				
				my @data = @{$hash_out{$fas_file_in}{$ref_name}{"data"}};
				
				my $wt_pair    = ''."\t".'';
				my $indel_pair = ''."\t".'';
				my $indel_out  = '';
				my $i = 0;
				
				for (@data){
					$i++;
					my @in= @{$_};
					
					
					if($in[7] !~ m/no\_indel/ and $in[7] ne '-'){
							if($indel_pair eq ''."\t".''){
									$indel_pair = $in[5]."\t".$in[6];
									my @indel_temp = split(/\s+/,$in[7]);
									
									my @ooo_indel = ();
									for my $ooo_in (@indel_temp ){
										if($ooo_in =~ m/I(\d+)/){
											push @ooo_indel, '+'.$1;		 	 	  	   	   
										}
										if($ooo_in =~ m/D(\d+)/){
											push @ooo_indel, '-'.$1;		 	  	   	   
										}
										if($ooo_in =~ m/strange/){
											push @ooo_indel, 'strange editing';		 	  	   	   
										}
									}
								$indel_out = 'Indel:'.join(',', @ooo_indel);
							}
					}
					else{
						# no indel
						if($wt_pair eq ''."\t".''){
								$wt_pair = $in[5]."\t".$in[6];
						}
					}
				}
				
				# first record
				
			
	
				# sum
				my  ($sub_string, $indel_string) = @{$hash_out{$fas_file_in}{$ref_name}{'sum'}} ;
				my $sub_num ;
				my $indel_num;	
				#print 		$sub_string. "sub summary\n";
				#print 		$indel_string ."sub summary\n";
				
				
				if($sub_string =~ m/\:[^\d]+(\d+)/){
						$sub_num = $1;
				}
				if($indel_string =~ m/\:[^\d]+(\d+)/){
						$indel_num = $1;
				}
				
				my $rate = int(($indel_num/$sub_num)*10000)/100;
				
				## indel 
													# two columns
				my @out_sum = ($fas_file_in,$ref_name, $indel_pair,$total_non_redun,$sub_num,$indel_num,$rate,$indel_out);
				print REPORT 'Sum:',join("\t", @out_sum),"\n";
				## wt like
				my $wt_num = $sub_num - $indel_num;
				my $wt_rate = 100-$rate;
				my @out_sum = ($fas_file_in,$ref_name, $wt_pair  ,$total_non_redun,$sub_num,$wt_num,$wt_rate,"WT_like");
				print REPORT 'Sum:',join("\t", @out_sum),"\n";
				
			} # records
			
			print REPORT "\n";
			
		} # file 
	} # sub 
	###################################################################
	
	
	
	sub get_editing_note{
		my ($address_in,$query_seq, $target_seq ,$full_read) = @_;
		# target: part of read
		# query : original sequence in design file
		
		my ($strand, $q_name,$q_size,$q_start,$q_end,
					$t_name,$t_size,$t_start,$t_end,
					$block_num,$block_size,$q_start_s,$t_start_s) =  @{$address_in};	
					
		my @out = ();	
		#print 		$strand ,"strand\n";
		if($strand eq '+'){
			@out = get_alignment($block_num,$block_size,$q_start_s,$t_start_s,$query_seq, $full_read);
		}
		else{
			@out = get_alignment($block_num,$block_size,$q_start_s,$t_start_s,rev_com($query_seq), $full_read);
			@out = get_reverse(@out);
		}
		
		$out[0]='B'.$out[0].'E';
		$out[1]='B'.$out[1].'E';
		if($out[2] !~ m/\w/){
			$out[2] = 'no_indel';
		}
		return join("\t",@out);
	}
	
	
	
	sub get_alignment{
		my ($blocks ,  $blockLengths,   $qStarts ,   $tStarts,$query_seq, $target_seq) =  @_;
		my @blockSizes = split(/\,/,$blockLengths);
		my @qArray = split(/\,/,$qStarts);
		my @tArray = split(/\,/,$tStarts);
		
		# target: part of read
		# query : original sequence in design file
		
		# print join("\t",@_),"\n";
		
		### before blocks     
		my $before_q = substr $query_seq, 0,$qArray[0];
		my $before_t = substr $target_seq,0,$tArray[0];
		
	
		($before_q,$before_t) = treat_sequence($before_q,$before_t,"before");
		
		### after blocks     
		my $after_q = substr $query_seq, $qArray[-1]+$blockSizes[-1];
		my $after_t = substr $target_seq,$tArray[-1]+$blockSizes[-1];
		
	
		($after_q,$after_t) =  treat_sequence($after_q,$after_t,"after") ;
	
		### blocks
		my @med_q = ();
		my @med_t = ();
		
		my @indel_out = ();
		my $out = '';
		## first block
			my $med_q_seq = substr $query_seq,$qArray[0],$blockSizes[0];
			my $med_t_seq = substr $target_seq,$tArray[0],$blockSizes[0];
			push @med_q,$med_q_seq;
			push @med_t,$med_t_seq;
	
		if($blocks>1){
			for (my $i= 0;$i<($blocks-1);$i++){
				####### interval 
				my $inter_q_seq = substr $query_seq,($qArray[$i]+$blockSizes[$i]), ($qArray[$i+1]-($qArray[$i]+$blockSizes[$i]));
				my $inter_t_seq = substr $target_seq,($tArray[$i]+$blockSizes[$i]),($tArray[$i+1]-($tArray[$i]+$blockSizes[$i]));
				
				### count deletion insertion
				if(length($inter_t_seq) - length($inter_q_seq)>0){
					# to fix the mismatch near deletion   length($inter_q_seq)
					$out = ($qArray[$i]+$blockSizes[$i]+length($inter_q_seq)).'I'.(length($inter_t_seq) - length($inter_q_seq));
				}
				if(length($inter_q_seq)-length($inter_t_seq)>0){
					$out = ($qArray[$i]+$blockSizes[$i]+length($inter_t_seq)).'D'.(length($inter_q_seq)-length($inter_t_seq));
				}     	      
				push @indel_out,$out;     	      
				
				($inter_q_seq,$inter_t_seq) = treat_inter ($inter_q_seq,$inter_t_seq) ;
				push @med_q, $inter_q_seq;
				push @med_t, $inter_t_seq;  
				
				######  block   after interval  
				$med_q_seq = substr $query_seq,$qArray[$i+1],$blockSizes[$i+1];
				$med_t_seq = substr $target_seq,$tArray[$i+1],$blockSizes[$i+1];
				push @med_q,$med_q_seq;
				push @med_t,$med_t_seq;     	 	
			}
		}
		
				push @med_q,$after_q;
				push @med_t,$after_t;
				
				unshift  @med_q,$before_q;
				unshift  @med_t,$before_t;
				my $q_string = join('',@med_q);
				my $t_string = join('',@med_t);
				
				if($q_string=~m/^(\-+)/){
					my $len = length($1);
					$q_string = substr $q_string,$len;
					$t_string = substr $t_string, $len;
				}
				if($q_string=~m/(\-+)$/){
					my $len = length($1);
					$q_string = substr $q_string,0,(-1)*$len;
					$t_string = substr $t_string,0, (-1)*$len;
				}
				
	
				return ($q_string,$t_string,join(' ',@indel_out));
	
	}
	
	
	sub get_reverse{
			my ($q_string,$t_string,$indel_string )= @_;
			$q_string = rev_com($q_string);
			$t_string = rev_com($t_string);
			
			my $string_test = $q_string;
			$string_test =~ s/\-//g;
			my $len = length($string_test);
			
			my @indel_in =  split(/\s/, $indel_string);
			my @indel_out = ();
			for (@indel_in){
				
				if(m/^(\d+)(.)(\d+)/){
					my $ooo = ($len-$1-$3).$2.$3;
					if($2 eq 'I'){
						$ooo = ($len-$1).$2.$3;
					}
					
					unshift @indel_out,$ooo;
				}
			}
			$indel_string = join(' ',@indel_out);
			return($q_string,$t_string, $indel_string);
	}
	
	
			
		
	
	sub treat_inter{
		my ($q,$t) = @_;
		my $q_len = length($q);
		my $t_len = length($t);
		my $dis = abs($q_len-$t_len);
		my @oooo = ();
		
		for(1..$dis){
				push @oooo,'-';
		}
		
		if($q_len-$t_len >0){
			$t = $t.join('',@oooo);
		}
		else{
			$q = $q.join('',@oooo);
		}
		return ($q,$t);
	}
	
	
	
	sub treat_sequence{
		my ($q,$t,$position) = @_;
		my $q_len= length($q);
		my $t_len = length($t);
		my $dis = abs($q_len-$t_len);
		
		if($dis > 0){
			my @oooo = ();
			for(1..$dis){
				push @oooo,'-';
			}
			if($q_len>$t_len){
				if($position eq 'before'){
					$t = join('',@oooo).$t ;
				}
				else{
					$t = $t.join('',@oooo) ;
				}
			}
			else{ # q small
				if($position eq 'before'){
					$q = join('',@oooo).$q ;
				}
				else{
					$q = $q.join('',@oooo) ;
				}
			}
		}
		return ($q,$t);
	}
	
	
	sub psl2bed  {
		my $psl_file = shift @_; # file in @psl_file_data
		my @psl_array = ();
		
		my $link_br = &check_line_break($psl_file);
		
		
		open (FILEINBLAT,"$psl_file") || die "Can not open PSL file, Please check BLAT software is working\n";
		
		## should be correct
		local $/ = $link_br;    
		
		for(1..5){
				my $line=<FILEINBLAT>;
				## remove head		   
		}
		## read psl
		while (<FILEINBLAT>)
		{
			my $lineIn = $_;
			chomp $lineIn;
			my @lineArr  = split ("\t", $lineIn);
			
			my $q_size   = $lineArr[10];
			my $m_size   = $lineArr[0];
			
			if( $m_size/$q_size > 0.20){
				push @psl_array,[@lineArr];
			}
		
		}# while file
		close(FILEINBLAT);
		## score large to small 
		
		#@psl_array = sort {$b->[0] <=>$a->[0]}@psl_array ;
		
		my @out_array = ();  
		my $psl_line_num = 0;
		my %psl_hash = ();   
		for  (@psl_array)
		{
				my @lineArr  = @{$_};
				$psl_line_num++;
				my $read_target  = $lineArr[13];
				my $query        = $lineArr[9];
				for(1..8){
					shift  @lineArr;
				}
	
				if(not exists $psl_hash{$read_target}{$query}){		## all alignment
				#if(not exists $psl_hash{$read_target}){  ## only keep the best hit			
					$psl_hash{$read_target}{$query} = \@lineArr;
				}
		}# while array
		print $psl_line_num." total psl lines\n";
		return (\%psl_hash);
	}
	
	
	
	
	sub single_fasta2array {
		my $input_file = shift @_;
			#@single_fas
		my $file_short = shift @_;
		
		my $link_br = &check_line_break($input_file);
		
		open (IN,"<$input_file" )or die "File $input_file  not found error here\n";
	
		my $id = $file_short;
		my $seq  = '';    
	
		local $/ = $link_br;
		while(<IN>){
			chomp;
			if($_ =~ m/^\>/ ){
				my @ttt = split(/\s+/,$_);
				$id = substr $ttt[0],1;
			}
			else{
				$seq = $seq.$_;
			}
		}
		# last records
		$seq =~ s/\s|\n|\r//g;
		
		$id =~ s/\_//g;
		$id = $id.'_1';
		
		push @single_fas,[($id,$seq)];
	
		close IN;
	}
	
	
	
	sub fasta2fasta {
		my $input_file = shift @_;
		my $out_file = shift @_;
		
		
		my $line_br =  &check_line_break($input_file);
		
		
		open IN,"<$input_file " or die "File $input_file  not found error here\n";
		open TGT,">$out_file" or die "File $out_file not found error here\n";
		my $c = 0;
		my @line = ();
		my $id ='';
		my $seq = '';
		my $processed = 0;
		
		local $/ = $line_br; 
		
		while(<IN>){
			chomp;
			if(m/^\>/){
				$processed ++;
				if($processed == 1){
					my @seq_id = split(/\s+/, $_);
					my $seq_id = $seq_id[0];
					$seq_id =~ s/\_//g;        		
					$seq_id = $seq_id.'_1';
					print TGT '>'.$seq_id,"\n";
				}
				else{
					$seq =~ s/\s|\n|\r//g;
					print TGT $seq,"\n";
					$seq = '';
					my @seq_id = split(/\s+/, $_);
					my $seq_id = $seq_id[0];
					$seq_id =~ s/\_//g;        		
					$seq_id = $seq_id.'_1';
					print TGT '>'.$seq_id,"\n";
				}
				
			}
			else{
				$seq = $seq.$_;
			}
		}
		# last records
		$seq =~ s/\s|\n|\r//g;
		print TGT $seq,"\n";
		
		close IN;
		close TGT;
		return $processed;
	}
		
		
		
	sub fastq2fasta {
		my $input_file = shift @_;
		my $out_file = shift @_;
		
		my $link_br = &check_line_break($input_file);
		open IN,"<$input_file " or die "File $input_file  not found error here\n";
		open TGT,">$out_file" or die "File $out_file not found error here\n";
		my $c = 0;
		my @line = ();
		my $id ='';
		my $seq = '';
		my $processed = 0;
		my %read_hash = ();
		
		##############################################################################################
		local $/ = $link_br;
		
		while(<IN>){
			chomp;
			$c++;
			if($c == 1){
				$processed++;
				#        print STDERR "$processed reads processed\r";
				@line = split();
				$id = $line[0];
				$id =~ s/\@//;
			}elsif($c == 2){
				$seq = $_;
			}elsif($c == 4){
				$c = 0;
				#print TGT ">$id\n$seq\n";
			# print TGT ">seq_$processed\n$seq\n";
				$read_hash{$seq}{$id} = 1;
				$id ='';
				@line =();
				$seq ='';
			}else{}
		}
		
		## write TGT
		my $count = 0;
		for my $seq_out (sort  keys %read_hash){
			$count++;
			my @hits = keys %{$read_hash{$seq_out}};
			my $num = @hits+0;
			my $id_out = "R".$count.'_'.$num;
			print TGT ">$id_out\n", $seq_out ,"\n";    	 
		}
		
		close IN;
		close TGT;
		return $processed;
	}
	
	
	sub check_read_num{
			my $file_in = shift @_;
			# check read number
			
			my $line_break = check_line_break($file_in);
			open FAS_in,"$file_in" || die "can not read fasta file\n";
			local $/ = $line_break;  
	
			my $seq_num = 0;
			while (<FAS_in>){
					if( m/^\>/){
						$seq_num ++;
					}
					if( $seq_num>1){
						last;
					}
			}
			
			return $seq_num;
		
	}
	
	sub check_line_break{
		my $file_in = shift @_;
		
		my $size ;
		if (-e $file_in){
			$size = -s $file_in;
		}
		if($size>1000){
			$size=1000;
		}
		
		open (FH,"<$file_in ")|| die "cannot open $file_in\n";
		my $buffer;
		my $value = read (FH, $buffer, $size, 0);
		my $out = "";
		if($buffer =~ m/\r\n/){
			$out = "\r\n";
		}
		else{
			if($buffer =~ m/\r/){
				$out = "\r";
			}
			if($buffer =~ m/\n/){
				$out = "\n";
			}
		}
		
		close(FH);
		return ($out);
	}
	
	
	sub rev_com {
	my $seq = shift @_;
	$seq = reverse($seq);
	$seq =~ tr/ACGT/TGCA/;
	return($seq);
	}
	
}

	

	sub handle_request {
		
		my ($self, $cgi) = @_;
		my $path_info = $cgi->path_info;
	
		if ($path_info eq '/ageseq') {
			my $directory_name = $cgi->param('directory');
			my $mismatch_cutoff = $cgi->param('mismatch_cutoff');
			my $min_cutoff = $cgi->param('min_cutoff');
			my $wt_like_report = $cgi->param('wt_like_report');
			my $indel_report = $cgi->param('indel_report');
			
			my $directory_reads = "/data/AGEseq/$directory_name/reads";
			my $targetsFilepath = "/data/AGEseq/$directory_name/targets.txt";  
			my $outputFilepath = "/data/AGEseq/$directory_name/output.txt";  # Construct the file path
		
			ageSeq($directory_reads, $targetsFilepath, $outputFilepath, $mismatch_cutoff, $min_cutoff, $wt_like_report, $indel_report);
			
			
			#print "HTTP/1.0 200 OK\r\n";
			#print "Content-Type: application/json\r\n\r\n";
			print '{"message":"Successful"}';

			
		}
		else{
			
			#my ($self, $cgi) = @_;
		
			# Get the value of the 'name' parameter from the request
			my $post_data = $cgi->param('POSTDATA');
			my $data = JSON::XS::decode_json($post_data);
		
			# Extract the value of the 'name' field
			my $alldata = $data->{data};
			my $timestamp = Time::Piece->localtime->strftime('%Y%m%d%H%M%S');# . "_reads";
			
			my $directory = "/data/AGEseq/$timestamp";
			mkdir $directory or die "Failed to create directory: $!";
			
			my $directory_reads = "/data/AGEseq/$timestamp/reads";
			mkdir $directory_reads or die "Failed to create directory: $!";
		
		
			# Send the response
			print "HTTP/1.0 200 OK\r\n";
			print "Content-Type: text/plain\r\n\r\n";
			
			# Iterate over each key-value pair
			foreach my $filename (keys %$alldata) {
				# Print the filename
				print "Filename: $filename\n";
				
				open(my $fh, '>', "$directory_reads/$filename") or die "Could not open file '$filename' $!";
				# Iterate over the lines in the array
				foreach my $line (@{$alldata->{$filename}}) {
					# Print each line
					print $fh "$line\n";
				}
				print "\n"; # Add a new line between files
		
		
			}
				
			my $targetData = $data->{targets};
			#my $targetsFilepath = "/data/AGEseq/$timestamp"."targets.txt";  # Construct the file path
			my $targetsFilepath = "/data/AGEseq/$timestamp/targets.txt";  # Construct the file path
			#my $targetsFilepath = "targets.txt";  # Construct the file path
			
			print $targetsFilepath;
			open(my $tr, '>', $targetsFilepath) or die "Could not open file '$targetsFilepath' $!";
			foreach my $line (@{$data->{'targets'}}) {
				print $tr "$line\n";
				
			}
			
			my $outputFilepath = "/data/AGEseq/$timestamp/output.txt";  # Construct the file path
		
			ageSeq($directory_reads, "targets.txt", $outputFilepath);
			#ageSeq($directory_reads, $targetsFilepath, $outputFilepath);
			
			
			
			
			#print "Content-Type: application/json\r\n\r\n";
			print "$directory\n";
			print "$targetsFilepath\n";
			print "$outputFilepath\n";
			print "Request URI: ", $cgi->request_uri, "\n";
			#print { message => "Successful!" };
			
		}
			
	}
	
	
}

# Create an instance of your subclassed server
my $server = MyHTTPServer->new($port);

# Start the server
print "Server running at http://$hostname:$port/\n";
$server->run($hostname, $port);






















##!/usr/bin/perl
#
#use strict;
#use warnings;
#use HTTP::Server::Simple::CGI;
##use CGI;
##use JSON::XS qw(decode_json);
#
#my $port = 8502;
#my $hostname = "172.30.18.104";
#
#{
#    package MyHTTPServer;
#    use base qw(HTTP::Server::Simple::CGI);
#	use CGI;
#
#    sub handle_request {
#		print "---------";
#        my ($self, $cgi) = @_;
#		# Debugging output: Print request method
#		#print "Request Method: ", $cgi->request_method, "\n";
#
#		# Debugging output: Print request URI
#		#print "Request URI: ", $cgi->request_uri, "\n";
#
#        ## Get the JSON data from the request
#        #my $json_data = $cgi->param('json_data') // '{}';
#		#
#		#
#		#print $json_data;
#		#
#        ## Decode the JSON data into a Perl data structure
#		#my $data = JSON::XS::decode_json($json_data);
#        ##my $data = decode_json($json_data);
#		#
#        ## Extract the value of the 'name' field
#        #my $name = $data->{name};
#		#
#        ## Generate the response
#        my $response = "Hello";#, $name";
#
#        # Send the response
#        print "HTTP/1.0 200 OK\r\n";
#        print "Content-Type: text/plain\r\n\r\n";
#        print $response;
#    }
#}
#
## Create an instance of your subclassed server
#my $server = MyHTTPServer->new($port);
#
## Start the server
#print "Server running at http://$hostname:$port/\n";
#$server->run($hostname, $port);
