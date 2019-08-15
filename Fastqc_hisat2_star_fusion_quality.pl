#!/usr/bin/perl -w
use strict;
use warnings;
use POSIX qw(strftime);
use Getopt::Std;
use Getopt::Long;
use File::Basename;
use Cwd 'abs_path';
my $verbose ="v1.0";

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my ($fq1,$fq2,$outdir,$hg38,$onlyprintcmd,$gatktype,$bed,$trim_q,$REMOVE_DUPLICATES);
my %opts;
GetOptions(
    'verbose' => \$verbose,
    'fq1=s' => \$fq1,
	'fq2=s' => \$fq2,
    'outdir=s' => \$outdir, 
	'hg38' => \$hg38,
#	'gatktype=i' => \$gatktype,
    'onlyprintcmd' => \$onlyprintcmd,
 	'bed=s' => \$bed,  
 	'trim_q=i'=>\$trim_q, 
 	'REMOVE_DUPLICATES'=>\$REMOVE_DUPLICATES,
#    'libs=i' => \@libs, ## 'libs=i@' => \$libs,
#    'define=s' => \%defines, ## 'define=s%' => \$defines,
) or die $!;
#&usage unless ( exists $fq1 && exists $outdir );
unless(defined $fq2 && defined $fq1 ){&usage();exit 0;}
my $fastqc="/home/fuzl/soft/FastQC/fastqc";
$bed||="/home/fuzl/bed/2.probe_covered_regions.bed";
$trim_q||="13";
my $qualimap="/home/fuzl/soft/qualimap_v2.2.1/qualimap" ;
my $hg19_root="/home/fuzl/soft/GATK/resources/bundle/hg19";

my ($genome,$INDEX,$dbSNP,$phasel_1KG,$Mills_and_1KG,$phasel_snp_1KG,$Fusion_database);
if (defined $hg38) { #dbSNP 库还有问题
	my $hg38_root="/home/fuzl/soft/GATK/resources/bundle/hg38/";
	$genome="$hg38_root/Homo_sapiens_assembly38.fasta";
	$INDEX="$hg38_root/bwa_index/gatk_hg38";
	$dbSNP="$hg38_root/dbsnp_146.hg38.vcf.gz";
	$phasel_1KG="$hg38_root/1000G_phase1.snps.high_confidence.hg38.vcf";
	$Mills_and_1KG="$hg38_root/Mills_and_1000G_gold_standard.indels.hg38.vcf";
}else {
	$genome="$hg19_root/ucsc.hg19.fasta";
	$INDEX="$hg19_root/bwa_index/gatk_hg19";
	$dbSNP="$hg19_root/dbsnp_138.hg19.vcf";	
	$phasel_snp_1KG="$hg19_root/1000G_phase1.snps.high_confidence.hg19.sites.vcf";
	$phasel_1KG="$hg19_root/1000G_phase1.indels.hg19.sites.vcf";
	$Mills_and_1KG="$hg19_root/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf";
	$Fusion_database="/home/fuzl/soft/Fusion/GRCh37_v19_CTAT_lib_Feb092018/ctat_genome_lib_build_dir";

}

my $picard="/home/fuzl/soft/picard.jar";
my $bin_trim_galore ="/home/fuzl/soft/TrimGalore-0.5.0/trim_galore";
#my $samtools="/home/wangb/samtools";
my $samtools="/home/fuzl/soft/samtools-1.9/samtools";
#my $GATK="/home/fuzl/soft/gatk-4.0.8.1/gatk";
my $GATK="/home/fuzl/soft/GATK/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar";
#my $bwa="/Share/home/tiangeng/.Work1/PERL/WGS/bwa-0.7.12/bwa";
my $bwa="/home/fuzl/soft/bwa-0.7.17/bwa";
my $STAR_Fusion="/home/fuzl/soft/Fusion/STAR-Fusion-v1.5.0/STAR-Fusion"; 
my $hisat2="/home/fuzl/soft/hisat2-2.1.0/hisat2";
$outdir||='.';
&MKDIR("$outdir");
$outdir=abs_path("$outdir");
$fq1=abs_path($fq1);
$fq2=abs_path($fq2);

my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is: \nfq1:$fq1\nfq2:$fq2\nOutput directory is $outdir\n";
print "Database file is $genome\n";

###############################################################################
my $cmd="";

my $genomename=basename($genome);
$genome=abs_path($genome);
my $genome_m=$genome;
$genome_m=~s/\.fa(sta)?$//;#绝对路径
my $genomename_m=$genomename;
$genomename_m=~s/\.fa(sta)?$//;#基本名

unless (-f "$genome.fai" && -f "$genome_m.dict" && -f "$INDEX.bwt")  {
#	system "mkdir $outdir/genome/ && ln -s $genome $outdir/genome/";
	&MKDIR("$outdir/genome/");
	`rm $outdir/genome/* && ln -s  $genome $outdir/genome/`;
	$genome="$outdir/genome/$genomename ";
	$cmd .="cd $outdir/genome/ \n";
	$cmd .="$samtools faidx $genome  \n" ;
	$cmd .="java -Xmx20G -jar $picard CreateSequenceDictionary R=$genome O=$genomename_m.dict \n" ;
	$cmd .="bwa index -a bwtsw -p $genomename_m $genome \n";
	$INDEX ="$outdir/genome/$genomename_m";
	&runcmd("Building database",$cmd);
}

#############################################
#QC
$cmd="";

my $fq_name=basename($fq1);
#$fq_name=~s/_(R)?1(.*)?.f(ast)?q(\.gz)?$//;
$fq_name=~s/_(R)?1(.*)?\.f(ast)?q(\.gz)?$//;
#$cmd .="mkdir $outdir/FASTQC \n" unless (-d "$outdir/FASTQC");
&MKDIR("$outdir/FASTQC");
$cmd .="$fastqc $fq1 -o $outdir/FASTQC/  && cd $outdir/FASTQC/ &&  unzip -o $outdir/FASTQC/${fq_name}_1_fastqc.zip \n";
$cmd .="$fastqc $fq2 -o $outdir/FASTQC/  && cd $outdir/FASTQC/ &&  unzip -o $outdir/FASTQC/${fq_name}_2_fastqc.zip \n ";

#$cmd .="mkdir $outdir/trim/ \n" unless (-d "$outdir/trim")
#$cmd .="$bin_trim_galore -q 25 --phred33 --length 50 -e 0.1 --stringency 3 --paired -o $outdir/trim  $fq1 $fq2 \n";

&MKDIR("$outdir/trim/");
#$cmd .="mkdir $outdir/trim/ \n" unless (-d "$outdir/trim");
$cmd .="java -jar /home/fuzl/soft/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 12 -phred33 -trimlog $outdir/trim/logfile ";
$cmd .="$fq1 $fq2 $outdir/trim/${fq_name}_1_paired.fq.gz  $outdir/trim/${fq_name}_1_unpaired.fq.gz  $outdir/trim/${fq_name}_2_paired.fq.gz $outdir/trim/${fq_name}_2_unpaired.fq.gz ";
#$cmd .="ILLUMINACLIP:/home/fuzl/soft/Trimmomatic-0.38/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 \n";
$cmd .=" HEADCROP:5 ";
#$cmd .=" CROP:147 "; # reads length 151 bp，切除3‘端3bp
$cmd .="ILLUMINACLIP:/home/fuzl/soft/Trimmomatic-0.38/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:$trim_q MINLEN:36 \n";# 修改LEADING 和 SLIDINGWINDOW

$fq1="$outdir/trim/${fq_name}_1_paired.fq.gz";
$fq2="$outdir/trim/${fq_name}_2_paired.fq.gz";
$fq_name=basename($fq1);
$fq_name=~s/_1_paired.f(ast)?q(\.gz)?$//;
&MKDIR("$outdir/FASTQC_trim/");
#$cmd .="mkdir $outdir/FASTQC_trim \n" unless (-d "$outdir/FASTQC_trim");
$cmd .="$fastqc $fq1 -o $outdir/FASTQC_trim/  && cd $outdir/FASTQC_trim/ && unzip -o $outdir/FASTQC_trim/${fq_name}_1_paired_fastqc.zip  & \n";
$cmd .="$fastqc $fq2 -o $outdir/FASTQC_trim/  && cd $outdir/FASTQC_trim/ && unzip -o $outdir/FASTQC_trim/${fq_name}_2_paired_fastqc.zip  & \n ";

&qsub("QC",$cmd,"10G");


#############################################################
#bwa
my $sample=$fq_name;
#`mkdir $outdir/GATK` unless (-d "$outdir/GATK");
&MKDIR("$outdir/hisat2");

$cmd="";
#$cmd .="mkdir $outdir/bwa \n" unless (-d "$outdir/bwa");
#$cmd .="$bwa mem -M -t 36 -R '\@RG\\tID:$sample\\tSM:$sample\\tLB:$sample\\tPL:Illumina' $INDEX  $fq1 $fq2 >$outdir/bwa/$sample.sam  \n";
$cmd .="$hisat2 -p 12 --dta -x $hg19_root/hisat_index/hg19 -1 $fq1 -2 $fq2 -S $outdir/hisat2/$sample.sam 2>$outdir/hisat2/$sample.alnstats \n";
$cmd .="$samtools view -Sb $outdir/hisat2/$sample.sam | $samtools sort >$outdir/hisat2/$sample.sort.bam  \n";
&qsub("Hisat2",$cmd,"20G");
=c
   $HISAT2 -p $NUMCPUS --dta -x ${GENOMEIDX} \
     -1 ${FASTQLOC}/${reads1[$i]} \
     -2 ${FASTQLOC}/${reads2[$i]} \
     -S ${TEMPLOC}/${sample}.sam 2>${ALIGNLOC}/${sample}.alnstats
=cut
################################
$cmd="";
my $bam="$outdir/hisat2/$sample.sort.bam";
#java  -Xmx20G -jar /home/fuzl/soft/picard.jar AddOrReplaceReadGroups I=/data2/fuzl/project/508/demo/analysis/bwa/lib0705-7-2_S19_100kread_1.fq.sam  O=/data2/fuzl/project/508/demo/analysis/bwa/lib0705-7-2_S19_100kread.sort.bam  SO=coordinate  RGLB="pe"  RGPU="HiSeq-2000" RGPL=illumina RGSM=lib0705-7-2_S19_100kread
#$cmd .="java -Xmx20G -Djava.io.tmpdir=./ -jar $picard AddOrReplaceReadGroups I=$outdir/bwa/$sample.sam  O=$bam  SO=coordinate  RGLB=\"pe\"  RGPU=\"HiSeq-2000\" RGPL=illumina RGSM=$sample  1>$outdir/bwa/log.sort 2>&1 \n";
$cmd .="$samtools flagstat $bam > $outdir/hisat2/${sample}.alignment.flagstat & \n ";
$cmd .="$samtools stats  $bam > $outdir/hisat2/${sample}.alignment.stat  & \n";
$cmd .="$samtools index $bam \n";
#$cmd .="sed \'s\/\$\/\\t$sample\\t0\\t\\+\/\'  $bed > $outdir/bwa/bed.m \n";
$cmd .="awk '{print \$1\"\\t\"\$2\"\\t\"\$3\"\\ttmp\\t0\\t+\"}' $bed >$outdir/hisat2/bed.6 \n" ;
 
$cmd .="$qualimap  bamqc  -bam $bam -gff $outdir/hisat2/bed.6   -outdir ${bam}_QC -outfile $sample -outformat PDF:HTML --java-mem-size=10G \n";

&qsub("bamQC",$cmd,"20G");

#Summary
open S ,">$outdir/Summary.txt" or die $! ;
#my $Summary="Total_Sequences_fq1\tTotal_Sequences_fq2\tTotal_Sequences_fq1_trim\tTotal_Sequences_fq2_trim\tReads_mapped\tOn_target\tcoverage50X\tMean_coverageData\tGC%\n";
my $Summary="RawData\tClearData\tReads_mapped\tOn_target\tCoverage50X\tMean_coverageData\tGC%\tDuplicate\n";
my $Total_Sequences_fq1=`grep "Total Sequences"  $outdir/FASTQC/${sample}_1_fastqc/fastqc_data.txt |cut -f 2 `;chomp ($Total_Sequences_fq1);
my $Total_Sequences_fq2=`grep "Total Sequences"  $outdir/FASTQC/${sample}_2_fastqc/fastqc_data.txt |cut -f 2 `;chomp ($Total_Sequences_fq2);
my $Total_Sequences_fq1_trim=`grep "Total Sequences"  $outdir/FASTQC_trim/${sample}_1_paired_fastqc/fastqc_data.txt |cut -f 2 `;chomp ($Total_Sequences_fq1_trim);
my $Total_Sequences_fq2_trim=`grep "Total Sequences"  $outdir/FASTQC_trim/${sample}_2_paired_fastqc/fastqc_data.txt |cut -f 2 `;chomp ($Total_Sequences_fq2_trim);

my $Reads_mapped=`grep "number of mapped reads" ${bam}_QC/genome_results.txt |cut -d "=" -f 2 `;chomp ($Reads_mapped);
my ${coverage50X}=`grep "coverageData >= 50X" ${bam}_QC/genome_results.txt |awk '{print \$4}'`;chomp (${coverage50X});
my $mean_coverageData=`grep "mean coverageData" ${bam}_QC/genome_results.txt |awk '{print \$4}'|sed 's/X//'`;chomp ($mean_coverageData);
my $On_target= `sed -n 173p  ${bam}_QC/qualimapReport.html|cut -d ">" -f 2|cut -d "<" -f 1|sed 's/\\s\\+\\/\\s\\+/ \\(/'|sed 's/\$/\\)/'`;chomp ($On_target);
my $GC=`grep "GC percentage" ${bam}_QC/genome_results.txt |cut -d "=" -f 2 `;chomp ($GC);
my $Duplicate=`grep "number of duplicated reads" ${bam}_QC/genome_results.txt |cut -d "=" -f 2 `;chomp ($Duplicate);
my $RawData=$Total_Sequences_fq1+$Total_Sequences_fq2;
my $ClearData=$Total_Sequences_fq1_trim+$Total_Sequences_fq2_trim;
my $Reads_mapped_reads=(split " ", $Reads_mapped)[0];
$Duplicate=~s/,//g;
$Reads_mapped_reads=~s/,//g;
my $Duplicat_freq=100*$Duplicate/$Reads_mapped_reads;
$Duplicat_freq=sprintf("%4.2f" , $Duplicat_freq);
$Summary .="$RawData\t$ClearData\t$Reads_mapped\t$On_target\t$coverage50X\t$mean_coverageData\t$GC\t$Duplicate\($Duplicat_freq%\)\n";
print S "$Summary";
print "\nSummary done.\n";
$cmd="";
&MKDIR("$outdir/star_fusion");
$cmd .="$STAR_Fusion --genome_lib_dir $Fusion_database --left_fq $fq1 --right_fq $fq2 --output_dir $outdir/star_fusion --CPU 12 \n";
#/home/fuzl/soft/Fusion/STAR-Fusion-v1.5.0/STAR-Fusion --genome_lib_dir /home/fuzl/soft/Fusion/GRCh37_v19_CTAT_lib_Feb092018/ctat_genome_lib_build_dir --left_fq /home/fuzl/project/demo/lib-HD-4_1.fq.gz --right_fq /home/fuzl/project/demo/lib-HD-4_2.fq.gz --output_dir /home/fuzl/project/demo/new_HD-4_star_fusion_outdir --CPU 36
&qsub("STAR Fusion",$cmd,"20G");

`rm $outdir/bwa/$sample.sam` if (-f "$outdir/bwa/$sample.sam");


###############################################################################
my $duration_time=time-$start_time;
print strftime("End time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "This compute totally consumed $duration_time s\.\n";

###############################################################################
sub usage {
    die(
        qq!
Usage:
	eg:
perl $0  -fq1  /data/fuzl/project/demo/demo_1.fq -fq2 /data/fuzl/project/demo/demo_2.fq  -outdir /data/fuzl/project/demo/analysis &
Function: Template for Perl FASTQC BWA GATK pipeline .
Command:	-fq1 str	fq1
			-fq2 str	fq2
			-outdir	outdir
			-hg38	refrence genome version. defalt [hg19]
			-onlyprintcmd
			-bed  bed file 
 			-trim_q trimmomatic Quality   [13]

Author:   Zhiliang Fu fuzl\@geneis.cn, QQ:594380908
Version:  v1.0
Update:   2018/8/8
Notes: 
自动统计summary
判断shell finish，避免重复执行  
切5‘端前5bp 
hisat2 比对
start fusion 做融合
qualimap 对bam QC
应用于晓爽的四个外泌体数据：/data2/fuzl/lixs/exosome_20190325/out_gff3_all
\n!
    )
}
=c
sub qsub {
	my $name=shift @_;
	my $cmd=shift @_;
	my $vf=shift @_;
	&MKDIR("$outdir/shell/");
	my $n=$name;
	$n=~s/\s+/_/g ;
	open CMD ,">$outdir/shell/$n.sh" or die $!;
	print CMD "$cmd";
	close CMD;
	if (-e "/dev/sg17dddddddddddddddddddd"){ 
	# `qsub -cwd -l vf=4G -q all.q@sge_c17   $cmd `;
		system "cd $outdir/shell/ && perl /home/fuzl/script/qsub.pl -l vf=$vf -q all.q\@sge_c17  -N  $n  -s 60  -b 100  $outdir/shell/$n.sh " unless (defined $onlyprintcmd);
	}else{
        	system "sh $outdir/shell/$n.sh && touch $outdir/shell/$n.sh.finish " unless (defined $onlyprintcmd || -f "$outdir/shell/$n.sh.finish");
	#	system "sh $outdir/shell/$n.sh " unless (defined $onlyprintcmd) ;
	}
}
=cut

sub qsub { # 
	my $name=shift @_;
	my $cmd=shift @_;
	&MKDIR("$outdir/shell/");
	my $start_time=time;
#	print strftime("Start $name analysis time is %Y-%m-%d %H:%M:%S\n", localtime(time));
	&log_current_time("Start $name analysis...\n");
#	print "Start $name analysis ... \n";
	my $n=$name;
	$n=~s/\s+/_/g ;
	my $log_file="$outdir/shell/$n.log";
	open CMD ,">$outdir/shell/$n.sh" or die $!;
	print CMD "$cmd";
	close CMD;
	unless (defined $onlyprintcmd || -f "$outdir/shell/$n.sh.finish"){
		my $flag = system("sh $outdir/shell/$n.sh > $log_file") ;
	    if ($flag != 0 ){
	        &log_current_time("Error: command failed: $cmd");
	        exit(1);
	    } else {
	        my $escaped_time = (time()-$start_time)."s";
	        &log_current_time("$name done, escaped time: $escaped_time.\n");
	       	system "touch $outdir/shell/$n.sh.finish ";
	    }
	}else{
		&log_current_time("$name result has been finish, no duplicate run.");
	}
}
sub MKDIR{
	my $dir=shift @_;
	system "mkdir  -p $dir " unless (-d "$dir");
}
sub log_current_time {
     # get parameter
     my ($info) = @_;
     my $curr_time = date_time_format(localtime(time()));
     print "[$curr_time] $info\n";
}
sub date_time_format {
    my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
    return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}


#perl /home/fuzl/pipeline/capture/Fastqc_bwa_gatk3.8_508_qsub_v2.pl -fq1 /data2/fuzl/project/QC_debug/capture_diff/18CF90161/18CF90161P_S5_1.fq.gz -fq2 /data2/fuzl/project/QC_debug/capture_diff/18CF90161/18CF90161P_S5_2.fq.gz -outdir /data2/fuzl/project/QC_debug/capture_diff/18CF90161/analysis -gatktype 4 -onlyprintcmd
#perl //home/fuzl/project/508/script/Fastqc_bwa_gatk3.8_508_qsub_v2.pl -fq1 lib0705-7-2_S19_100kread_1.fq -fq2 lib0705-7-2_S19_100kread_2.fq -outdir /data2/fuzl/project/508/demo/analysis2
