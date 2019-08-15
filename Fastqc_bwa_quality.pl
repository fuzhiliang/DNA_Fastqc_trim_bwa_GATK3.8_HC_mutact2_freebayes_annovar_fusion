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

my ($genome,$INDEX,$dbSNP,$phasel_1KG,$Mills_and_1KG,$phasel_snp_1KG);
if (defined $hg38) { #dbSNP 库还有问题
	my $hg38_root="/home/fuzl/soft/GATK/resources/bundle/hg38/";
	$genome="$hg38_root/Homo_sapiens_assembly38.fasta";
	$INDEX="$hg38_root/bwa_index/gatk_hg38";
	$dbSNP="$hg38_root/dbsnp_146.hg38.vcf.gz";
	$phasel_1KG="$hg38_root/1000G_phase1.snps.high_confidence.hg38.vcf";
	$Mills_and_1KG="$hg38_root/Mills_and_1000G_gold_standard.indels.hg38.vcf";
}else {
	my $hg19_root="/home/fuzl/soft/GATK/resources/bundle/hg19";
	$genome="$hg19_root/ucsc.hg19.fasta";
	$INDEX="$hg19_root/bwa_index/gatk_hg19";
	$dbSNP="$hg19_root/dbsnp_138.hg19.vcf";	
	$phasel_snp_1KG="$hg19_root/1000G_phase1.snps.high_confidence.hg19.sites.vcf";
	$phasel_1KG="$hg19_root/1000G_phase1.indels.hg19.sites.vcf";
	$Mills_and_1KG="$hg19_root/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf";
}

my $picard="/home/fuzl/soft/picard.jar";
my $bin_trim_galore ="/home/fuzl/soft/TrimGalore-0.5.0/trim_galore";
#my $samtools="/home/wangb/samtools";
my $samtools="/home/fuzl/soft/samtools-1.9/samtools";
#my $GATK="/home/fuzl/soft/gatk-4.0.8.1/gatk";
my $GATK="/home/fuzl/soft/GATK/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar";
#my $bwa="/Share/home/tiangeng/.Work1/PERL/WGS/bwa-0.7.12/bwa";
my $bwa="/home/fuzl/soft/bwa-0.7.17/bwa";

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
$cmd .="java -jar /home/fuzl/soft/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 36 -phred33 -trimlog $outdir/trim/logfile ";
$cmd .="$fq1 $fq2 $outdir/trim/${fq_name}_1_paired.fq.gz  $outdir/trim/${fq_name}_1_unpaired.fq.gz  $outdir/trim/${fq_name}_2_paired.fq.gz $outdir/trim/${fq_name}_2_unpaired.fq.gz ";
#$cmd .="ILLUMINACLIP:/home/fuzl/soft/Trimmomatic-0.38/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 \n";
#$cmd .=" HEADCROP:3 ";
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
&MKDIR("$outdir/bwa");

$cmd="";
#$cmd .="mkdir $outdir/bwa \n" unless (-d "$outdir/bwa");
$cmd .="$bwa mem -M -t 36 -R '\@RG\\tID:$sample\\tSM:$sample\\tLB:$sample\\tPL:Illumina' $INDEX  $fq1 $fq2 >$outdir/bwa/$sample.sam  \n";
&qsub("BWA",$cmd,"50G");

################################
&MKDIR("$outdir/bwa/");
$cmd="";
my $bam="$outdir/bwa/$sample.sort.bam";
#java  -Xmx20G -jar /home/fuzl/soft/picard.jar AddOrReplaceReadGroups I=/data2/fuzl/project/508/demo/analysis/bwa/lib0705-7-2_S19_100kread_1.fq.sam  O=/data2/fuzl/project/508/demo/analysis/bwa/lib0705-7-2_S19_100kread.sort.bam  SO=coordinate  RGLB="pe"  RGPU="HiSeq-2000" RGPL=illumina RGSM=lib0705-7-2_S19_100kread
$cmd .="java -Xmx20G -Djava.io.tmpdir=./ -jar $picard AddOrReplaceReadGroups I=$outdir/bwa/$sample.sam  O=$bam  SO=coordinate  RGLB=\"pe\"  RGPU=\"HiSeq-2000\" RGPL=illumina RGSM=$sample  1>$outdir/bwa/log.sort 2>&1 \n";
$cmd .="$samtools flagstat $bam > $outdir/bwa/${sample}.alignment.flagstat & \n ";
$cmd .="$samtools stats  $bam > $outdir/bwa/${sample}.alignment.stat  & \n";
$cmd .="$samtools index $bam \n";
my $vf=50;
$cmd .="java -Xmx${vf}G -Djava.io.tmpdir=./ -jar $picard   MarkDuplicates  I=$bam  O=$outdir/bwa/${sample}_marked.bam M=$outdir/bwa/${sample}.metrics ";
$cmd .=" REMOVE_DUPLICATES=true "if ($REMOVE_DUPLICATES);
$cmd .=" 1>$outdir/bwa/log.mark  2>&1\n";
$bam="$outdir/bwa/${sample}_marked.bam";
$cmd .="$samtools index  $bam \n";

#$cmd .="sed \'s\/\$\/\\t$sample\\t0\\t\\+\/\'  $bed > $outdir/bwa/bed.m \n";
$cmd .="awk '{print \$1\"\\t\"\$2\"\\t\"\$3\"\\ttmp\\t0\\t+\"}' $bed >$outdir/bwa/bed.6 \n" ;
 
$cmd .="$qualimap  bamqc  -bam $bam -gff $outdir/bwa/bed.6   -outdir ${bam}_QC -outfile $sample -outformat PDF:HTML --java-mem-size=10G \n";

&qsub("bamQC",$cmd,"50G");

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
`rm $outdir/bwa/$sample.sam` if (-f "$outdir/bwa/$sample.sam");
=c
#################################
&MKDIR("$outdir/GATK_pretreatment");
$cmd ="";
$cmd .="java -Xmx20G -Djava.io.tmpdir=./ -jar $picard   MarkDuplicates  I=$bam  O=$outdir/GATK_pretreatment/${sample}_marked.bam M=$outdir/GATK_pretreatment/${sample}.metrics ";
$cmd .=" REMOVE_DUPLICATES=true "if ($REMOVE_DUPLICATES);
$cmd .=" 1>$outdir/GATK_pretreatment/log.mark  2>&1 \n";
$bam="$outdir/GATK_pretreatment/${sample}_marked.bam";
$cmd .="$samtools index  $bam \n";

#################################
$cmd .="java  -Xmx30G -Djava.io.tmpdir=./ -jar $GATK -T RealignerTargetCreator -R  $genome -I $bam -o $outdir/GATK_pretreatment/${sample}.intervals ";
$cmd .=" -known  $phasel_1KG " if ($phasel_1KG);
$cmd .=" -known  $Mills_and_1KG " if ($Mills_and_1KG);
$cmd .=" 1>$outdir/GATK_pretreatment/log.fix 2>&1 \n";

$cmd .="java -Xmx30G -Djava.io.tmpdir=./ -jar $GATK  -T IndelRealigner -R $genome -I $bam -o $outdir/GATK_pretreatment/${sample}.Realgn.bam  -targetIntervals $outdir/GATK_pretreatment/${sample}.intervals ";
$cmd .=" 1>$outdir/GATK_pretreatment/log.Realgn 2>&1 \n";
$bam="$outdir/GATK_pretreatment/${sample}.Realgn.bam";
#$cmd .="$samtools index  $bam \n";

##################################recal
my $recal=1;
if ($recal) {
	$cmd .="java -Xmx30G -Djava.io.tmpdir=./ -jar $GATK -T  BaseRecalibrator -R $genome  -I $bam -o $outdir/GATK_pretreatment/${sample}_recal_data.table  ";
	$cmd .="--knownSites   $dbSNP " ;
	$cmd .="--knownSites   $phasel_1KG " if ($phasel_1KG);
	$cmd .="--knownSites   $Mills_and_1KG " ;
	$cmd .=" 1>$outdir/GATK_pretreatment/log.recal 2>&1 \n";

	$cmd .="java -Xmx30G -Djava.io.tmpdir=./ -jar $GATK -T  PrintReads  -R $genome  -I $bam "; # 
	$cmd .="-BQSR $outdir/GATK_pretreatment/${sample}_recal_data.table "; 
	$cmd .="-o  $outdir/GATK_pretreatment/${sample}_recal.bam 1>$outdir/GATK_pretreatment/log.recal 2>&1 \n";

	$bam="$outdir/GATK_pretreatment/${sample}_recal.bam";
	$cmd .="$samtools index  $bam \n";
}
#################################
#$cmd .="$samtools depth $bam >$bam.detph\n";
#$cmd .="perl /home/fuzl/script/depth.pl  $bam.detph $bed 100  ${bam}.depth_windows_100  501 &\n";
#$cmd .="sed \'s\/\$\/\\t508\\t0\\t\\+\/\'  $bed > $outdir/GATK_pretreatment/bed.m \n";
#$cmd .="$qualimap  bamqc  -bam $bam -gff $outdir/GATK_pretreatment/bed.m   -outdir ${bam}_QC -outfile $sample -outformat PDF:HTML \n";
#perl /home/fuzl/script/depth.pl  *.bam.depth 38-cf1-genepanel.bed 10 *_windows_10  500

&qsub("GATK pretreatment",$cmd,"50G");

&MKDIR("$outdir/vcf/");
$gatktype||=4;
my $vcf ;
$cmd="";
if($gatktype == "1"){
	$cmd .="java -Xmx30G -Djava.io.tmpdir=./  -jar $GATK    -T   HaplotypeCaller  -R $genome -I $bam ";
	#$cmd .=" --dbsnp $dbSNP "; 
	$cmd .=" -o  $outdir/vcf/${sample}_gatk_HC_raw.vcf 1>$outdir/vcf/log.HC  2>&1 \n";
	$vcf="$outdir/vcf/${sample}_gatk_HC_raw.vcf ";
}elsif($gatktype == "2"){
	$cmd .="java -Xmx30G -Djava.io.tmpdir=./  -jar $GATK  -T  MuTect2  -R $genome -I:tumor $bam ";
#	$cmd .="-I:normal $normal_bam ;
	$cmd .=" --dbsnp $dbSNP "; 
	$cmd .=" -o $outdir/vcf/${sample}_Mutect2_raw.vcf 1>$outdir/vcf/log.SM  2>&1 \n";
	$vcf="$outdir/vcf/${sample}_Mutect2_raw.vcf ";
}elsif($gatktype == "3"){
	$cmd .="java -Xmx30G -Djava.io.tmpdir=./  -jar $GATK  -T UnifiedGenotyper  -R $genome -I $bam ";
	#$cmd .=" --dbsnp $dbSNP "; 
	$cmd .=" -o $outdir/vcf/${sample}_gatk_UG.raw.vcf 1>$outdir/vcf/log.UG  2>&1 \n";
	$vcf="$outdir/vcf/${sample}_gatk_UG.raw.vcf";
}elsif($gatktype == "4"){
	$cmd .="freebayes -j -m 10 -q 30 -F 0.001 -C 1 -t $bed -f $genome $bam > $outdir/vcf/$sample.freebayesp0.001.vcf \n";
	$cmd .="freebayes -j -m 10 -q 30 -F 0.01 -C 1 -t $bed -f $genome $bam > $outdir/vcf/$sample.freebayesp0.01.vcf \n";
#	freebayesp 4 freebayes -j -m 10 -q 20 -F 0.001 -C 1 -t /web/geneis/panel/scripts/38-cf1-genepanel.bed -f /database/ucsc.hg19.fa ./sort_merge/18HE25724F.rmdup.bam > ./analysis/18HE25724F---0.001.rmdup.vcf
#	freebayesp 4 freebayes -j -m 10 -q 20 -F 0.01 -C 1 -t /web/geneis/panel/scripts/38-cf1-genepanel.bed -f /database/ucsc.hg19.fa ./sort_merge/18HE25724F.rmdup.bam > ./analysis/18HE25724F---0.01.rmdup.vcf
	$cmd .= "perl /home/fuzl/soft/annovar/table_annovar.pl  $outdir/vcf/$sample.freebayesp0.001.vcf   /data2/fuzl/soft/annovar/humandb/ -buildver hg19 -out  $outdir/vcf/${sample}_freebayes0.001  -remove -protocol refGene,cytoBand,genomicSuperDups,esp6500siv2_all,1000g2014oct_all,1000g2014oct_afr,1000g2014oct_eas,1000g2014oct_eur,snp138,ljb26_all -operation g,r,r,f,f,f,f,f,f,f -nastring .  --vcfinput --otherinfo --dot2underline  & \n";
	$cmd .= "perl /home/fuzl/soft/annovar/table_annovar.pl  $outdir/vcf/$sample.freebayesp0.01.vcf   /data2/fuzl/soft/annovar/humandb/ -buildver hg19 -out  $outdir/vcf/${sample}_freebayes0.01  -remove -protocol refGene,cytoBand,genomicSuperDups,esp6500siv2_all,1000g2014oct_all,1000g2014oct_afr,1000g2014oct_eas,1000g2014oct_eur,snp138,ljb26_all -operation g,r,r,f,f,f,f,f,f,f -nastring .  --vcfinput --otherinfo --dot2underline &\n";
#	$cmd .="table_annovar.pl $outdir/vcf/$sample.freebayesp0.001.vcf /home/fuzl/soft/annovar/humandb/ -out $outdir/vcf/$sample.freebayesp0.001 -remove -protocol refGene,cytoBand,dbnsfp30a,cosmic83,snp144,1000g2015aug_all,1000g2015aug_afr,1000g2015aug_amr,1000g2015aug_eas,1000g2015aug_eur,1000g2015aug_sas,clinvar_20170905,esp6500siv2_all -operation g,r,f,f,f,f,f,f,f,f,f,f,f --buildver hg19 --nastring . --vcfinput --otherinfo --dot2underline\n";
#	$cmd .="table_annovar.pl $outdir/vcf/$sample.freebayesp0.01.vcf /home/fuzl/soft/annovar/humandb/ -out $outdir/vcf/$sample.freebayesp0.01 -remove -protocol refGene,cytoBand,dbnsfp30a,cosmic83,snp144,1000g2015aug_all,1000g2015aug_afr,1000g2015aug_amr,1000g2015aug_eas,1000g2015aug_eur,1000g2015aug_sas,clinvar_20170905,esp6500siv2_all -operation g,r,f,f,f,f,f,f,f,f,f,f,f --buildver hg19 --nastring . --vcfinput --otherinfo --dot2underline\n";
#table_annovar.pl ./analysis/18HE25724F---0.001.rmdup.vcf /database/humandb_hg19/ -out ./analysis/18HE25724F---0.001.rmdup -remove -protocol refGene,cytoBand,dbnsfp30a,cosmic83,snp144,1000g2015aug_all,1000g2015aug_afr,1000g2015aug_amr,1000g2015aug_eas,1000g2015aug_eur,1000g2015aug_sas,clinvar_20170905,esp6500siv2_all -operation g,r,f,f,f,f,f,f,f,f,f,f,f --buildver hg19 --nastring . --vcfinput --otherinfo --dot2underline
#table_annovar.pl ./analysis/18HE25724F---0.01.rmdup.vcf /database/humandb_hg19/ -out ./analysis/18HE25724F---0.01.rmdup -remove -protocol refGene,cytoBand,dbnsfp30a,cosmic83,snp144,1000g2015aug_all,1000g2015aug_afr,1000g2015aug_amr,1000g2015aug_eas,1000g2015aug_eur,1000g2015aug_sas,clinvar_20170905,esp6500siv2_all -operation g,r,f,f,f,f,f,f,f,f,f,f,f --buildver hg19 --nastring . --vcfinput --otherinfo --dot2underline
}else{
	die "Please select call mutation tools. \n"; 
}
################################
&runcmd("Genotype calling",$cmd,"50G");

unless ($gatktype == "4"){
	$cmd ="";
	$cmd .="java -jar $GATK -T  SelectVariants -R $genome -V $vcf  -selectType SNP  -o $outdir/vcf/${sample}_gatk_raw.snp.vcf  \n";
	$cmd .="java -jar $GATK -T  VariantFiltration -R $genome -V $outdir/vcf/${sample}_gatk_raw.snp.vcf   -filter  \"QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0\" -filterName \"my_snp_filter\"  ";
	$cmd .=" -o $outdir/vcf/${sample}.filter.snp.vcf 1>$outdir/vcf/log.snp.filter  2>&1 \n";
	$cmd .="grep -w -v my_snp_filter $outdir/vcf/${sample}.filter.snp.vcf > $outdir/vcf/${sample}.filter.PASS.snp.vcf \n";
	$cmd .="java -jar $GATK -T  SelectVariants -R $genome -V $vcf  -selectType INDEL  -o $outdir/vcf/${sample}_gatk_raw.InDel.vcf \n";
	$cmd .="java -jar $GATK -T  VariantFiltration -R $genome  -V $outdir/vcf/${sample}_gatk_raw.InDel.vcf    -filter \"QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0\"  -filterName \"my_indel_filter\"  ";
	$cmd .="-o $outdir/vcf/${sample}.filter.InDel.vcf 1>$outdir/vcf/log.InDel.filter  2>&1 \n";
	$cmd .="grep -w -v my_indel_filter $outdir/vcf/${sample}.filter.InDel.vcf > $outdir/vcf/${sample}.filter.PASS.InDel.vcf \n";

	#grep -v "^#"  lib-FZ18-04229F_S8.filter.PASS.InDel.vcf |cut -f 1-5,10-|sed 's/:/\t/g' |sed 's/,/\t/'|awk '{print $8/($7+$8)"\t"$0}' > lib-FZ18-04229F_S8.filter.PASS.InDel.vcf_freq
	#awk '{if ($7+$8 ==0){}else{print $8/($7+$8)"\t"$0}}

	#$cmd .= "grep -v \"^#\"  $outdir/vcf/${sample}.filter.PASS.snp.vcf |cut -f 1-5,10-|sed \'s\/:\/\\t\/g\' |sed \'s\/,\/\\t/\'|awk \'\{print \$8\/\(\$7+\$8\)\"\\t\"\$0\}\' \> $outdir/vcf/${sample}.filter.PASS.snp.vcf_freq \n";
	$cmd .= "grep -v \"^#\"  $outdir/vcf/${sample}.filter.PASS.snp.vcf |cut -f 1-5,10-|sed \'s\/:\/\\t\/g\' |sed \'s\/,\/\\t/\'|awk '\{if \(\$7\+\$8==0\)\{\}else\{print \$8\/\(\$7\+\$8\)\"\\t\"\$0\}\}\' \> $outdir/vcf/${sample}.filter.PASS.snp.vcf_freq \n";
	$cmd .= "grep -v \"^#\"  $outdir/vcf/${sample}.filter.PASS.InDel.vcf |cut -f 1-5,10-|sed \'s\/:\/\\t\/g\' |sed \'s\/,\/\\t/\'|awk '\{if \(\$7\+\$8==0\)\{\}else\{print \$8\/\(\$7\+\$8\)\"\\t\"\$0\}\}\' \> $outdir/vcf/${sample}.filter.PASS.InDel.vcf_freq \n";
	&qsub("Filter",$cmd,"10G");

	#&MKDIR("$outdir/ANNOVAR/");
	$cmd ="";
	$cmd .= "perl  /home/fuzl/soft/annovar/convert2annovar.pl -format  vcf4 $outdir/vcf/${sample}.filter.PASS.snp.vcf >$outdir/vcf/${sample}.filter.PASS.snp.vcf.avinput \n ";
	$cmd .= "perl /home/fuzl/soft/annovar/table_annovar.pl $outdir/vcf/${sample}.filter.PASS.snp.vcf.avinput   /data2/fuzl/soft/annovar/humandb/ -buildver hg19 -out  $outdir/vcf/${sample}_filter_PASS_snp_gatk -remove -protocol refGene,cytoBand,genomicSuperDups,esp6500siv2_all,1000g2014oct_all,1000g2014oct_afr,1000g2014oct_eas,1000g2014oct_eur,snp138,ljb26_all -operation g,r,r,f,f,f,f,f,f,f -nastring . \n";

	$cmd .= "perl  /home/fuzl/soft/annovar/convert2annovar.pl -format  vcf4  $outdir/vcf/${sample}.filter.PASS.InDel.vcf  >  $outdir/vcf/${sample}.filter.PASS.InDel.vcf.avinput \n ";
	$cmd .= "perl /home/fuzl/soft/annovar/table_annovar.pl  $outdir/vcf/${sample}.filter.PASS.InDel.vcf.avinput  /data2/fuzl/soft/annovar/humandb/ -buildver hg19 -out  $outdir/vcf/${sample}_filter_PASS_InDel_gatk  -remove -protocol refGene,cytoBand,genomicSuperDups,esp6500siv2_all,1000g2014oct_all,1000g2014oct_afr,1000g2014oct_eas,1000g2014oct_eur,snp138,ljb26_all -operation g,r,r,f,f,f,f,f,f,f -nastring . \n";

	$cmd .="awk \'NR==FNR\{a[\$2\"\\t\"\$3]=\$1\"\\t\"\$8\"\\t\"\$9\}NR>FNR\{if \(FNR==1\)\{print \"Freq\\tRefDP\\tAltDP\\t\"\$0\}else\{if(\$1\"\\t\"\$2 in a ){print a[\$1\"\\t\"\$2]\"\\t\"\$0}}}' " ;
	$cmd .=" $outdir/vcf/${sample}.filter.PASS.snp.vcf_freq  $outdir/vcf/${sample}_filter_PASS_snp_gatk.hg19_multianno.txt > $outdir/vcf/${sample}_filter_PASS_snp_gatk.hg19_multianno.txt.freq \n";

	#$cmd .="awk \'NR==FNR{if\(length\(\$6\)\<length\(\$5\)\){a[\$2\"\\t\"\$3+1]=\$1\"\\t\"\$8\"\\t\"\$9}else{a[\$2\"\\t\"\$3]=\$1\"\\t\"\$8\"\\t\"\$9}}NR>FNR\{if \(FNR==1\)\{print \"Freq\\tRefDP\\tAltDP\\t\"\$0\}else\{if(\$1\"\\t\"\$2 in a ){print a[\$1\"\\t\"\$2]\"\\t\"\$0}}}' " ;
	#$cmd .=" $outdir/vcf/${sample}.filter.PASS.InDel.vcf_freq   $outdir/vcf/${sample}_filter_PASS_InDel_gatk.hg19_multianno.txt > $outdir/vcf/${sample}_filter_PASS_InDel_gatk.hg19_multianno.txt.freq \n";
	$cmd .="awk \'NR==FNR{if\(length\(\$6\)\<length\(\$5\)\){a[\$2\"\\t\"\$3+1]=\$1\"\\t\"\$8\"\\t\"\$9}else{a[\$2\"\\t\"\$3]=\$1\"\\t\"\$8\"\\t\"\$9}}NR>FNR\{if \(FNR==1\)\{print \"Freq\\tRefDP\\tAltDP\\t\"\$0\}else\{if(\$1\"\\t\"\$2 in a ){print a[\$1\"\\t\"\$2]\"\\t\"\$0}}}' " ;
	$cmd .=" $outdir/vcf/${sample}.filter.PASS.InDel.vcf_freq   $outdir/vcf/${sample}_filter_PASS_InDel_gatk.hg19_multianno.txt > $outdir/vcf/${sample}_filter_PASS_InDel_gatk.hg19_multianno.txt.freq \n";

	&qsub("ANNOVAR",$cmd,"10G");
}
=cut

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
			#-gatktype gatk type,  1=HaplotypeCaller, 2=Mutect2, 3=UnifiedGenotype, 4=freebayes  [4p]
			-onlyprintcmd
			-bed  bed file 
 			-trim_q trimmomatic Quality   [13]
 			-REMOVE_DUPLICATES  picard REMOVE DUPLICATES  []

Author:   Zhiliang Fu fuzl\@geneis.cn, QQ:594380908
Version:  v1.0
Update:   2018/8/8
Notes: 自动统计summary，
判断shell finish，避免重复执行 
去重  
\n!
    )
}

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
	if (-e "/dev/sg17"){ 
	# `qsub -cwd -l vf=4G -q all.q@sge_c17   $cmd `;
		system "cd $outdir/shell/ && perl /home/fuzl/script/qsub.pl -l vf=$vf -q all.q\@sge_c17  -N  $n  -s 60  -b 100  $outdir/shell/$n.sh " unless (defined $onlyprintcmd);
	}else{
        	system "sh $outdir/shell/$n.sh && touch $outdir/shell/$n.sh.finish " unless (defined $onlyprintcmd || -f "$outdir/shell/$n.sh.finish");
	#	system "sh $outdir/shell/$n.sh " unless (defined $onlyprintcmd) ;
	}
}


sub runcmd { # 
	my $name=shift @_;
	my $cmd=shift @_;
	&MKDIR("$outdir/shell/");
	print "Start $name analysis ... \n";
	my $n=$name;
	$n=~s/\s+/_/g ;
	open CMD ,">$outdir/shell/$n.sh" or die $!;
	print CMD "$cmd";
	close CMD;
        system "sh $outdir/shell/$n.sh && touch $outdir/shell/$n.sh.finish " unless (defined $onlyprintcmd || -f "$outdir/shell/$n.sh.finish");
	#system "sh $outdir/shell/$n.sh " unless (defined $onlyprintcmd);
	print "End $name analysis !\n";
}

sub MKDIR{
	my $dir=shift @_;
	system "mkdir  -p $dir " unless (-d "$dir");
}
#perl /home/fuzl/pipeline/capture/Fastqc_bwa_gatk3.8_508_qsub_v2.pl -fq1 /data2/fuzl/project/QC_debug/capture_diff/18CF90161/18CF90161P_S5_1.fq.gz -fq2 /data2/fuzl/project/QC_debug/capture_diff/18CF90161/18CF90161P_S5_2.fq.gz -outdir /data2/fuzl/project/QC_debug/capture_diff/18CF90161/analysis -gatktype 4 -onlyprintcmd
#perl //home/fuzl/project/508/script/Fastqc_bwa_gatk3.8_508_qsub_v2.pl -fq1 lib0705-7-2_S19_100kread_1.fq -fq2 lib0705-7-2_S19_100kread_2.fq -outdir /data2/fuzl/project/508/demo/analysis2
