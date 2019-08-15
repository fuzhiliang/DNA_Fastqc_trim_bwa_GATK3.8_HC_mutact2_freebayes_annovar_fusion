#!/usr/bin/perl -w
use strict;
use warnings;
use POSIX qw(strftime);
use Getopt::Std;
use Getopt::Long;
use File::Basename;
use Cwd 'abs_path';
my $BEGIN_TIME=time();
my $verbose ="v1.0";


###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my ($fq1,$fq2,$outdir,$hg38,$onlyprintcmd,$gatktype,$bed,$SE_bed,$trim_q,$fusion,$REMOVE_DUPLICATES,$Amplica,$huada_adapt,$threads,$vf);
my %opts;
GetOptions(
    'verbose' => \$verbose,
    'fq1=s' => \$fq1,
	'fq2=s' => \$fq2,
    'outdir=s' => \$outdir, 
	'hg38' => \$hg38,
	'gatktype=i' => \$gatktype,
    'onlyprintcmd' => \$onlyprintcmd,
 	'bed=s' => \$bed,
 	'SE_bed=s' => \$SE_bed,  
 	'trim_q=i'=>\$trim_q, 
 	'REMOVE_DUPLICATES'=>\$REMOVE_DUPLICATES,
 	'Amplica'=>\$Amplica,
 	'huada_adapt'=>\$huada_adapt,
 	'vf=i'=>\$vf,
 	'fusion'=>\$fusion,
 	'threads=i'=>\$threads,
) or die $!;
unless(defined $fq2 && defined $fq1 ){&usage();exit 0;}
##########################
my $fastqc="/home/fuzl/soft/FastQC/fastqc";
$bed||="/home/fuzl/bed/Illumina_WES.bed"; #Exon
$trim_q||="13";
$vf||="10";
$threads||="36";
my $qualimap="/home/fuzl/soft/qualimap_v2.2.1/qualimap" ;
$SE_bed||="/home/fuzl/script/SE-GeneFusion-Evan/database/fusion_38cf.bed";

my ($genome,$INDEX,$dbSNP,$phasel_1KG,$Mills_and_1KG,$phasel_snp_1KG,$Fusion_database,$cosmic);
if (defined $hg38) { #dbSNP 库还有问题
	my $hg38_root="/home/fuzl/soft/GATK/resources/bundle/hg38/";
	$genome="$hg38_root/Homo_sapiens_assembly38.fasta";
	$INDEX="$hg38_root/bwa_index/gatk_hg38";
	$dbSNP="$hg38_root/dbsnp_146.hg38.vcf.gz";
	$phasel_1KG="$hg38_root/1000G_phase1.snps.high_confidence.hg38.vcf";
	$Mills_and_1KG="$hg38_root/Mills_and_1000G_gold_standard.indels.hg38.vcf";
	$Fusion_database="/home/fuzl/soft/Fusion/GRCh38_v27_CTAT_lib_Feb092018/ctat_genome_lib_build_dir/";
}else {
	my $hg19_root="/home/fuzl/soft/GATK/resources/bundle/hg19";
	$genome="$hg19_root/ucsc.hg19.fasta";
	$INDEX="$hg19_root/bwa_index/gatk_hg19";
	$dbSNP="$hg19_root/dbsnp_138.hg19.vcf";	
	$cosmic="$hg19_root/CosmicCodingMuts.chr.sort.head.vcf";
	$phasel_snp_1KG="$hg19_root/1000G_phase1.snps.high_confidence.hg19.sites.vcf";
	$phasel_1KG="$hg19_root/1000G_phase1.indels.hg19.sites.vcf";
	$Mills_and_1KG="$hg19_root/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf";
	$Fusion_database="/home/fuzl/soft/Fusion/GRCh37_v19_CTAT_lib_Feb092018/ctat_genome_lib_build_dir";
}

my $picard="/home/fuzl/soft/picard.jar";
my $bin_trim_galore ="/home/fuzl/soft/TrimGalore-0.5.0/trim_galore";
my $samtools="/home/fuzl/soft/samtools-1.9/samtools";
#my $GATK="/home/fuzl/soft/gatk-4.0.8.1/gatk";
my $GATK="/home/fuzl/soft/GATK/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar";
#my $bwa="/Share/home/tiangeng/.Work1/PERL/WGS/bwa-0.7.12/bwa";
my $bwa="/home/fuzl/soft/bwa-0.7.17/bwa";
my $STAR_Fusion="/home/fuzl/soft/Fusion/STAR-Fusion-v1.5.0/STAR-Fusion"; 

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
$fq_name=~s/_(R)?1\.f(ast)?q(\.gz)?$//;
&MKDIR("$outdir/FASTQC");
$cmd .="$fastqc $fq1 -o $outdir/FASTQC/  && cd $outdir/FASTQC/ &&  unzip -o $outdir/FASTQC/${fq_name}_1_fastqc.zip \n";
$cmd .="$fastqc $fq2 -o $outdir/FASTQC/  && cd $outdir/FASTQC/ &&  unzip -o $outdir/FASTQC/${fq_name}_2_fastqc.zip \n ";

&MKDIR("$outdir/trim/");
$cmd .="java -jar /home/fuzl/soft/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads $threads -phred33 -trimlog $outdir/trim/logfile ";
$cmd .="$fq1 $fq2 $outdir/trim/${fq_name}_1_paired.fq.gz  $outdir/trim/${fq_name}_1_unpaired.fq.gz  $outdir/trim/${fq_name}_2_paired.fq.gz $outdir/trim/${fq_name}_2_unpaired.fq.gz ";
$cmd .=" HEADCROP:3 " if ($Amplica);
$cmd .=" CROP:147 " if ($Amplica) ; # reads length 151 bp，切除3‘端3bp
my $adapt_seq="/home/fuzl/soft/Trimmomatic-0.38/adapters/TruSeq3-PE.fa";
$adapt_seq="/home/fuzl/soft/Trimmomatic-0.38/adapters/huada-PE.fa" if ($huada_adapt); 
$cmd .="ILLUMINACLIP:${adapt_seq}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:$trim_q MINLEN:36 \n";# 修改LEADING 和 SLIDINGWINDOW
#$cmd .="ILLUMINACLIP:/home/fuzl/soft/Trimmomatic-0.38/adapters/huada-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:$trim_q MINLEN:36 \n";# 修改LEADING 和 SLIDINGWINDOW
$fq1="$outdir/trim/${fq_name}_1_paired.fq.gz";
$fq2="$outdir/trim/${fq_name}_2_paired.fq.gz";
$fq_name=basename($fq1);
$fq_name=~s/_1_paired.f(ast)?q(\.gz)?$//;
&MKDIR("$outdir/FASTQC_trim/");
$cmd .="$fastqc $fq1 -o $outdir/FASTQC_trim/  && cd $outdir/FASTQC_trim/ && unzip -o $outdir/FASTQC_trim/${fq_name}_1_paired_fastqc.zip \n";
$cmd .="$fastqc $fq2 -o $outdir/FASTQC_trim/  && cd $outdir/FASTQC_trim/ && unzip -o $outdir/FASTQC_trim/${fq_name}_2_paired_fastqc.zip \n ";

&runcmd("QC",$cmd);

#############################################################
#bwa
my $sample=$fq_name;
&MKDIR("$outdir/bwa");

$cmd="";
$cmd .="$bwa mem -M -t $threads -R '\@RG\\tID:$sample\\tSM:$sample\\tLB:$sample\\tPL:Illumina' $INDEX  $fq1 $fq2 >$outdir/bwa/$sample.sam  \n";
&runcmd("BWA",$cmd);

################################
&MKDIR("$outdir/GATK_pretreatment/");
$cmd="";
my $bam="$outdir/GATK_pretreatment/$sample.sort.bam";
#java  -Xmx20G -jar /home/fuzl/soft/picard.jar AddOrReplaceReadGroups I=/home/fuzl/project/508/demo/analysis/bwa/lib0705-7-2_S19_100kread_1.fq.sam  O=/home/fuzl/project/508/demo/analysis/bwa/lib0705-7-2_S19_100kread.sort.bam  SO=coordinate  RGLB="pe"  RGPU="HiSeq-2000" RGPL=illumina RGSM=lib0705-7-2_S19_100kread
$cmd .="java -Xmx${vf}G -Djava.io.tmpdir=./ -jar $picard AddOrReplaceReadGroups I=$outdir/bwa/$sample.sam  O=$bam  SO=coordinate  RGLB=\"pe\"  RGPU=\"HiSeq-2000\" RGPL=illumina RGSM=$sample  1>$outdir/GATK_pretreatment/log.sort 2>&1 \n";

$cmd .="$samtools index $bam \n";

&runcmd("GATK_1 Picard sort",$cmd);

#################################
unless ($Amplica){
	$cmd ="";
	$cmd .="java -Xmx${vf}G -Djava.io.tmpdir=./ -jar $picard   MarkDuplicates  I=$bam  O=$outdir/GATK_pretreatment/${sample}_marked.bam M=$outdir/GATK_pretreatment/${sample}.metrics ";
	$cmd .=" REMOVE_DUPLICATES=true "if ($REMOVE_DUPLICATES);
	$cmd .=" 1>$outdir/GATK_pretreatment/log.mark  2>&1\n";
	$bam="$outdir/GATK_pretreatment/${sample}_marked.bam";
	$cmd .="$samtools index  $bam \n";

	&runcmd("GATK_2 MarkDuplicates",$cmd);

#################################
	$cmd="";
	$cmd .="java  -Xmx${vf}G -Djava.io.tmpdir=./ -jar $GATK -T RealignerTargetCreator -R  $genome -I $bam -o $outdir/GATK_pretreatment/${sample}.intervals ";
	$cmd .=" -known  $phasel_1KG " if ($phasel_1KG);
	$cmd .=" -known  $Mills_and_1KG " if ($Mills_and_1KG);
	$cmd .=" 1>$outdir/GATK_pretreatment/log.fix 2>&1  \n";

	$cmd .="java -Xmx${vf}G -Djava.io.tmpdir=./ -jar $GATK  -T IndelRealigner -R $genome -I $bam -o $outdir/GATK_pretreatment/${sample}.Realgn.bam  -targetIntervals $outdir/GATK_pretreatment/${sample}.intervals ";
	$cmd .=" 1>$outdir/GATK_pretreatment/log.Realgn  2>&1 \n";
	$bam="$outdir/GATK_pretreatment/${sample}.Realgn.bam";
#$cmd .="$samtools index  $bam \n";

	&runcmd("GATK_3 IndelRealigner",$cmd);

##################################recal
	$cmd="";
	$cmd .="java -Xmx${vf}G -Djava.io.tmpdir=./ -jar $GATK -T  BaseRecalibrator -R $genome  -I $bam -o $outdir/GATK_pretreatment/${sample}_recal_data.table  ";
	$cmd .="--knownSites   $dbSNP " ;
	$cmd .="--knownSites   $phasel_1KG " if ($phasel_1KG);
	$cmd .="--knownSites   $Mills_and_1KG " ;
	$cmd .=" 1>$outdir/GATK_pretreatment/log.recal_table 2>&1 \n";

	$cmd .="java -Xmx${vf}G -Djava.io.tmpdir=./ -jar $GATK -T  PrintReads  -R $genome  -I $bam "; # 
	$cmd .="-BQSR $outdir/GATK_pretreatment/${sample}_recal_data.table "; 
	$cmd .="-o  $outdir/GATK_pretreatment/${sample}_recal.bam 1>$outdir/GATK_pretreatment/log.recal  2>&1 \n";

	$bam="$outdir/GATK_pretreatment/${sample}_recal.bam";
	$cmd .="$samtools index  $bam \n";

	&runcmd("GATK_4 Recal",$cmd);
}
#################################
$cmd="";
#$cmd .="$samtools depth $bam >$bam.detph\n";
#$cmd .="perl /home/fuzl/script/depth.pl  $bam.detph $bed 100  ${bam}.depth_windows_100  500 &\n";
#$cmd .="sed \'s\/\$\/\\ttmp\\t0\\t\\+\/\'  $bed > $outdir/GATK_pretreatment/bed.m \n";
$cmd .="$samtools flagstat $bam > $outdir/GATK_pretreatment/${sample}.alignment.flagstat & \n ";
$cmd .="$samtools stats -d $bam > $outdir/GATK_pretreatment/${sample}.alignment.stat   \n";
$cmd .="awk \'\{print \$1\"\\t\"\$2\"\\t\"\$3\"\\ttmp\\t0\\t\+\"\}\'  $bed > $outdir/GATK_pretreatment/bed.m \n";
$cmd .="$qualimap  bamqc  -bam $bam -gff $outdir/GATK_pretreatment/bed.m   -outdir ${bam}_QC -outfile $sample -outformat PDF:HTML --java-mem-size=${vf}G -sd \n";
&runcmd("Qualimap",$cmd);

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
my $On_target= `sed -n 176p  ${bam}_QC/qualimapReport.html|cut -d ">" -f 2|cut -d "<" -f 1|sed 's/\\s\\+\\/\\s\\+/ \\(/'|sed 's/\$/\\)/'`;chomp ($On_target);
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
######################################################

&MKDIR("$outdir/vcf/");
$gatktype||=4;
my $vcf ;
my $freebayes="/home/fuzl/script/freebayesp.py  $threads /home/fuzl/miniconda2/bin/freebayes ";
my $humandb="/home/fuzl/soft/annovar/humandb_hg19_20190121 "; # 新版本注释数据库,注意修改后注释命令行对应也要改

$cmd="";
if($gatktype == "1"){
	$cmd .="java -Xmx${vf}G -Djava.io.tmpdir=./  -jar $GATK    -T   HaplotypeCaller  -R $genome -I $bam ";
	#$cmd .=" --dbsnp $dbSNP "; 
	$cmd .=" -o  $outdir/vcf/${sample}_gatk_HC_raw.vcf 1>$outdir/vcf/log.HC  2>&1 \n";
	$vcf="$outdir/vcf/${sample}_gatk_HC_raw.vcf ";
}elsif($gatktype == "2"){
	$cmd .="java -Xmx${vf}G -Djava.io.tmpdir=./  -jar $GATK  -T  MuTect2  -R $genome -I:tumor $bam ";
#	$cmd .="-I:normal $normal_bam ;
	$cmd .=" --dbsnp $dbSNP "; 
	$cmd .=" --cosmic $cosmic ";
	$cmd .=" --intervals $bed ";
	#$cmd .=" -PON /path/to/EGA/normal/MuTect2_PON.vcf "
	$cmd .=" -o $outdir/vcf/${sample}_Mutect2_raw.vcf 1>$outdir/vcf/log.SM  2>&1 \n";
	$vcf="$outdir/vcf/${sample}_Mutect2_raw.vcf ";
}elsif($gatktype == "3"){
	$cmd .="java -Xmx${vf}G -Djava.io.tmpdir=./  -jar $GATK  -T UnifiedGenotyper  -R $genome -I $bam ";
	#$cmd .=" --dbsnp $dbSNP "; 
	$cmd .=" -o $outdir/vcf/${sample}_gatk_UG.raw.vcf 1>$outdir/vcf/log.UG  2>&1 \n";
	$vcf="$outdir/vcf/${sample}_gatk_UG.raw.vcf";
}elsif($gatktype == "4"){
	$cmd .="$freebayes -j -m 10 -q 30 -F 0.001 -C 1 -t $bed -i --no-indels -f $genome $bam > $outdir/vcf/$sample.snp.freebayesp_q30.vcf \n";
	$cmd .="$freebayes -j -m 10 -q 20 -F 0.001 -C 1 -t $bed -I  -f $genome $bam > $outdir/vcf/$sample.freebayesp_q20.vcf \n";

	$cmd .="cat $outdir/vcf/$sample.snp.freebayesp_q30.vcf $outdir/vcf/$sample.freebayesp_q20.vcf  >$outdir/vcf/$sample.vcf \n";
	$cmd .= "perl /home/fuzl/soft/annovar/table_annovar.pl $outdir/vcf/$sample.vcf  $humandb -buildver hg19 "; 
	$cmd .= "-out  $outdir/vcf/${sample}  -remove -protocol refGene,cytoBand,dbnsfp30a,cosmic83,snp144,1000g2015aug_all,1000g2015aug_afr,1000g2015aug_amr,1000g2015aug_eas,1000g2015aug_eur,1000g2015aug_sas,clinvar_20170905,esp6500siv2_all ";
	$cmd .= "-operation g,r,f,f,f,f,f,f,f,f,f,f,f --buildver hg19 --nastring . --vcfinput --otherinfo --dot2underline \n";
	$cmd .= "perl /home/fuzl/script/split_multi-mut-with-header-new-20181114.pl $outdir/vcf/${sample}.hg19_multianno.txt > $outdir/vcf/${sample}.hg19_multianno.txt.freq \n";

}else{
	die "Please select call mutation tools. \n"; 
}
&runcmd("GATK_5 Genotype calling",$cmd);
################################

unless ($gatktype == "4"){

	$cmd ="";
	$cmd .="java -jar $GATK -T  SelectVariants -R $genome -V $vcf  -selectType SNP  -o $outdir/vcf/${sample}_gatk_raw.snp.vcf  \n";
	$cmd .="java -jar $GATK -T  SelectVariants -R $genome -V $vcf  -selectType INDEL  -o $outdir/vcf/${sample}_gatk_raw.InDel.vcf \n";
	$cmd .= "perl /home/fuzl/soft/annovar/table_annovar.pl  $outdir/vcf/${sample}_gatk_raw.snp.vcf   $humandb  -buildver hg19 ";
	$cmd .=" -out  $outdir/vcf/${sample}_gatk_raw.snp  -remove ";
	$cmd .= " -protocol refGene,cytoBand,dbnsfp30a,cosmic83,snp144,1000g2015aug_all,1000g2015aug_afr,1000g2015aug_amr,1000g2015aug_eas,1000g2015aug_eur,1000g2015aug_sas,clinvar_20170905,esp6500siv2_all ";
	$cmd .= "perl /home/fuzl/soft/annovar/table_annovar.pl  $outdir/vcf/${sample}_gatk_raw.InDel.vcf   $humandb -buildver hg19 -out  $outdir/vcf/${sample}_gatk_raw.InDel ";
	$cmd .= " -remove -protocol refGene,cytoBand,dbnsfp30a,cosmic83,snp144,1000g2015aug_all,1000g2015aug_afr,1000g2015aug_amr,1000g2015aug_eas,1000g2015aug_eur,1000g2015aug_sas,clinvar_20170905,esp6500siv2_all ";
	$cmd .= "perl /home/fuzl/script/split_multi-mut-with-header-new-20181114.pl $outdir/vcf/${sample}_gatk_raw.snp.hg19_multianno.txt > $outdir/vcf/${sample}_gatk_raw.snp.hg19_multianno.txt.freq \n";
	$cmd .= "perl /home/fuzl/script/split_multi-mut-with-header-new-20181114.pl $outdir/vcf/${sample}_gatk_raw.InDel.hg19_multianno.txt > $outdir/vcf/${sample}_gatk_raw.indel.hg19_multianno.txt.freq \n";

	&runcmd("ANNOVAR",$cmd);
}

######################
$cmd="";
$cmd.="perl /home/fuzl/script/SE-GeneFusion-Evan/SEGF.pl -fq1 $fq1 -fq2 $fq2  -odir $outdir/SE_$sample -trim_len 10 -remain_len 35 ";
$cmd.="-bed $SE_bed -process $threads \n";
#perl /data2/fuzl/script/SE-GeneFusion-Evan/SEGF.pl -fq1 /data2/fuzl/project/SE/rawdata/lib-FZ19-00886F_S8_1.fq.gz -fq2 /data2/fuzl/project/SE/rawdata/lib-FZ19-00886F_S8_2.fq.gz -odir /data2/fuzl/project/SE/rawdata/lib-FZ19-00886F_S8 -trim_len 10 -remain_len 35 -bed /home/fuzl/script/SE-GeneFusion-Evan/database/fusion_38cf.bed -process 36 &
$cmd .="$STAR_Fusion --genome_lib_dir $Fusion_database --left_fq $fq1 --right_fq $fq2 --output_dir $outdir/star_fusion_$sample --CPU $threads \n" if (defined $fusion);
#/home/fuzl/soft/Fusion/STAR-Fusion-v1.5.0/STAR-Fusion --genome_lib_dir /home/fuzl/soft/Fusion/GRCh37_v19_CTAT_lib_Feb092018/ctat_genome_lib_build_dir --left_fq /home/fuzl/project/demo/lib-HD-4_1.fq.gz --right_fq /home/fuzl/project/demo/lib-HD-4_2.fq.gz --output_dir /home/fuzl/project/demo/new_HD-4_star_fusion_outdir --CPU 36
&runcmd("Fusion",$cmd);

&MKDIR("$outdir/CNV");
$cmd ="cnvkit.py batch $bam  --method hybrid --targets $bed ";
$cmd .=" --normal /home/fuzl/project/boke519_wangyy_v2/exon_intron/lib0dian9percenttest_L4/GATK_pretreatment/lib0dian9percenttest_L4_recal.bam  ";
$cmd .=" --annotate /home/zhangsw/reference/annotation/hg19/gencode.allgene.position.bed ";
$cmd .=" --fasta $genome ";
$cmd .=" --output-reference $outdir/CNV/lib0dian9percenttest_cnvkit.cnn --output-dir $outdir/CNV/ \n";
#$cmd .="--diagram --scatter \n"; # 192.168.10.2 上图形界面bug，否则需加上该参

&runcmd("CNV",$cmd) if (defined $fusion);

`rm $outdir/bwa/$sample.sam ` if (-f "$outdir/bwa/$sample.sam"); 
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
perl $0   -fq1 /home/fuzl/demo/demo1_1.fq.gz -fq2 /home/fuzl/demo/demo1_2.fq.gz -outdir /home/fuzl/demo/demo1 -gatktype 1 &
Function: Template for Perl FASTQC BWA GATK pipeline .
Command:	-fq1 str	fq1
			-fq2 str	fq2
			-outdir	outdir
			-hg38	refrence genome version. defalt [hg19]
			-gatktype gatk type,  1=HaplotypeCaller, 2=Mutect2, 3=UnifiedGenotype, 4=freebayes  [4]
			-onlyprintcmd
			-bed   bed file 
			-SE_bed  bed of SE fusion   [fusion_38cf.bed]
 			-trim_q trimmomatic Quality   [13]
 			-REMOVE_DUPLICATES  picard REMOVE DUPLICATES  []
 			-Amplica    Amplica seq    []
 			-huada_adapt  huada MGI adapers   [Illumina]
 			-threads		[36]
 			-vf  			[10]G		
Author:   Zhiliang Fu fuzl\@geneis.cn, QQ:594380908
Version:  v1.0
Update:   2018/8/8
Version:  v3.0
Update:   2019/2/27

Notes:    192.168.10.11
修改snp 和 indel 参数 #2019.01.10
snp q 30 -F 0.001
InDel q 20 -F 0.001
切首端3bp 加参数-Amplica
华大的接头 加参数-huada_adapt
添加vf 参数，默认10G
添加线程参数
freebayes 多线程
单步添加finish
gatk 拆分多个shell
annovar 注释数据库版本更新
添加报错程序中断机制
star_fusion 在11 上跑不了
添加SE融合检测，用的是trim后的reads
修复Haplotype 注释，数据库版本的bug，只能用新的annovar 数据库
#v3.2
freebayes 单独call snp 和indel，跟生产流程一致，但是存在bug，有的位点q30call不出来，q20能call出来的，结果中不会有该位点
#v3.3
针对519panel添加cnv分析
freebayes q20 0.01 call indel
freebyes 多线程改用/home/fuzl/script/freebayesp.py
qualimap 加-sd ,考虑dup
\n!
    )
}

sub runcmd { # 
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
	        log_current_time("Error: command failed: $cmd");
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

     # get current time with string
     my $curr_time = date_time_format(localtime(time()));

     # print info with time
     print "[$curr_time] $info\n";
}

#############################################################################################################
sub date_time_format {
    my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
    return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

#perl /home/fuzl/pipeline/capture/Fastqc_bwa_gatk3.8_508_qsub_v2.pl -fq1 /home/fuzl/project/QC_debug/capture_diff/18CF90161/18CF90161P_S5_1.fq.gz -fq2 /home/fuzl/project/QC_debug/capture_diff/18CF90161/18CF90161P_S5_2.fq.gz -outdir /home/fuzl/project/QC_debug/capture_diff/18CF90161/analysis -gatktype 4 -onlyprintcmd
#perl //home/fuzl/project/508/script/Fastqc_bwa_gatk3.8_508_qsub_v2.pl -fq1 lib0705-7-2_S19_100kread_1.fq -fq2 lib0705-7-2_S19_100kread_2.fq -outdir /home/fuzl/project/508/demo/analysis2
