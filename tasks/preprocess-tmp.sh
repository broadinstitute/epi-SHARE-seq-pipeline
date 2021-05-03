if [ "$Start" = Fastq ]; then
    echo "Skip bcltofastq"
    if [ -f $dir/fastqs/filesi1.xls ]; then rm $dir/fastqs/filler*; fi
    cd $rawdir
    if ls *L001_R1_001.fastq.gz 1> /dev/null 2>&1; then
	temp=$(ls *_S1_L001_R1_001.fastq.gz)
	Run=$(echo $temp | sed -e 's/\_S1\_L001\_R1\_001.fastq.gz//')
	singlelane=F
	temp=$(ls *_S1_L00*_R1_001.fastq.gz)
	VAR=( $temp )
	nolane=${#VAR[@]}
	echo "Detected $nolane lanes"
    elif ls *S1_R1_001.fastq.gz 1> /dev/null 2>&1; then
	echo "Detected single lane"
	temp=$(ls *S1_R1_001.fastq.gz)
	Run=$(echo $temp | sed -e 's/\_\S1\_\R1\_\001.fastq.gz//')
	singlelane=T
	nolane=1
    else
	echo "No fastq with matched naming format detected; exit..."
	exit
    fi
    
    echo "Run number is:" $Run

    # split fastqs
    mkdir $dir/smallfastqs/
    if [ -f $dir/smallfastqs/0001.1.$Run.R2.fastq.gz ]; then
        echo "Found 0001.$Run.R2.fastq, skip split fastqs"
    else
	if [ ! -f $dir/R1.1.fastq.gz ]; then
            echo "Link fastqs"
	    if [ $singlelane == T ]; then
		ln -s $rawdir/"$Run"_S1_R1_001.fastq.gz $dir/R1.1.fastq.gz
		ln -s $rawdir/"$Run"_S1_R4_001.fastq.gz $dir/R2.1.fastq.gz
		ln -s $rawdir/"$Run"_S1_R2_001.fastq.gz $dir/I1.1.fastq.gz
		ln -s $rawdir/"$Run"_S1_R3_001.fastq.gz $dir/I2.1.fastq.gz
	    else
		parallel 'ln -s '$rawdir'/'$Run'_S1_L00{}_R1_001.fastq.gz '$dir'/R1.{}.fastq.gz' ::: $(seq $nolane)
		parallel 'ln -s '$rawdir'/'$Run'_S1_L00{}_R4_001.fastq.gz '$dir'/R2.{}.fastq.gz' ::: $(seq $nolane)
		parallel 'ln -s '$rawdir'/'$Run'_S1_L00{}_R2_001.fastq.gz '$dir'/I1.{}.fastq.gz' ::: $(seq $nolane)
		parallel 'ln -s '$rawdir'/'$Run'_S1_L00{}_R3_001.fastq.gz '$dir'/I2.{}.fastq.gz' ::: $(seq $nolane)
	    fi
	fi
	
	if [ "$Runtype" = full ]; then
	    # Runing full pipeline
	    echo "Split fastqs to small files"
#	    $fastpPath/fastp -i $dir/R1.fastq.gz -o $dir/smallfastqs/$Run.1.R1.fastq.gz -S 12000000 --thread 1 -d 4 -A -G -L -Q 2>>$dir/split.log &
#	    $fastpPath/fastp -i $dir/R2.fastq.gz -o $dir/smallfastqs/$Run.1.R2.fastq.gz -S 12000000 --thread 1 -d 4 -A -G -L -Q 2>>$dir/split.log &
#	    $fastpPath/fastp -i $dir/I1.fastq.gz -o $dir/smallfastqs/$Run.1.I1.fastq.gz -S 12000000 --thread 1 -d 4 -A -G -L -Q 2>>$dir/split.log &
#	    $fastpPath/fastp -i $dir/I2.fastq.gz -o $dir/smallfastqs/$Run.1.I2.fastq.gz -S 12000000 --thread 1 -d 4 -A -G -L -Q 2>>$dir/split.log &
	    wait
	    dosplitfull(){
                /mnt/users/sai/Package/fastp/fastp -i $2/R1.$1.fastq.gz -o $2/smallfastqs/$1.$3.R1.fastq.gz -S 12000000 --thread 1 -d 4 -A -G -L -Q 2>>$2/split.log &
                /mnt/users/sai/Package/fastp/fastp -i $2/R2.$1.fastq.gz -o $2/smallfastqs/$1.$3.R2.fastq.gz -S 12000000 --thread 1 -d 4 -A -G -L -Q 2>>$2/split.log &
                /mnt/users/sai/Package/fastp/fastp -i $2/I1.$1.fastq.gz -o $2/smallfastqs/$1.$3.I1.fastq.gz -S 12000000 --thread 1 -d 4 -A -G -L -Q 2>>$2/split.log &
                /mnt/users/sai/Package/fastp/fastp -i $2/I2.$1.fastq.gz -o $2/smallfastqs/$1.$3.I2.fastq.gz -S 12000000 --thread 1 -d 4 -A -G -L -Q 2>>$2/split.log &
		wait
            }
            export -f dosplitfull
            parallel --delay 1 dosplitfull {} $dir $Run ::: $(seq $nolane)
	elif [ "$Runtype" = QC ]; then
	    # Runing QC pipeline
	    echo "Split fastqs to small files"
	    # dosplitQC() works for fastqs that were not trimmed; it is fast
	    dosplitQC(){
		/mnt/users/sai/Package/fastp/fastp -i $2/R1.$1.fastq.gz -o $2/smallfastqs/$1.$3.R1.fastq.gz --reads_to_process $4 -S 2000000 --thread 1 -d 4 -A -G -L -Q 2>>$2/split.log &
		/mnt/users/sai/Package/fastp/fastp -i $2/R2.$1.fastq.gz -o $2/smallfastqs/$1.$3.R2.fastq.gz --reads_to_process $4 -S 2000000 --thread 1 -d 4 -A -G -L -Q 2>>$2/split.log &
		/mnt/users/sai/Package/fastp/fastp -i $2/I1.$1.fastq.gz -o $2/smallfastqs/$1.$3.I1.fastq.gz --reads_to_process $4 -S 2000000 --thread 1 -d 4 -A -G -L -Q 2>>$2/split.log &
		/mnt/users/sai/Package/fastp/fastp -i $2/I2.$1.fastq.gz -o $2/smallfastqs/$1.$3.I2.fastq.gz --reads_to_process $4 -S 2000000 --thread 1 -d 4 -A -G -L -Q 2>>$2/split.log &
		wait
            }
            export -f dosplitQC
	    # dosplitQC2() works for accidentially trimmed fastq, may have empty lines
	    dosplitQC2(){
		zcat $2/R1.$1.fastq.gz | head -n $4 | awk 'BEGIN { file = "1" } { print | "pigz --fast -p 18 >" file ".'$1'.'$3'.R1.fastq.gz" } NR % 2000000 == 0 {close("pigz --fast -p 18  > " file".'$1'.'$3'.R1.fastq.gz"); file = file + 1}' &
		zcat $2/R2.$1.fastq.gz | head -n $4 | awk 'BEGIN { file = "1" } { print | "pigz --fast -p 18 >" file ".'$1'.'$3'.R2.fastq.gz" } NR % 2000000 == 0 {close("pigz --fast -p 18  > " file".'$1'.'$3'.R2.fastq.gz"); file = file + 1}' &
		zcat $2/I1.$1.fastq.gz | head -n $4 | awk 'BEGIN { file = "1" } { print | "pigz --fast -p 18 >" file ".'$1'.'$3'.I1.fastq.gz" } NR % 2000000 == 0 {close("pigz --fast -p 18  > " file".'$1'.'$3'.I1.fastq.gz"); file = file + 1}' &
		zcat $2/I2.$1.fastq.gz | head -n $4 | awk 'BEGIN { file = "1" } { print | "pigz --fast -p 18 >" file ".'$1'.'$3'.I2.fastq.gz" } NR % 2000000 == 0 {close("pigz --fast -p 18  > " file".'$1'.'$3'.I2.fastq.gz"); file = file + 1}' &
		wait
	    }
            export -f dosplitQC2
	    let reads=12100000/$nolane
	    parallel --delay 1 dosplitQC {} $dir $Run $reads ::: $(seq $nolane)
#	    cd $dir/smallfastqs/
#	    parallel --delay 1 dosplitQC2 {} $dir $Run $reads ::: $(seq $nolane)
	else
            echo "Unknown sequencer type, exiting" && exit
        fi
    fi
    
    # trim and index fastq
    if [ -f $dir/fastp.json ]; then rm $dir/fastp.json $dir/fastp.html; fi
    ls $dir/smallfastqs | grep R1 > $dir/filesr1.xls
    ls $dir/smallfastqs | grep R2 > $dir/filesr2.xls
    ls $dir/smallfastqs | grep I1 > $dir/filesi1.xls
    ls $dir/smallfastqs | grep I2 > $dir/filesi2.xls
    cd $dir/
    if [ -f $dir/fastqs/Sub.0001.1.discard.R1.fq.gz ]; then
        echo "Found Sub.0001.1.discard.R1.fq.gz, skip updating index"
    else
        echo "Update index and trim fastqs"
	noreadfile=`ls $dir/smallfastqs | grep R1 2>/dev/null | wc -l`
	noindexfile=`ls $dir/smallfastqs | grep I1 2>/dev/null | wc -l`
	if [ $noreadfile == $noindexfile ]; then
	    paste filesr1.xls filesr2.xls filesi1.xls filesi2.xls | awk -v OFS='\t' '{print $1, $2, $3, $4, substr($1,1,7)}'> Filelist2.xls
	    parallel --jobs $cores --colsep '\t' 'if [ -f '$dir'/fastqs/Sub.{5}.discard.R1.fq.gz ]; then echo "found Sub.{5}.discard.R1.fq.gz"; \
	    	     	    	   	  else python3 '$myPATH'/fastq.process.py3.v0.6.py \
                                          -a '$dir'/smallfastqs/{1} -b '$dir'/smallfastqs/{2} \
                                          --c '$dir'/smallfastqs/{3} --d '$dir'/smallfastqs/{4} \
					  --out '$dir'/fastqs/Sub.{5} \
                                          -y '$yaml' && pigz --fast -p 4 '$dir'/fastqs/Sub.{5}*fq; fi' :::: Filelist2.xls
	else
	    paste filesr1.xls filesr2.xls | awk -v OFS='\t' '{print $1, $2, substr($1,1,7)}'> Filelist2.xls
	    parallel --jobs $cores --colsep '\t' 'if [ -f '$dir'/fastqs/Sub.{3}.discard.R1.fq.gz ]; then echo "found Sub.{3}.discard.R1.fq.gz"; \ 
	    	     	    	   	  else python3 '$myPATH'/fastq.process.py3.v0.6.py -a '$dir'/smallfastqs/{1} -b '$dir'/smallfastqs/{2} \
                                          --out '$dir'/fastqs/Sub.{3} \
					  -y '$yaml' && pigz --fast -p 4 '$dir'/fastqs/Sub.{3}*fq; fi' :::: Filelist2.xls
	    # --qc # option is available to process even smaller number of reads
	fi
    fi
    rm filesr1.xls filesr2.xls filesi1.xls filesi2.xls
    if [ -f Filelist2.xls ]; then
	rm Filelist2.xls
    fi
fi 

# merge fastq
echo "Merge fastqs"
parallel --jobs 4 'if [ -f '$dir'/fastqs/{}.R1.fastq.gz ]; then echo "found {}.R1.fastq.gz"; \
	 	      	   else ls '$dir'/fastqs/Sub*{}*R1.fq.gz | xargs cat > '$dir'/fastqs/{}.R1.fastq.gz; fi' ::: ${Project[@]} discard
parallel --jobs 4 'if [ -f '$dir'/fastqs/{}.R2.fastq.gz ]; then echo "found {}.R2.fastq.gz"; \
                           else ls '$dir'/fastqs/Sub*{}*R2.fq.gz | xargs cat > '$dir'/fastqs/{}.R2.fastq.gz; fi' ::: ${Project[@]} discard
