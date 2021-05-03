index=0
for Name in ${Project[@]}; do
    echo "project $index : $Name" 
    if [ -d $dir/$Name ]; then
	echo "Found $Name dir, skip this project"
    else    
	if [ ${Genomes[$index]} == "both" ]; then
	    if [ ${Type[$index]} == "ATAC" ] || [ ${Type[$index]} == "DipC" ] || [ ${Type[$index]} == "TAPS" ]; then
		Genome1=(hg19 mm10) # Genome1 for alignment, Genome2 for other steps
	    else
		Genome1=(both) # for RNA and cellhash
	    fi
            Genome2=(hg19 mm10)
	else
	    Genome1=${Genomes[$index]}
	    Genome2=${Genomes[$index]}
	fi
	
	for Species in ${Genome1[@]}; do
	    if [ -d $dir/$Name ]; then
		echo "Found $Name folder, skip alignment"
	    else
		echo "Align $Name to $Species ${Type[$index]} library"
		cd $dir/fastqs/
		if [ -f $dir/fastqs/$Name.$Species.bam ]; then
		    echo "Found $Name.$Species.bam, skip alignment"
		elif [ ${Type[$index]} == "ATAC" ]; then
		    (bowtie2 -X2000 -p $cores --rg-id $Name \
			     -x $bowtieGenome/$Species/$Species \
			     -1 $dir/fastqs/$Name.R1.fastq.gz \
			     -2 $dir/fastqs/$Name.R2.fastq.gz | \
			 samtools view -bS -@ $cores - -o $Name.$Species.bam) 2>$Name.$Species.align.log
		elif [ ${Type[$index]} == "crop" ]; then
                    (bowtie2 -X2000 -p $cores --rg-id $Name \
			     -N 0 -L 30 \
                             -x $bowtieGenome/viralRNA/viralRNA \
			     -U $dir/fastqs/$Name.R1.fastq.gz | \
                         samtools view -bS -@ $cores - -o $Name.$Species.bam) 2>$Name.$Species.align.log
		elif [ ${Type[$index]} == "cite" ]; then
                    (bowtie2 -X2000 -p $cores --rg-id $Name \
                             -N 0 -L 30 \
                             -x $bowtieGenome/cite/cite \
                             -U $dir/fastqs/$Name.R1.fastq.gz | \
                         samtools view -bS -@ $cores - -o $Name.$Species.bam) 2>$Name.$Species.align.log
		elif [ ${Type[$index]} == "cellhash" ]; then
                    echo "Match CITE-seq barcodes, convert fastq to bam"
                    (bowtie2 -p $cores --rg-id $Name \
                             -x $bowtieGenome/cellhash/$Species.cellhash/$Species.cellhash \
                             -U $dir/fastqs/$Name.R1.fastq.gz | \
                         samtools view -bS - -o $Name.$Species.bam) 2>$Name.$Species.align.log
                elif [ ${Type[$index]} == "TAPS" ]; then
		    (bismark --parallel 4  \
			     --multicore 4 \
			     --rg_id $Name \
			     --non_directional \
                             --genome $bismarkGenome/$Species \
			     -X 2000 \
                             -1 $dir/fastqs/$Name.R1.fastq.gz \
                             -2 $dir/fastqs/$Name.R2.fastq.gz) 1> $Name.$Species.bismark.log 2>>$Name.$Species.bismark.log
		    mv $Name.R1_bismark_bt2_pe.bam $Name.$Species.bam
		    mv $Name.R1_bismark_bt2_PE_report.txt $Name.$Species.align.log			
		elif [ ${Type[$index]} == "DipC" ]; then
		    ## align and process with HiC-pro
		    ## without removing duplicates
		    if [ -d $dir/fastqs/$Name.$Species.processed/ ]; then
			echo "found $Name.$Species.processed, skip HiC-pro"
		    else
			echo "running HiC-pro pipeline"
			~/Package/HiCpro/HiC-Pro_2.11.4/bin/HiC-Pro -c $myPATH/Config_MboI_$Species.txt -i ./$Name.rawdata/ -o $Name.$Species.processed/
		    fi
		elif [ ${Type[$index]} == "RNA" ]; then
		    if [ -f $dir/fastqs/$Name.$Species.align.log ]; then
			echo "found $dir/fastqs/$Name.$Species.align.log, skip alignment for UMI reads"
		    else
			echo "Align UMI reads"
			STAR --chimOutType WithinBAM \
			     --runThreadN $cores \
			     --genomeDir $starGenome/$Species/ \
			     --readFilesIn $dir/fastqs/$Name.R1.fastq.gz  \
			     --outFileNamePrefix $dir/fastqs/$Name.$Species. \
			     --outFilterMultimapNmax 20 \
			     --outFilterMismatchNoverLmax 0.06 \
			     --limitOutSJcollapsed 2000000 \
			     --outSAMtype BAM Unsorted \
			     --limitIObufferSize 400000000 \
			     --outReadsUnmapped Fastx \
			     --readFilesCommand zcat
		    
			mv $dir/fastqs/$Name.$Species.Aligned.out.bam $dir/fastqs/$Name.$Species.bam
			rm -r *_STARtmp *Log.progress.out *SJ.out.tab *Unmapped.out.mate* *_STARpass1 *_STARgenome  
			mv $dir/fastqs/$Name.$Species.Log.final.out $dir/fastqs/$Name.$Species.align.log
		    fi
		else
		    echo "Unknown parameter"
		    exit
		fi
		if [ -f $dir/fastqs/$Name.$Species.st.bam ]; then
		    echo "Found $Name.$Species.st.bam, skip sorting bam"
		else
		    echo "Sort $Name.$Species.bam"
		    cd $dir/fastqs/
		    if [ $mode == "regular" ]; then
			java -Xmx24g -Djava.io.tmpdir=$dir/temp/ -jar $picardPATH SortSam SO=coordinate I=$Name.$Species.bam O=$Name.$Species.st.bam VALIDATION_STRINGENCY=SILENT TMP_DIR=$dir/temp/ 2>>$dir/Run.log
                        samtools index -@ $cores $Name.$Species.st.bam
                       rm $Name.$Species.bam
                    else
			samtools sort -@ $cores -m 2G $Name.$Species.bam > $Name.$Species.st.bam
			samtools index -@ $cores $Name.$Species.st.bam
			rm $Name.$Species.bam
		    fi
		fi
				
		# Update RGID's and Headers
		if [ -f $dir/fastqs/$Name.$Species.rigid.reheader.st.bam.bai ]; then
		    echo "Found $Name.$Species.rigid.reheader.st.bam, skip update RGID"
		elif [ ${Type[$index]} == "TAPS" ]; then
		    if [ -f $Name.$Species.rigid.reheader.st.bam ]; then
                        echo "Found $Name.$Species.rigid.reheader.st.bam, skip update RG tag"
                    else
                        echo "Update RGID for $Name.$Species.st.bam"
                        samtools view -H $Name.$Species.st.bam > $Name.$Species.st.header.sam
                        samtools view -@ $cores $Name.$Species.st.bam | cut -f1 | \
                            sed 's/_/\t/g' | cut -f2 | sort --parallel=$cores -S 10G | \
                            uniq | awk -v OFS='\t' '{print "@RG", "ID:"$1, "SM:Barcode"NR}' > header.temp.sam
                        sed -e '/\@RG/r./header.temp.sam' $Name.$Species.st.header.sam > $Name.$Species.rigid.st.header.sam
                        cat $Name.$Species.rigid.st.header.sam <(samtools view $Name.$Species.st.bam | \
                                                                     awk -v OFS='\t' '{$1=substr($1,1,length($1)-23)""substr($1,length($1)-22,23); print $0}') |
                            samtools view -@ $cores -bS > $Name.$Species.rigid.reheader.st.bam
                        samtools index -@ $cores $Name.$Species.rigid.reheader.st.bam
                        rm header.temp.sam $Name.$Species.rigid.st.header.sam
                    fi
		elif [ ${Type[$index]} == "RNA" ] || [ ${Type[$index]} == "crop" ] || [ ${Type[$index]} == "cite" ] || [ ${Type[$index]} == "cellhash" ]; then
		    if [ -f $Name.$Species.rigid.reheader.unique.st.bam.bai ]; then
                        echo "Found $Name.$Species.rigid.reheader.st.bam, skip update RG tag"
		    else
			echo "Update RGID for $Name.$Species.st.bam"
			echo "Remove low quality reads and unwanted chrs"
			samtools view -H $Name.$Species.st.bam | sed 's/chrMT/chrM/g' > $Name.$Species.st.header.sam
			if [ ${keepMultiMapping[$index]} == "T" ]; then
			    # keep primary aligned reads only
			    # modify chrMT to chrM
			    cat $Name.$Species.st.header.sam <(samtools view -@ $cores $Name.$Species.st.bam | \
				   sed 's/chrMT/chrM/g' | \
				   awk -v OFS='\t' '{$1=substr($1,1,length($1)-34)""substr($1,length($1)-22,23)""substr($1,length($1)-34,11); print $0}') |
				samtools view -@ $cores -bS -F 256 > $Name.$Species.rigid.reheader.unique.st.bam
			else
                            cat $Name.$Species.st.header.sam <(samtools view -@ $cores $Name.$Species.st.bam | \
				   sed 's/chrMT/chrM/g'	| \
				   awk -v OFS='\t' '{$1=substr($1,1,length($1)-34)""substr($1,length($1)-22,23)""substr($1,length($1)-34,11); print $0}') |
				samtools view -@ $cores -bS -q 30  > $Name.$Species.rigid.reheader.unique.st.bam
			fi
			samtools index -@ $cores $Name.$Species.rigid.reheader.unique.st.bam
			rm $Name.$Species.st.header.sam
		    fi
		fi
		
		# ATAC processing convert bam to bed
		if [ ${Type[$index]} == "ATAC" ]; then
		    if [ -f $Name.$Species.rmdup.bam ]; then
			echo "Skip processing ATAC bam"
		    else
			if [ -f $Name.$Species.namesort.bam ]; then
		            echo "Skip sort bam on name"
			else
			    echo "Sort bam on name"
			    echo "Remove low quality reads, unwanted chrs & namesort" $Name.$Species.st.bam
			    chrs=`samtools view -H $Name.$Species.st.bam | grep chr | cut -f2 | sed 's/SN://g' | grep -v chrM | grep -v Y | awk '{if(length($0)<6)print}'`
			    samtools view -b -q 30 -f 0x2 $Name.$Species.st.bam  `echo $chrs` | samtools sort -@ $cores -m 2G -n -o $Name.$Species.namesort.bam
			fi
			if [ -f $Name.$Species.bed.gz ]; then
			    echo "Skip converting bam to bed.gz"
			else
                            echo "Convert rmdup.bam to bed.gz & mark duplicates"
                            bedtools bamtobed -i $Name.$Species.namesort.bam -bedpe | \
				sed 's/_/\t/g' | \
				awk -v OFS="\t" '{if($10=="+"){print $1,$2+4,$6+4,$8}else if($10=="-"){print $1,$2-5,$6-5,$8}}' |\
				sort --parallel=$cores -S 40G  -k4,4 -k1,1 -k2,2 -k3,3 | \
				uniq -c | \
				awk -v OFS="\t" '{print $2, $3, $4, $5, $1}' | \
				pigz --fast -p $cores > $Name.$Species.bed.gz
			    rm $Name.$Species.namesort.bam
			fi
			# convert to a bam file for QC
                        bedToBam -i <(zcat $Name.$Species.bed.gz) -g /mnt/users/sai/Script/$Species.chrom.sizes | samtools sort -@ $cores -m 2G - > $Name.$Species.rmdup.bam
                        samtools index -@ $cores $Name.$Species.rmdup.bam
                    fi
		fi
			
		# Update RGID's and Headers for dipc bedpe
		## this part needs to be updated
		if [ ${Type[$index]} == "DipC" ] ; then
                    if [ -f $dir/fastqs/$Name.$Species/$Name.$Species.rg.allPairs ]; then
			echo "Found $Name.$Species.rg.allPairs, skip update RGID"
			echo "Update RGID for $Name.$Species.allValidPairs"
			cd $dir/fastqs/
			if [ "$Sequencer" = Novaseq ]; then
                            /usr/bin/python $myPATH/updateRGID_dipC_novaseq_V1.py --input $dir/fastqs/$Name.$Species.processed/hic_results/data/$Name/$Name.allValidPairs --out $Name.$Species.rg.allValidPairs
			else
                            /usr/bin/python $myPATH/updateRGID_dipC_nextseq_V1.py --input $dir/fastqs/$Name.$Species.processed/hic_results/data/$Name/$Name.allValidPairs --out $Name.$Species.rg.allValidPairs
			    tempdir=$dir/fastqs/$Name.$Species.processed/hic_results/data/$Name/
			    cat $tempdir/$Name\_*DEPairs $tempdir/$Name\_*.DumpPairs $tempdir/$Name\_*.FiltPairs $tempdir/$Name\_*.REPairs \
				$tempdir/$Name\_*.SCPairs $tempdir/$Name\_*.SinglePairs $tempdir/$Name\_*.validPairs > $dir/fastqs/$Name.$Species.allPairs
			    /usr/bin/python $myPATH/updateRGID_dipC_nextseq_V1.py --input $dir/fastqs/$Name.$Species.allPairs --out $Name.$Species.rg.allPairs
			fi
		    fi
                
		    if [ -f $dir/fastqs/$Name.$Species.rg.allValidPairs.hic ]; then
			echo "found $Name.$Species.rg.allValidPairs, skip convert to .hic"
		    else
			echo "convert Valid pairs to hic"
			~/Package/HiCpro/HiC-Pro_2.11.4/bin/utils/hicpro2juicebox.sh \
			    -i $Name.$Species.rg.allValidPairs -g $Species \
			    -j /mnt/users/sai/Package/juicer/CPU/common/juicer_tools.1.9.9_jcuda.0.8.jar
		    fi
		fi
	        
		# remove dups, mito and secondary alignment
		if [ -f $dir/fastqs/$Name.$Species.rigid.reheader.unique.st.bam ] || [ -f $dir/fastqs/$Name.$Species.wdup.all.bam ]; then
		    echo "Found $Name.$Species.rigid.reheader.nomito.st.bam, skip mito removal"
		else
		    if [ ${Type[$index]} == "TAPS" ]; then
			echo "Remove unwanted chrs:" $Name.$Species.rigid.reheader.st.bam
			chrs=`samtools view -H $Name.$Species.st.bam | grep chr | cut -f2 | sed 's/SN://g' | grep -v chrM | grep -v Y | awk '{if(length($0)<6)print}'`
			samtools view -@ $cores -b $Name.$Species.rigid.reheader.st.bam -o $Name.$Species.wdup.all.bam `echo $chrs`
			samtools index -@ $cores $Name.$Species.wdup.all.bam
		    fi
		fi

		if [ ${Type[$index]} == "DipC" ]; then
		    cd $dir/fastqs/
                    # count unfiltered reads
                    if [ ! -f $Name.$Species.unfiltered.counts.csv ]; then
                        echo "Count unfiltered reads" $Name.$Species.rg.allPairs
                        if [ ! -f $Name.$Species.wdup.RG.bed ]; then
			    cat $Name.$Species.rg.allPairs | cut -f13 > $Name.$Species.wdup.RG.bed
                        fi
                        if [ ! -f $Name.$Species.wdup.RG.freq.bed ]; then
                            cat $Name.$Species.wdup.RG.bed | sort --parallel=$cores -S 24G | uniq -c > $Name.$Species.wdup.RG.freq.bed
                        fi
                        Rscript $myPATH/sum_reads.R $dir/fastqs/ $Name.$Species.wdup.RG.freq.bed --save
                        mv $Name.$Species.wdup.RG.freq.bed.csv $Name.$Species.unfiltered.counts.csv
                        rm $Name.$Species.wdup.RG.bed $Name.$Species.wdup.RG.freq.bed
                    fi

                    # remove barcode combination that has less then N unfiltered reads
                    if [ -f $dir/fastqs/$Name.$Species.rg.cutoff.allPairs ]; then
                        echo "Skip removing barcode combinations"
                    else
                        echo "Remove barcode combinations that have less then ${ReadsPerBarcode[$index]} reads"
                        sed -e 's/,/\t/g' $Name.$Species.unfiltered.counts.csv | awk -v OFS=',' 'NR>=2 {if($5 >= '${ReadsPerBarcode[$index]}') print $1,$2,$3,$4}'  > $Name.$Species.barcodes.txt
                        grep -wFf $Name.$Species.barcodes.txt $Name.$Species.rg.allPairs | awk -v OFS='\t' '{print $9,$10,$13}' | \
			    awk -v OFS='\t' '{if($1 >=  $2) print $2, $1, $3; else print $1, $2, $3}' > $Name.$Species.rg.cutoff.allPairs
			grep -wFf $Name.$Species.barcodes.txt $Name.$Species.rg.allValidPairs | awk -v OFS='\t' '{print $9,$10,$13}' | \
                            awk -v OFS='\t' '{if($1 >=  $2) print $2, $1, $3; else print $1, $2, $3}' > $Name.$Species.rg.cutoff.allValidPairs
		    fi

		    # remove duplicates
                    if [ -f $dir/fastqs/$Name.$Species.rg.cutoff.st.allValidPairs ]; then
                        echo "Skip split $Name.$Species.rg.cutoff.allValidPairs and mark duplicate"
                    else
                        echo "Remove duplicates by cell barcode and cutting sites"
			cat $Name.$Species.rg.cutoff.allPairs | sort -k1,1 -k2,2 -k3,3 | uniq -c | \
			    awk -v OFS='\t' '{print $2,$3,$4,$1}' > $Name.$Species.rg.cutoff.st.allPairs
			cat $Name.$Species.rg.cutoff.allValidPairs | sort -k1,1 -k2,2 -k3,3 | uniq -c | \
                            awk -v OFS='\t' '{print $2,$3,$4,$1}' > $Name.$Species.rg.cutoff.st.allValidPairs
			rm $Name.$Species.rg.cutoff.allPairs $Name.$Species.rg.cutoff.allValidPairs
                    fi
                fi
	       
		# count reads for ATAC
		if [ ${Type[$index]} == "ATAC" ] ; then
		    if [ -f $Name.$Species.filtered.counts.csv ] || [ -f $Name.$Species.filtered.counts.csv.gz ]; then
			echo "found $Name.$Species.filtered.counts.csv"
		    else
			echo "count unfiltered reads"
			zcat $Name.$Species.bed.gz | \
			    awk -v OFS='\t' '{a[$4] += $5} END{for (i in a) print a[i], i}'| \
			    awk -v OFS='\t' '{if($1 >= '${ReadsPerBarcode[$index]}') print }'> $Name.$Species.wdup.RG.freq.bed
			Rscript $myPATH/sum_reads.R $dir/fastqs/ $Name.$Species.wdup.RG.freq.bed --save
			mv $Name.$Species.wdup.RG.freq.bed.csv $Name.$Species.unfiltered.counts.csv

			echo "count filtered reads"
			zcat $Name.$Species.bed.gz | cut -f4 | uniq -c | \
			    awk -v OFS='\t' '{if($1 >= '${ReadsPerBarcode[$index]}') print }' > $Name.$Species.rmdup.RG.freq.bed
			Rscript $myPATH/sum_reads.R $dir/fastqs/ $Name.$Species.rmdup.RG.freq.bed --save
			mv $Name.$Species.rmdup.RG.freq.bed.csv $Name.$Species.filtered.counts.csv

			rm $Name.$Species.wdup.RG.freq.bed $Name.$Species.rmdup.RG.freq.bed 
		    fi
		fi

		# remove barcode with low counts from the fragment file for ATAC
                if [ ${Type[$index]} == "ATAC" ] ; then
		    if [ -f $Name.$Species.cutoff.bed.gz ]; then
			echo "Skip removing low counts barcode combination"
                    else
			echo "Remove low counts barcode combination"
			sed -e 's/,/\t/g' $Name.$Species.filtered.counts.csv | \
                            awk -v OFS=',' 'NR>=2 {if($5 >= '${ReadsPerBarcode[$index]}') print $1,$2,$3,$4} '  > $Name.$Species.barcodes.txt
			grep -wFf $Name.$Species.barcodes.txt  <(zcat $Name.$Species.bed.gz) | pigz --fast -p $cores  > $Name.$Species.cutoff.bed.gz
		    fi
		fi
		# count reads for TAPS
		if [ ${Type[$index]} == "TAPS" ]; then
		    # count unfiltered reads
		    if [ ! -f $Name.$Species.unfiltered.counts.csv ] || [ -f $Name.$Species.unfiltered.counts.csv.gz ]; then
			echo "Count unfiltered reads" $Name.$Species.wdup.all.bam
			if [ ! -f $Name.$Species.wdup.RG.bed ]; then
			    samtools view -@ $cores $Name.$Species.wdup.all.bam | cut -f1 |  sed 's/_/\t/g' | cut -f2 > $Name.$Species.wdup.RG.bed
			fi
			if [ ! -f $Name.$Species.wdup.RG.freq.bed ]; then
			    cat $Name.$Species.wdup.RG.bed | sort --parallel=$cores -S 24G | uniq -c > $Name.$Species.wdup.RG.freq.bed
			fi
			Rscript $myPATH/sum_reads.R $dir/fastqs/ $Name.$Species.wdup.RG.freq.bed --save
			mv $Name.$Species.wdup.RG.freq.bed.csv $Name.$Species.unfiltered.counts.csv
			rm $Name.$Species.wdup.RG.bed $Name.$Species.wdup.RG.freq.bed
		    fi
		    
		    # remove barcode combination that has less then N unfiltered reads
		    if [ -f $dir/fastqs/$Name.$Species.wdup.bam.bai ]; then
			echo "Skip removing barcode combinations"
		    else
			echo "Remove barcode combinations that have less then ${ReadsPerBarcode[$index]} reads"
			sed -e 's/,/\t/g' $Name.$Species.unfiltered.counts.csv | awk -v OFS=',' 'NR>=2 {if($5 >= '${ReadsPerBarcode[$index]}') print $1,$2,$3,$4}'  > $Name.$Species.barcodes.txt
			# too many barcode combinations will result in picard failure
			samtools view -@ $cores -R $Name.$Species.barcodes.txt $Name.$Species.wdup.all.bam -o $Name.$Species.wdup.bam
			samtools index -@ $cores $Name.$Species.wdup.bam	
		    fi
		    # remove duplicates
		    if [ -f $dir/fastqs/$Name.$Species.rmdup.bam ]; then
			echo "Skip split $Name.$Species.wdup.bam and mark duplicate"
		    else
			echo "Mark duplicates: $Name.$Species.wdup.bam"
			samtools sort -@ $cores -m 2G -n -o $Name.$Species.wdup.namest.bam $Name.$Species.wdup.bam
			(deduplicate_bismark --bam $Name.$Species.wdup.namest.bam -o $Name.$Species) 1>> $Name.$Species.bismark.log 2>> $Name.$Species.bismark.log
			java -Xmx24g -Djava.io.tmpdir=$dir/temp/ -jar $picardPATH SortSam SO=coordinate I=$Name.$Species.deduplicated.bam O=$Name.$Species.rmdup.bam VALIDATION_STRINGENCY=SILENT TMP_DIR=$dir/temp/ 2>>$dir/Run.log
			samtools index -@ $cores $Name.$Species.rmdup.bam
			rm $Name.$Species.wdup.namest.bam
		    fi
		fi
		
		# RNA processing
		if [ ${Type[$index]} == "RNA" ] && [ ! -f $Name.mm10.rigid.reheader.unique.st.bam ]; then
		    cd $dir/fastqs/
		    if [ ${Genomes[$index]} == "both" ]; then
			echo "Split into hg and mm"
			chrs1=`samtools view -H $Name.$Species.rigid.reheader.unique.st.bam | grep hg19 | cut -f2 | sed 's/SN://g' | awk '{if(length($0)<8)print}'`
			chrs2=`samtools view -H $Name.$Species.rigid.reheader.unique.st.bam | grep mm10 | cut -f2 | sed 's/SN://g' | awk '{if(length($0)<8)print}'`
			samtools view -@ $cores -b $Name.$Species.rigid.reheader.unique.st.bam -o $Name.$Species.temp1.bam `echo ${chrs1[@]}`
			samtools view -@ $cores -b $Name.$Species.rigid.reheader.unique.st.bam -o $Name.$Species.temp2.bam `echo ${chrs2[@]}`
			samtools view -@ $cores -h $Name.$Species.temp1.bam | sed 's/hg19_/chr/g' | sed 's/chrMT/chrM/g'| samtools view -@ $cores -b -o $Name.hg19.rigid.reheader.unique.st.bam
			samtools view -@ $cores -h $Name.$Species.temp2.bam | sed 's/mm10_/chr/g' | sed 's/chrMT/chrM/g'| samtools view -@ $cores -b -o $Name.mm10.rigid.reheader.unique.st.bam

			samtools index -@ $cores $Name.hg19.rigid.reheader.unique.st.bam &
			samtools index -@ $cores $Name.mm10.rigid.reheader.unique.st.bam &
			wait
			rm $Name.$Species.temp1.bam $Name.$Species.temp2.bam
		    else
			echo "Single species is aligned"
		    fi		
		fi

		# cellhashing processing
		if [ ${Type[$index]} == "cellhash" ] && [ ! -f $Name.$Species.rigid.reheader.unique.st.bam ]; then
		    cd $dir/fastqs/
                    echo "Split into hg and mm"
                    chrs1=`samtools view -H $Name.$Species.rigid.reheader.unique.st.bam | grep hg19 | cut -f2 | sed 's/SN://g'`
                    chrs2=`samtools view -H $Name.$Species.rigid.reheader.unique.st.bam | grep mm10 | cut -f2 | sed 's/SN://g'`
                    samtools view -@ $cores -b $Name.$Species.rigid.reheader.unique.st.bam -o temp1.bam `echo ${chrs1[@]}`
                    samtools view -@ $cores -b $Name.$Species.rigid.reheader.unique.st.bam -o temp2.bam `echo ${chrs2[@]}`
                    mv temp1.bam $Name.hg19.rigid.reheader.unique.st.bam
                    mv temp2.bam $Name.mm10.rigid.reheader.unique.st.bam
                    samtools index -@ $cores $Name.hg19.rigid.reheader.unique.st.bam &
                    samtools index -@ $cores $Name.mm10.rigid.reheader.unique.st.bam &
		    wait
                    rm temp1.bam temp2.bam
		fi

		# assign feature to reads
		cd  $dir/fastqs/
		if [ ${Type[$index]} == "RNA" ]; then
		    for Species2 in ${Genome2[@]}; do
			if [ -f $dir/fastqs/$Name.$Species2.exon.featureCounts.bam ]; then
			    echo "Skip exon feasure count"
			else
		    	    # excliude multimapping, uniquely mapped reads only, -Q 30, for real sample, might consider include multi-mapping
			    echo "Feature counting on exons"
			    # count exon
                            if [ ${keepMultiMapping[$index]} == "T" ]; then
				featureCounts -T $cores -Q 0 -M -a $myPATH/gtf/$Species2.$refgene.gtf -t exon -g $genename \
                                    -o $Name.$Species2.feature.count.txt -R BAM $Name.$Species2.rigid.reheader.unique.st.bam >>$dir/Run.log
                            else
				featureCounts -T $cores -Q 30 -a $myPATH/gtf/$Species2.$refgene.gtf -t exon -g $genename \
                                    -o $Name.$Species2.feature.count.txt -R BAM $Name.$Species2.rigid.reheader.unique.st.bam >>$dir/Run.log
                            fi
			    # Extract reads that assigned to genes
			    mv $Name.$Species2.rigid.reheader.unique.st.bam.featureCounts.bam $Name.$Species2.exon.featureCounts.bam
			fi
			if [ -f $dir/fastqs/$Name.$Species2.wdup.bam.bai ]; then
			    echo "Skip intron and exon feasure count"
			else
			    # count both intron and exon
			    echo "Count feature on both intron and exon"
			    if [ $keepIntron == "T" ]; then
				if [ ${keepMultiMapping[$index]} == "T" ]; then
				    featureCounts -T $cores -Q 0 -M -a $myPATH/gtf/$Species2.$refgene.gtf -t gene -g $genename \
						  -o $Name.$Species2.feature.count.txt -R BAM $Name.$Species2.exon.featureCounts.bam >>$dir/Run.log
				# then for reads mapped to multiple genes,  only keep if all its alignments are within a single gene
				else
				    featureCounts -T $cores -Q 30 -a $myPATH/gtf/$Species2.$refgene.gtf -t gene -g $genename \
					-o $Name.$Species2.feature.count.txt -R BAM $Name.$Species2.exon.featureCounts.bam >>$dir/Run.log
				fi
				samtools sort -@ $cores -m 2G -o $Name.$Species2.wdup.bam  $Name.$Species2.exon.featureCounts.bam.featureCounts.bam
				rm $Name.$Species2.exon.featureCounts.bam.featureCounts.bam
			    else
				samtools sort -@ $cores -m 2G -o $Name.$Species2.wdup.bam $Name.$Species2.exon.featureCounts.bam
			    fi
			    samtools index -@ $cores $Name.$Species2.wdup.bam
			fi
		    done
		fi
		# add crop barcode to reads
                if [ ${Type[$index]} == "crop" ] || [ ${Type[$index]} == "cite" ] || [ ${Type[$index]} == "cellhash" ]; then
                    for Species2 in ${Genome2[@]}; do
                        if [ -f $dir/fastqs/$Name.$Species2.gene.wdup.bam.bai ]; then
                            echo "Skip adding aligned guide barcode to reads"
                        else
			    samtools view -H $Name.$Species2.rigid.reheader.unique.st.bam > $Name.$Species2.header.sam
			    samtools view $Name.$Species2.rigid.reheader.unique.st.bam | awk '{print $0"\tXT:Z:"$3}' > $Name.$Species2.temp.sam
			    cat $Name.$Species2.header.sam $Name.$Species2.temp.sam | samtools view -b -@ $cores -o $Name.$Species2.gene.wdup.bam
                            samtools index $Name.$Species2.gene.wdup.bam
                            rm $Name.$Species2.header.sam $Name.$Species2.temp.sam
			fi
		    done
		fi
		
		# group UMIs
		if [ ${Type[$index]} == "RNA" ] || [ ${Type[$index]} == "cellhash" ] || [ ${Type[$index]} == "crop" ] || [ ${Type[$index]} == "cite" ]; then
		    for Species2 in ${Genome2[@]}; do
			if [ -f $dir/fastqs/$Name.$Species2.grouped.bam ] || [ -f $dir/fastqs/$Name.$Species2.groups.tsv.gz ] || [ -f $dir/fastqs/$Name.$Species2.groups.tsv ]; then
			    echo "Found $Name.$Species.groups.bam, skip grouping UMIs"
			else
			    echo "Group reads to unique UMIs"
			    if [ $mode == "regular" ]; then
				## this is previously used. Seems to get more slant and fewer UMIs, but get accurate lib size estimation
				umi_tools group --extract-umi-method=read_id \
                                          --per-gene --gene-tag=XT --per-cell \
                                          -I $dir/fastqs/$Name.$Species2.wdup.bam \
                                          --output-bam -S $dir/fastqs/$Name.$Species2.grouped.bam \
                                          --group-out=$dir/fastqs/$Name.$Species2.groups.tsv --skip-tags-regex=Unassigned >>$dir/Run.log
			    else
				## own UMI dedup by matching bc-umi-align position
                                samtools view -@ $cores $Name.$Species2.wdup.bam | grep XT:Z: | \
                                    sed 's/Unassigned_Ambiguity/discard/g' | \
                                    sed 's/Unassigned_MappingQuality/discard/g' | \
                                    sed 's/\_/\t/g' | \
                                    awk -v OFS='\t' '{if($NF ~/discard/){$NF=$(NF-1)} print $5, $6, $2, $3, $NF}' | \
                                    sed 's/XT:Z://g' > $Name.$Species2.wdup.bed
                                # remove dup reads
                                python3 $myPATH/rm_dup_barcode_UMI_v3.py \
                                        -i $Name.$Species2.wdup.bed -o $Name.$Species2.groups.tsv --m 1
#                                pigz --fast -p $cores $dir/fastqs/$Name.$Species2.groups.tsv
                                rm $Name.$Species2.wdup.bed
			    fi
			fi
			# convert groupped UMI to bed file
			if [ -f $dir/fastqs/$Name.$Species2.bed.gz ]; then
			    echo "Skip removing single-read UMIs and GGGG"
			else
			    if [ $removeSingelReadUMI == "T" ]; then
				echo "Filter UMI that has only 1 read and convert to bed file"
				if [ $mode == "regular" ]; then
				    less $dir/fastqs/$Name.$Species2.groups.tsv | \
                                       	sed 's/_/\t/g' | awk 'FNR>1 {if($10 >1){if($9 != "GGGGGGGGGG"){print}}}' | \
                                        awk 'NR==1{ id="N"} {if(id != $11 ) {id = $11; print $2, $6, $10}}' | \
                                       	awk -v OFS="\t" 'NR==1{ t1=$1;t2=$2;readsum=0; umisum=0} {if(t1==$1 && t2==$2) {readsum+=$3; umisum+=1} \
                                            else {print t1, t2, umisum, readsum; t1=$1;t2=$2;umisum=1;readsum=$3}}' | \
                                        pigz --fast -p $cores > $Name.$Species2.bed.gz
				else
				    less $dir/fastqs/$Name.$Species2.groups.tsv | \
					awk -v OFS="\t" '{if($3 >1){print}}' |
                                        awk -v OFS="\t" 'NR==1{ t1=$1;t2=$2; t3=$3} {if(t1==$1 && t2==$2) {readsum+=$3} else {print t1, t2, t3, readsum; t1=$1;t2=$2;t3=$3;readsum=$3}}' | \
                                        pigz --fast -p $cores > $Name.$Species2.bed.gz
				fi
			    else
				echo "Convert to bed file"
				if [ $mode == "regular" ]; then
				    ## count reads per gene
				    less $dir/fastqs/$Name.$Species2.groups.tsv | \
					sed 's/_/\t/g' | awk 'FNR>1 {if($9 != "GGGGGGGGGG"){print}}' | \
					awk 'NR==1{ id="N"} {if(id != $11 ) {id = $11; print $2, $6, $10}}' | \ 
				        awk -v OFS="\t" 'NR==1{ t1=$1;t2=$2;readsum=0; umisum=0} {if(t1==$1 && t2==$2) {readsum+=$3; umisum+=1} \
					    else {print t1, t2, umisum, readsum; t1=$1;t2=$2;umisum=1;readsum=$3}}' | \	   
					pigz --fast -p $cores > $Name.$Species2.bed.gz
					## 3rd column is umi, 4th column is read
					## note: difference between these two versions of groups.tsv
					## 1) umitools output keep all the barcode-UMI and don't collapse them
					## 2) my script already collpsed them at alignment position level
				else
				    less $dir/fastqs/$Name.$Species2.groups.tsv | \
                                        awk -v OFS="\t" 'NR==1{ t1=$1;t2=$2; t3=$3} {if(t1==$1 && t2==$2) {readsum+=$3} else {print t1, t2, t3, readsum; t1=$1;t2=$2;t3=$3;readsum=$3}}' | \
                                        pigz --fast -p $cores > $Name.$Species2.bed.gz
				fi
			    fi
			fi

			# count reads for RNA
			if [ -f $Name.$Species2.filtered.counts.csv ] || [ -f $Name.$Species.filtered.counts.csv.gz ]; then
                            echo "found $Name.$Species.filtered.counts.csv"
			else
                            echo "count unfiltered reads"
			    zcat $Name.$Species2.bed.gz | \
				awk -v OFS='\t' '{a[$1] += $4} END{for (i in a) print a[i], i}' | \
				awk -v OFS='\t' '{if($1 >= '${ReadsPerBarcode[$index]}') print }'> $Name.$Species2.wdup.RG.freq.bed
                            Rscript $myPATH/sum_reads.R $dir/fastqs/ $Name.$Species2.wdup.RG.freq.bed --save
                            mv $Name.$Species2.wdup.RG.freq.bed.csv $Name.$Species2.unfiltered.counts.csv

                            echo "count filtered reads"
                            zcat $Name.$Species2.bed.gz | \
				awk -v OFS='\t' '{a[$1] += $3} END{for (i in a) print a[i], i}' | \
				awk -v OFS='\t' '{if($1 >= '${ReadsPerBarcode[$index]}') print }' > $Name.$Species2.rmdup.RG.freq.bed
                            Rscript $myPATH/sum_reads.R $dir/fastqs/ $Name.$Species2.rmdup.RG.freq.bed --save
                            mv $Name.$Species2.rmdup.RG.freq.bed.csv $Name.$Species2.filtered.counts.csv

                            rm $Name.$Species2.wdup.RG.freq.bed $Name.$Species2.rmdup.RG.freq.bed
			fi		
			
			# remove barcode combination that has less then N reads
			if [ -f $Name.$Species2.cutoff.bed.gz ]; then
			    echo "Skip removing low counts barcode combination"
			else
			    echo "Remove low counts barcode combination"
			    sed -e 's/,/\t/g' $Name.$Species2.filtered.counts.csv | \
				awk -v OFS=',' 'NR>=2 {if($5 >= '${ReadsPerBarcode[$index]}') print $1,$2,$3,$4} '  > $Name.$Species2.barcodes.txt
				grep -wFf $Name.$Species2.barcodes.txt <(zcat $Name.$Species2.bed.gz) | pigz --fast -p $cores > $Name.$Species2.cutoff.bed.gz
			fi
			# deduplicats for non-UMI reads
		    done
		fi
		
		# Gene body coverage and reads distribution
		if [ ${Type[$index]} == "RNA" ]; then
		    # split bam to hg and mm
		    cd $dir/fastqs/
		    for Species2 in ${Genome2[@]}; do
			#			if [ -f $dir/fastqs/$Name.$Species2.geneBodyCoverage.curves.pdf ]; then
			#			    echo "Skip calculate gene coverage"
			#			else
			#			    echo "Calculate genebody converage and reads distribution of" $Name.$Species2.wdup.bam
			#			    if [ $Runtype = QC ]; then
			#				geneBody_coverage.py -i $dir/fastqs/$Name.$Species2.st.bam -r <(head -n 5000 $myPATH/$Species2.UCSC_RefSeq.bed) -o $Name.$Species2 2>>$dir/Run.log
			#			    else
			#				samtools view -s 0.01 -o $dir/fastqs/$Name.$Species2.temp.bam $dir/fastqs/$Name.$Species2.wdup.bam
			#				geneBody_coverage.py -i $dir/fastqs/$Name.$Species2.temp.bam -r <(head -n 5000 $myPATH/$Species2.UCSC_RefSeq.bed) -o $Name.$Species2 2>>$dir/Run.log
			#				rm $dir/fastqs/$Name.$Species2.temp.bam
			#			    fi
			#			fi
			if [ -f $dir/fastqs/$Name.$Species2.read_distribution.txt ]; then
			    echo "Skip calculate read disbution"
			else
			    echo "Calculate read distribution"
			    if [ $Runtype = QC ]; then
				read_distribution.py -i $dir/fastqs/$Name.$Species2.wdup.bam -r $myPATH/$Species2.UCSC_RefSeq.bed > $Name.$Species2.read_distribution.txt 2>>$dir/Run.log
			    else
				# only use 1% of reads
				samtools view -s 0.01 -o $dir/fastqs/temp.bam $dir/fastqs/$Name.$Species2.wdup.bam
				read_distribution.py -i $dir/fastqs/temp.bam -r $myPATH/$Species2.UCSC_RefSeq.bed > $Name.$Species2.read_distribution.txt 2>>$dir/Run.log
				rm $dir/fastqs/temp.bam
			    fi
			    # plot reads disbution
			    tail -n +5 $Name.$Species2.read_distribution.txt | head -n -1 > temp1.txt
			    head -n 3  $Name.$Species2.read_distribution.txt | grep Assigned | sed 's/Total Assigned Tags//g' | sed 's/ //g' > temp2.txt
			    Rscript $myPATH/Read_distribution.R $dir/fastqs/ $Name.$Species2 --save
			    rm temp1.txt temp2.txt
			fi
		    done
		fi 

		if [ ${Type[$index]} == "ATAC" ] || [ ${Type[$index]} == "TAPS" ]; then
		    # get final quality stats
		    echo -e "Chromosome\tLength\tProperPairs\tBadPairs:Raw" > $Name.$Species.stats.log
		    samtools idxstats $Name.$Species.st.bam >> $Name.$Species.stats.log
		    echo -e "Chromosome\tLength\tProperPairs\tBadPairs:Filtered" >> $Name.$Species.stats.log
		    samtools idxstats $Name.$Species.rmdup.bam >> $Name.$Species.stats.log
		    
		    # get insert-sizes
		    if [ -f $Name.$Species.rmdup.hist_data.pdf ]; then
			echo "Skip checking insert size"
		    else
			echo '' > $Name.$Species.rmdup.hist_data.log
			if [ ${Type[$index]} == "TAPS" ]; then
			    java -jar $picardPATH CollectInsertSizeMetrics VALIDATION_STRINGENCY=SILENT I=$Name.$Species.rmdup.bam O=$Name.$Species.rmdup.hist_data.log H=$Name.$Species.rmdup.hist_data.pdf W=1000  2>>$dir/Run.log
			else
			    java -jar $picardPATH CollectInsertSizeMetrics VALIDATION_STRINGENCY=SILENT I=$Name.$Species.st.bam O=$Name.$Species.rmdup.hist_data.log H=$Name.$Species.rmdup.hist_data.pdf W=1000  2>>$dir/Run.log
			fi
		    fi
		    if [ ! -f $Name.$Species.RefSeqTSS ]; then		
			# make TSS pileup fig
			echo "Create TSS pileup"; set +e
			$toolPATH/pyMakeVplot.py -a $Name.$Species.rmdup.bam -b $tssFilesPATH/$Species.TSS.bed -e 2000 -p ends -v -u -o $Name.$Species.RefSeqTSS
		    fi
		fi
						
       		if [ -d "$dir/fastqs/tmp" ]; then rm -r $dir/fastqs/tmp; fi
	    fi
	done
	
	# count reads for TAPS
	for Species2 in ${Genome2[@]}; do
	    # echo "remove shared reads"
	    #  if [ -f $Name.$Species2.rmdup.cutoff.rmduper.bam.bai ]; then
	    #	echo "found $Name.$Species2.rmdup.cutoff.rmduper.bam.bai, skip filter shared reads"
	    #  else
	    #	$toolPATH/rmduper2.py --bam $Name.$Species2.rmdup.cutoff.bam --out $Name.$Species2.rmdup.cutoff.rmduper.bam
	    #  fi
	    #   get promoter bed file
	    #  if [ -f $Name.$Species2.rmdup.cutoff.2000.500.bed.gz ]; then
	    #	echo "Skip extract promoter reads from rmdup.cutoff.bed.gz"
	    #  else
	    #	echo "Extract promoter reads from rmdup.cutoff.bed.gz"
	    #	bedtools intersect -u -a <(zcat $Name.$Species2.rmdup.cutoff.bed.gz) -b $tssFilesPATH/$Species2.TSS.2000.500bp.bed | \
		#	gzip > $Name.$Species2.rmdup.cutoff.2000.500.bed.gz
	    #  fi
		
	    if [ ${Type[$index]} == "RNA" ] || [ ${Type[$index]} == "crop" ] || [ ${Type[$index]} == "cite" ] || [ ${Type[$index]} == "cellhash" ]; then
	  	# plot UMI/cell or gene/cell
		echo "Convert bed to UMI count matrix, plot UMI_gene_perCell, generate h5"
		/usr/local/bin/Rscript $myPATH/UMI_gene_perCell_plot_v2.R $dir/fastqs/ $Name.$Species2 --save
	    fi
	done

	# count reads for dipC
        for Species2 in ${Genome2[@]}; do
            if [ ${Type[$index]} == "DipC" ]; then
		# count filtered reads
		if [ ! -f $Name.$Species2.filtered.counts.csv ]; then
                    echo "Count filtered reads" $Name.$Species2.rg.cutoff.st.allPairs
		    if [ ! -f $Name.$Species2.rmdup.RG.bed ]; then
			cat $Name.$Species2.rg.cutoff.st.allPairs | cut -f3 > $Name.$Species2.rmdup.RG.bed
                    fi
                    if [ ! -f $Name.$Species2.rmdup.RG.freq.bed ]; then
			cat $Name.$Species2.rmdup.RG.bed | sort --parallel=$cores -S 24G | uniq -c > $Name.$Species2.rmdup.RG.freq.bed
                    fi
                    Rscript $myPATH/sum_reads.R $dir/fastqs/ $Name.$Species2.rmdup.RG.freq.bed --save
                    mv $Name.$Species2.rmdup.RG.freq.bed.csv $Name.$Species2.filtered.counts.csv
                    rm $Name.$Species2.rmdup.RG.bed $Name.$Species.rmdup.RG.freq.bed
		fi
		# remove barcode combinations with less then N reads
		if [ -f $Name.$Species2.final.allValidPairs ]; then
                    echo "Skip removing low counts barcode combination"
		else
                    echo "Remove low counts barcode combination"
		    sed -e 's/,/\t/g' $Name.$Species2.filtered.counts.csv | \
			awk -v OFS=',' 'NR>=2 {if($5 >= '${ReadsPerBarcode[$index]}') print $1,$2,$3,$4} '  > $Name.$Species2.barcodes.txt
		    grep -wFf $Name.$Species2.barcodes.txt  $Name.$Species2.rg.cutoff.st.allPairs > $Name.$Species2.final.allPairs
		    grep -wFf $Name.$Species2.barcodes.txt  $Name.$Species2.rg.cutoff.st.allValidPairs > $Name.$Species2.final.allValidPairs
		fi
	    fi
	done

	# methylation call for TAPS
	for Species2 in ${Genome2[@]}; do
            if [ ${Type[$index]} == "TAPS" ]; then
		if [ -f $Name.$Species2.CpG.bed.gz ]; then
                    echo "Skip removing low counts barcode combination"
		else
		    samtools view -h $Name.$Species2.rmdup.cutoff.namesort.bam | sed 's/\.//g' | samtools view -bS > $Name.$Species2.temp.bam

		    ## mask methylation call in 9bp reads on the left and 1bp on the right
		    cat <(samtools view -H $Name.$Species2.temp.bam) <(samtools view $Name.$Species2.temp.bam |  \
                        awk -v OFS='\t' '{if ($9>0) {gsub("Z:[ZzHhXx.][ZzHhXx.][ZzHhXx.][ZzHhXx.][ZzHhXx.][ZzHhXx.][ZzHhXx.][ZzHhXx.][ZzHhXx.]","Z:.........", $14); print} \
                  		    else {gsub("[ZzHhXx.][ZzHhXx.][ZzHhXx.][ZzHhXx.][ZzHhXx.][ZzHhXx.][ZzHhXx.][ZzHhXx.][ZzHhXx.]$",".........", $14);print}}') | \
			samtools view -bS > $Name.$Species2.bismark.bam
		    rm $Name.$Species2.temp.bam
		    if [ ${GpC[$index]} == "T" ]; then
			(bismark_methylation_extractor --multicore 10 -p $Name.$Species2.bismark.bam --mbias_off --no_header --no_overlap \
						       --gzip --bedGraph --CX) 1>> $Name.$Species2.bismark.log 2>> $Name.$Species2.bismark.log
			## GC and CG report
			(coverage2cytosine --nome-seq --genome_folder ~/Script/Split-seq_Sai/refGenome/fasta/$Species2/ -o $Name.$Species2 \
					   --gzip --split_by_chromosome --gc $Name.$Species2.bismark.bismark.cov.gz) 1>> $Name.$Species2.bismark.log 2>> $Name.$Species2.bismark.log

			## calculate mean methylation% and GpC MTase effciency
			zcat $Name.$Species2.NOMe.CpG_report.txt.gz | \
			    awk -v OFS='\t' '{unme+=$4}{me+=$5}END{print "CpG", me/(unme+me)}' > $Name.$Species2.methy.pct.xls
			zcat $Name.$Species2.NOMe.GpC_report.txt.gz | \
                            awk -v OFS='\t' '{unme+=$4}{me+=$5}END{print "GpC", me/(unme+me)}' >> $Name.$Species2.methy.pct.xls

		    else
			(bismark_methylation_extractor --multicore 10 -p $Name.$Species2.bismark.bam --mbias_off --no_header --gzip --bedGraph --no_overlap) 1>> $Name.$Species2.bismark.log 2>> $Name.$Species2.bismark.log
			zcat $Name.$Species2.bismark.bismark.cov.gz | \
                            awk -v OFS='\t' '{unme+=$5}{me+=$6}END{print "CpG", me/(unme+me)}' > $Name.$Species2.methy.pct.xls
		    fi
		    zcat $Name.$Species2.bismark.bedGraph.gz | awk -v OFS='\t' 'NR > 1 {print $1, $2, $3, 100-$4}' | pigz --fast -p $cores > $Name.$Species2.corrected.bedGraph.gz
		    
		    ## convert to fragments file and merge top and bottom
		    ## convert Z to "unmethylation", 0
		    zcat CpG_OT_$Name.$Species2.bismark.txt.gz CpG_CTOT_$Name.$Species2.bismark.txt.gz| \
			sed 's/_/\t/g' | awk -v OFS='\t' '{if($6 == "z") print $4, $5, $5+1, $2, 1, "+"; else print $4, $5, $5+1, $2, 0, "+"}' | \
			sed 's/R1/R1./g' | sed 's/R2/R2./g' | sed 's/R3/R3./g' | sed 's/P1/P1./g' | \
			pigz --fast -p $cores > $Name.$Species2.CpG_OT.bed.gz
		    zcat CpG_OB_$Name.$Species2.bismark.txt.gz CpG_CTOB_$Name.$Species2.bismark.txt.gz| \
		        sed 's/_/\t/g' | awk -v OFS='\t' '{if($6 == "z") print $4, $5, $5+1, $2, 1, "-"; else print $4, $5, $5+1, $2, 0, "-"}' | \
			sed 's/R1/R1./g' | sed 's/R2/R2./g' | sed 's/R3/R3./g' | sed 's/P1/P1./g' | \
			pigz --fast -p $cores > $Name.$Species2.CpG_OB.bed.gz
		    cat $Name.$Species2.CpG_OT.bed.gz $Name.$Species2.CpG_OB.bed.gz > $Name.$Species2.CpG.bed.gz
		    rm $Name.$Species2.CpG_OT.bed.gz $Name.$Species2.CpG_OB.bed.gz
		fi
	    fi
	done

	# estimate lib size
	if [ -f $Name.counts.csv ]; then
	    echo "Found $Name.counts.csv, skip calculate lib size"
	else
	    echo "Estimate lib size"
            if [ ${Genomes[$index]} == "both" ]; then
		if [ -f $Name.hg19.unfiltered.counts.csv ] && [ ! -f $Name.hg19.filtered.counts.csv ]; then
		    cp $Name.hg19.unfiltered.counts.csv $Name.hg19.filtered.counts.csv
		    echo "Error: Could locate $Name.hg19.unfiltered.counts.csv"
		fi
		if [ -f $Name.mm10.unfiltered.counts.csv ] && [ ! -f $Name.mm10.filtered.counts.csv ]; then
		    cp $Name.mm10.unfiltered.counts.csv $Name.mm10.filtered.counts.csv
		    echo "Error: Could locate $Name.mm10.unfiltered.counts.csv"
		fi
		if [ -f $Name.hg19.unfiltered.counts.csv ] && [ -f $Name.mm10.unfiltered.counts.csv ]; then
                    echo "Calcuating library size for $Name"
		    Rscript $myPATH/lib_size_sc_V5_species_mixing.R ./ $Name ${ReadsPerBarcode[$index]} ${Type[$index]} --save
		fi
	    else
		if [ -f $Name.${Genomes[$index]}.unfiltered.counts.csv ] && [ ! -f $Name.${Genomes[$index]}.filtered.counts.csv ]; then
		    cp $Name.${Genomes[$index]}.unfiltered.counts.csv $Name.${Genomes[$index]}.filtered.counts.csv
		else
		    Rscript $myPATH/lib_size_sc_V5_single_species.R ./ $Name ${ReadsPerBarcode[$index]} ${Genomes[$index]} ${Type[$index]} --save
		fi
	    fi
	fi
	if [ ! -d $dir/$Name/ ]; then 
	    mkdir $dir/$Name/ && mv $dir/fastqs/$Name.* $dir/$Name
	    if [ ${Type[$index]} == "TAPS" ]; then
		mv $dir/fastqs/*$Name.*txt.gz $dir/$Name
	    fi
	fi
    fi
    index=$((index + 1))
done
