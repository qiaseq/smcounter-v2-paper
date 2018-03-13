#!/bin/bash

# output directory with smCounter results
echo $1,$2
runPath=$1
result=$runPath/result
vcf=$runPath/vcf
intermediate=$runPath/intermediate
processed=$runPath/processed
tmp=$runPath/tmp
misc=$runPath/misc
final=$runPath/final

# output prefix
caller=$2
caller2=$3

# common files
filePath=/mnt/webserver/datadisk/resources/varcall/frequentlyUsedFiles
RefSDF=$filePath/ucsc.hg19.sdf


# create directory to keep processed vcfs
if [ ! -s $processed ]; then
   mkdir $processed
fi
# create results directory
if [ ! -s $result ]; then
   mkdir $result
fi
# create tmp directory
if [ ! -s $tmp ]; then
   mkdir $tmp
fi
# create final directory
if [ ! -s $final ]; then
   mkdir $final
fi


for region in all coding noncoding; do
  subcaller=$caller.$region
  # create results cutoff directory
  if [ ! -s $result/$subcaller ]; then
      mkdir $result/$subcaller
  fi

  if [ $region = all ]; then
    GroundTruth=$misc/na12878.uniq.all.het.vcf.gz
    bkgVcf=$misc/na24385.plus.na12878.hom.all.vcf
    targetRegionSize=850316
    finalBed=$misc/$caller2.hc.$region.na12878.na24385.v3.3.2.bed
    noheaderVcf=$misc/na12878.uniq.$region.het.noheader.vcf
    overlap=1.0
  fi

  if [ $region = coding ]; then
    GroundTruth=$misc/na12878.uniq.coding.het.vcf.gz
    bkgVcf=$misc/na24385.plus.na12878.hom.coding.vcf
    targetRegionSize=591154
    finalBed=$misc/$caller2.hc.$region.na12878.na24385.v3.3.2.bed
    noheaderVcf=$misc/na12878.uniq.$region.het.noheader.vcf
    overlap=1.0
  fi

  if [ $region = noncoding ]; then
    GroundTruth=$misc/na12878.uniq.noncoding.het.vcf.gz
    bkgVcf=$misc/na24385.plus.na12878.hom.noncoding.vcf
    targetRegionSize=259162
    finalBed=$misc/$caller2.hc.$region.na12878.na24385.v3.3.2.bed
    noheaderVcf=$misc/na12878.uniq.$region.het.noheader.vcf
    overlap=1E-9
  fi

  # remove bkg variants and out-of-HC variants
  grep -v LM $intermediate/$caller.VariantList.long.txt > tmp.long.txt
  python /mnt/webserver/datadisk/resources/varcall/smCounter.v2.paper_backup/code/HCandDilute.fast.alt.py  $runPath tmp.long.txt  $intermediate/$subcaller.VariantList.long.hc.dil.txt  $finalBed  $overlap $bkgVcf
  rm tmp.long.txt
  
  for cutoff in 1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5 5.0 5.5 6.0 6.5 7.0 7.5 8.0 8.5 9.0 9.5 10.0; do
    # create results cutoff directory. if already exists, delete and re-create
    if [ ! -s $result/$subcaller/$cutoff ]; then
      mkdir $result/$subcaller/$cutoff
    else
      rm -r $result/$subcaller/$cutoff
      mkdir $result/$subcaller/$cutoff
    fi

    python  /mnt/webserver/datadisk/resources/varcall/smCounter.v2.paper_backup/code/HCandDilute.fast.py $runPath  $vcf/$caller.$cutoff.vcf $vcf/$subcaller.$cutoff.hc.dil.vcf  $finalBed  $overlap   $bkgVcf
    cat $filePath/VCFheader_dup.txt $vcf/$subcaller.$cutoff.hc.dil.vcf > $vcf/$subcaller.$cutoff.hc.dil.header.vcf

    file=$vcf/$subcaller.$cutoff.hc.dil.header.vcf
    awk '{if (/^#/) {print} else {if (/^chr/) {print} else {print "chr"$0}}}' $file > $processed/$subcaller.$cutoff.hc.dil.header.vcf
    bedtools sort -header -i $processed/$subcaller.$cutoff.hc.dil.header.vcf | bgzip -c > $processed/$subcaller.$cutoff.hc.dil.header.vcf.gz
    tabix -f -p vcf $processed/$subcaller.$cutoff.hc.dil.header.vcf.gz > /dev/null
    
    # Run RTG evaluation
    file=$processed/$subcaller.$cutoff.hc.dil.header.vcf.gz
    # check if RTG results not already exist
    if [ -s $result/$subcaller/$cutoff/filter.applied ]; then
      echo "RTG tools evaluation results already exist for ${file:10}"
    else
      cmd="/mnt/webserver/datadisk/resources/varcall/software/RTG/rtg-tools-3.5/rtg vcfeval --all-records -b $GroundTruth -c $file --squash-ploidy --no-gzip -t $RefSDF -f \"QUAL\" -o $result/$subcaller/$cutoff/filter.ignored &> $runPath/rtg.log"
      echo $cmd
      eval $cmd
      cmd="/mnt/webserver/datadisk/resources/varcall/software/RTG/rtg-tools-3.5/rtg vcfeval --baseline-tp -b $GroundTruth -c $file --squash-ploidy --no-gzip -t $RefSDF -f \"QUAL\" -o $result/$subcaller/$cutoff/filter.applied &>> $runPath/rtg.log"
      echo $cmd
      eval $cmd
    fi

    ## create tables by merging info from TP, FP, and FN vcf files
    awk '!/^#/ {OFS="\t"; if ($7 == "PASS") {Pred="TP"} else {Pred="FN"}; print $1":"$2";"$4"/"$5,Pred}' $result/$subcaller/$cutoff/filter.ignored/tp.vcf > $tmp/tp.txt 2>>$runPath/rtg.log
    awk '!/^#/ {OFS="\t"; if ($7 == "PASS") {Pred="FP"} else {Pred="TN"}; print $1":"$2";"$4"/"$5,Pred}' $result/$subcaller/$cutoff/filter.ignored/fp.vcf > $tmp/fp.txt 2>>$runPath/rtg.log
    awk '!/^#/ {OFS="\t"; print $1":"$2";"$4"/"$5,"FN"}' $result/$subcaller/$cutoff/filter.ignored/fn.vcf > $tmp/fn.txt 2>>$runPath/rtg.log
    awk '!/^#/ {OFS="\t"; print $1":"$2";"$4"/"$5,"TP"}' $result/$subcaller/$cutoff/filter.applied/tp-baseline.vcf > $result/$subcaller/$cutoff/tp-baseline.txt 2>>$runPath/rtg.log

    # merge tables to one for one cutoff and remove the duplicates (caused by RTG tools limitations/bugs)
    cat $tmp/tp.txt $tmp/fp.txt $tmp/fn.txt | sort | uniq | awk -v duplicate=0 '{OFS="\t"; indb=ind; typeb=type; ind=$1; type=$2; if (ind==indb) {if (typeb ~ /T/) {print indb,typeb} else if (type ~ /T/) {print ind,type}; duplicate=1} else if (duplicate==0) {print indb,typeb} else {duplicate=0}}' | awk 'NF>0' > $result/$subcaller/$cutoff/table.txt 2>>$runPath/rtg.log
  done    # end of cutoff loop
#
  # integrate tables to one Table for each group -- this is at pct level, after looping all cutoffs
  if [ ! -s $result/$subcaller ]; then
    echo "warning: no RTG tools results found for group: " $pct
  else
  # initialize
    count=1
    OutField="0"
    OutHeader="#Variant"

    # remove Table.txt and Table-baseline.txt from previous runs
    if [ -f $result/$subcaller/Table.txt ]; then
      rm $result/$subcaller/Table.txt
    fi
    if [ -f $result/$subcaller/Table-baseline.txt ]; then
      rm $result/$subcaller/Table-baseline.txt
    fi

    # folders for each cutoff
    cutOffs=$(ls $result/$subcaller)
    for cutOff in $cutOffs; do
      if [ $count -eq 1 ]; then
        cp $result/$subcaller/${cutOff}/table.txt $result/$subcaller/Table.txt
        cp $result/$subcaller/${cutOff}/tp-baseline.txt $result/$subcaller/Table-baseline.txt
      else
        OutField=${OutField}",1."$count
        # join tables from all cutoff values
        join -a1 -a2 -1 1 -2 1 -e "TN" -o ${OutField}",2.2" -t $'\t' <(sort -t $'\t' -k1,1 $result/$subcaller/Table.txt | uniq) <(sort -t $'\t' -k1,1 $result/$subcaller/${cutOff}/table.txt | uniq) > $result/$subcaller/Table.tmp
        mv $result/$subcaller/Table.tmp $result/$subcaller/Table.txt
        # same for baseline
        join -a1 -a2 -1 1 -2 1 -e "FN" -o ${OutField}",2.2" -t $'\t' <(sort -t $'\t' -k1,1 $result/$subcaller/Table-baseline.txt | uniq) <(sort -t $'\t' -k1,1 $result/$subcaller/${cutOff}/tp-baseline.txt | uniq) > $result/$subcaller/Table.tmp
        mv $result/$subcaller/Table.tmp $result/$subcaller/Table-baseline.txt
      fi

      count=$((count+1))
      OutHeader=${OutHeader}$'\t'"Thr="$cutOff
    done # end of cutoff folder loop

    echo "$OutHeader" > $result/$subcaller/Header.txt
  fi

  # add header info (cutoff values)
  sort $result/$subcaller/Table.txt | uniq > $result/$subcaller/Table.tmp
  cat $result/$subcaller/Header.txt $result/$subcaller/Table.tmp > $result/$subcaller/Table.txt
  sort $result/$subcaller/Table-baseline.txt | uniq > $result/$subcaller/Table-baseline.tmp
  cat $result/$subcaller/Header.txt $result/$subcaller/Table-baseline.tmp > $result/$subcaller/Table-baseline.txt
  rm $result/$subcaller/Header.txt $result/$subcaller/Table.tmp $result/$subcaller/Table-baseline.tmp

  # get summary results
  Rscript /mnt/webserver/datadisk/resources/varcall/smCounter.v2.paper_backup/code/validateWithGIB.withRTG.v2.2.R $runPath  $intermediate/$subcaller.VariantList.long.hc.dil.txt $noheaderVcf $result/$subcaller/Table.txt $result/$subcaller/Table-baseline.txt $targetRegionSize $final/roc.$subcaller.png $final/summary.$subcaller.csv $final/details.$subcaller.csv $final/tpfpfn.$subcaller.csv $subcaller 1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5 5.0 5.5 6.0 6.5 7.0 7.5 8.0 8.5 9.0 9.5 10.0 >/dev/null 2>$runPath/validateR.log

done
  
