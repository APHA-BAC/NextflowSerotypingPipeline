#!/usr/bin/env nextflow

/*
 * STEP 0 - define the input path to the sequences that will be analysed 
*/ 

params.reads = "$projectDir/test_isolates/*_{R1,R2}.fastq.gz"

/*
 * STEP 1 - fastqc
*/ 
reads = Channel.fromFilePairs(params.reads)

process fastqc {
    publishDir "/home/WGS_Results/TestIsolates/${sample_id}/FASTQC_Reports", mode: 'move'

    input:
    tuple sample_id, file(reads_file) from reads   

    output:  
    file("${sample_id}_{R1,R2}_fastqc.{html,zip}")
    file("${sample_id}_1.txt") into out1_ch
    file("${sample_id}_1.txt") into out1_ch_rem

    
    script:
    """
    fastqc  -f fastq -q ${reads_file}
   > ${sample_id}_1.txt 
    """
}


/*
 * STEP 2 - shovill
*/ 
reads1 = Channel.fromFilePairs(params.reads)

process shovill {
    publishDir "/home/WGS_Results/TestIsolates/${sample_id}/shovill", mode: 'copy'
    
    input:
    tuple sample_id, file(reads_file) from reads1
  
    output:
    file "*" into shovill_all_ch
    set sample_id, file("${sample_id}_contigs.fa") into quast_ch  
    set sample_id, file("${sample_id}_contigs.fa") into sistr_ch  
    file("${sample_id}_2.txt") into out2_ch
    file("${sample_id}_2.txt") into out2_ch_rem

  
    script:
    """    
    /opt/conda/bin/shovill --ram 7 --R1 ${sample_id}_R1.fastq.gz --R2 ${sample_id}_R2.fastq.gz
    mv contigs.fa ${sample_id}_contigs.fa
   > ${sample_id}_2.txt 
   """
}


/*
 * STEP 3 - quast 
*/ 

reads2 = Channel.fromFilePairs(params.reads)
quast_ch
.join(reads2)
.set { quast_in }

process quast {
    publishDir "/home/WGS_Results/TestIsolates/${sample_id}/quast",  mode: 'copy'
    
    input:
    set sample_id, file("${sample_id}_contigs.fa"), file (reads_file) from quast_in

    output:
    file("${sample_id}_3.txt") into out3_ch
    file("${sample_id}_3.txt") into out3_ch_rem

   
    script:
    """
    python /NextflowSerotypingPipeline/quast-5.1.0rc1/quast.py -o /home/WGS_Results/TestIsolates/${sample_id}/quast "${sample_id}_contigs.fa"
    > ${sample_id}_3.txt
    """
}


/*
 * STEP 4 - kmerid 
*/ 
reads3 = Channel.fromFilePairs(params.reads)

process kmerid {
    publishDir "/home/WGS_Results/TestIsolates/${sample_id}/Kmerid",  mode: 'copy'

    input:
    tuple sample_id, file(reads_file) from reads3

    output:  
    file("${sample_id}_R1.tsv")
    file("${sample_id}_4.txt") into out4_ch
    file("${sample_id}_4.txt") into out4_ch_rem
 
    script:
    """     
    python /opt/kmerid/kmerid_python3.py -f /home/WGS_Data/TestIsolates/${sample_id}_R1.fastq.gz  -c /opt/kmerid/config/config.cnf -n > ${sample_id}_R1.tsv
   > ${sample_id}_4.txt 
    """
}


/*
 * STEP 5 - seqsero2 
*/ 
reads4 = Channel.fromFilePairs(params.reads)

process seqsero2 {
   publishDir "/home/WGS_Results/TestIsolates/${sample_id}/SeqSero2", mode: 'copy'


    input:
    tuple sample_id, file(reads_file) from reads4

    output:
    file("${sample_id}_5.txt") into out5_ch   
    file("${sample_id}_5.txt") into out5_ch_rem   
   
    script:
    """     
   #/opt/conda/bin/conda init bash
   /opt/conda/bin/SeqSero2_package.py -m a -b mem -t 2 -i \$PWD/${sample_id}_{R1,R2}.fastq.gz -d /home/WGS_Results/TestIsolates/${sample_id}/SeqSero2 > ${sample_id}_5.txt
    
    """
}


/*
 * STEP 6 - sistr 
*/ 
reads5 = Channel.fromFilePairs(params.reads)
sistr_ch
.join(reads5)
.set { sistr_in }


process sistr {
   publishDir "/home/WGS_Results/TestIsolates/${sample_id}/sistr", mode: 'move'


    input:
    set sample_id, file("${sample_id}_contigs.fa"), file (reads_file) from sistr_in
    
    output:
    file 'sistr_prediction.csv' into sistr_out_ch  
    file("${sample_id}_6.txt") into out6_ch   
    file("${sample_id}_6.txt") into out6_ch_rem   

       
    script:
    """     
    /opt/conda/bin/sistr -i "${sample_id}_contigs.fa" ${sample_id} -f csv -o sistr_prediction.csv --qc
    > ${sample_id}_6.txt
    """
}


/*
 * STEP 7 - MOST 
*/ 
reads6 = Channel.fromFilePairs(params.reads)

process most {
   publishDir "/home/WGS_Results/TestIsolates/${sample_id}/MOST", mode: 'copy'


    input:
    tuple sample_id, file(reads_file) from reads6

    output:
    file("${sample_id}_7.txt") into out7_ch   
    file("${sample_id}_7.txt") into out7_ch_rem   
    set sample_id, file("${sample_id}_serovar.tsv") into most_out_ch

   
    script:
    """
    MOST_DIR=\$PWD/MOST     
    python /opt/most/MOST-master/MOST.py -1 \$PWD/${sample_id}_R1.fastq.gz  -2 \$PWD/${sample_id}_R2.fastq.gz -st /opt/most/MOST-master/MLST_data/salmonella --output_directory \$MOST_DIR -serotype True --bowtie /opt/most/bowtie2-2.1.0/bowtie2 --samtools /opt/most/samtools-0.1.18/samtools
    if grep "predicted_serotype" \$MOST_DIR/MOST/${sample_id}_R1.fastq.results.xml
    then
    grep "predicted_serotype" \$MOST_DIR/${sample_id}_R1.fastq.results.xml >> serovar1.txt
    if grep -q "ST-serotype" serovar1.txt
    then
    awk '{print substr(\$2,1,5); }' serovar1.txt > serovar2.txt
    mv serovar2.txt  ${sample_id}_serovar.tsv 
    else
    awk '{print substr(\$3,10); }' serovar1.txt > serovar2.txt   
    mv serovar2.txt  ${sample_id}_serovar.tsv 
    fi
    else
    grep "profile" \$MOST_DIR/${sample_id}_R1.fastq.results.xml >> serovar1.txt
    awk '{print substr(\$3,1,5); }' serovar1.txt > serovar2.txt
    mv serovar2.txt  ${sample_id}_serovar.tsv 
    fi
    > ${sample_id}_7.txt
    """
}


/*
 * STEP 8 - srst2
*/ 

reads7 = Channel.fromFilePairs(params.reads)
most_out_ch
.join(reads7)
.set { sero_in }


process srst2 {
   publishDir "/home/WGS_Results/TestIsolates/${sample_id}/srst2", mode: 'copy'


    input:
    set sample_id, file("${sample_id}_serovar.tsv"), file (reads_file) from sero_in



    output:
    file("${sample_id}_8.txt") into out8_ch   
    file("${sample_id}_8.txt") into out8_ch_rem   


    script:    

          """
    if grep -E "(Typhimurium|Enteritidis|Gallinarum|Pullorum|Idikan|Kedougou|Java|Paratyphi)" ${sample_id}_serovar.tsv 
    then
    export SRST2_BOWTIE2=/opt/srst2/bowtie2-2.2.3/bowtie2
    export SRST2_BOWTIE2_BUILD=/opt/srst2/bowtie2-2.2.3/bowtie2-build
    srst2.py  --input_pe ${sample_id}_R1.fastq.gz ${sample_id}_R2.fastq.gz --forward _R1 --reverse _R2 --output /home/WGS_Results/TestIsolates/${sample_id}/srst2/ --log --gene_db /opt/srst2/VaccineDifferentiation/allVacDB9h1_clustered.fasta
    > ${sample_id}_8.txt    
    else 
    > ${sample_id}_8.txt
    fi
        """
    
}


/*
 * STEP 9 - summary 
*/ 
reads_summ = Channel.fromFilePairs(params.reads)

process summary {

  input:  
  val read from reads_summ.count()
  .view()

  val c1 from out1_ch.count()
  .view()

  val c2 from out2_ch.count()
  .view()

  val c3 from out3_ch.count()
  .view()

  val c4 from out4_ch.count()
  .view()

  val c5 from out5_ch.count()
  .view()

  val c6 from out6_ch.count()
  .view()

  val c7 from out7_ch.count()
  .view()

 val c8 from out8_ch.count()
  .view()


  when: 
  c1 + c2 + c3 + c4+ c5+ c6 + c7 + c8 == read*8
   
  script: 
  """ 
  python /NextflowSerotypingPipeline/summaryTable_Python3.py \$PWD 
  """
}


/*
 * STEP 10 - remove text files 
*/ 

reads8 = Channel.fromFilePairs(params.reads)
reads_summ1 = Channel.fromFilePairs(params.reads)

process remove {

  input:  
  tuple sample_id, file(reads_file) from reads8   

  input:  
  val read_rem from reads_summ1.count()
  .view()

  val c1_rem from out1_ch_rem.count()
  .view()

  val c2_rem from out2_ch_rem.count()
  .view()

  val c3_rem from out3_ch_rem.count()
  .view()

  val c4_rem from out4_ch_rem.count()
  .view()

  val c5_rem from out5_ch_rem.count()
  .view()

  val c6_rem from out6_ch_rem.count()
  .view()

  val c7_rem from out7_ch_rem.count()
  .view()

  val c8_rem from out8_ch_rem.count()
  .view()

  when: 
  c1_rem + c2_rem + c3_rem + c4_rem + c5_rem+ c6_rem + c7_rem + c8_rem == read_rem*8
   
  script: 
  """ 
  rm /home/WGS_Results/TestIsolates/${sample_id}/FASTQC_Reports/${sample_id}_1.txt
  rm /home/WGS_Results/TestIsolates/${sample_id}/shovill/${sample_id}_2.txt 
  rm /home/WGS_Results/TestIsolates/${sample_id}/quast/${sample_id}_3.txt
  rm /home/WGS_Results/TestIsolates/${sample_id}/Kmerid/${sample_id}_4.txt
  rm /home/WGS_Results/TestIsolates/${sample_id}/SeqSero2/${sample_id}_5.txt
  rm /home/WGS_Results/TestIsolates/${sample_id}/sistr/${sample_id}_6.txt
  rm /home/WGS_Results/TestIsolates/${sample_id}/MOST/${sample_id}_7.txt
  rm /home/WGS_Results/TestIsolates/${sample_id}/srst2/${sample_id}_8.txt
  """
}
