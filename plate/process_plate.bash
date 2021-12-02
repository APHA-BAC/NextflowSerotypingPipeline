

# TODO: Download from s3


plate_name=new-plate
reads=/home/aaronfishman/temp/salmonella-reads/$plate_name/
results=/home/aaronfishman/temp/salmonella-results/$plate_name/

sudo docker run -it \
    -v $reads:/root/WGS_Data/$plate_name/ \
    -v $results:/root/WGS_Results/$plate_name/ \
    jguzinski/salmonella-seq:prod \
    /root/nextflow/nextflow SCE3_pipeline_update.nf \
    --runID $plate_name
