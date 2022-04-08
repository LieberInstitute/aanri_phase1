#!/bin/bash

SOFTWARE="/ceph/opt/PGDSpider_2.1.1.5"
TISSUE="caudate"

echo $TISSUE
mkdir $TISSUE
ls ../../_m/$TISSUE/*ped | \
    parallel "java -Xmx1G -Xms512m -jar $SOFTWARE/PGDSpider2-cli.jar \
                   -inputfile {} -inputformat PED \
                   -outputfile $TISSUE/{/.}.txt -outputformat GENEPOP \
                   -spid ../_h/PED_to_GENEPOP.spid"
