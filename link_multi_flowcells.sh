#!/bin/bash
chmod u+w ./data/raw/run/BSF_0666_HCW2YDRXX
if [ -d ./data/raw/run/BSF_0666_HCW2YDRXX/BSF_0666_HCW2YDRXX_3_samples/ ];then
  rm -r ./data/raw/run/BSF_0666_HCW2YDRXX/BSF_0666_HCW2YDRXX_3_samples/
fi
mkdir -p ./data/raw/run/BSF_0666_HCW2YDRXX/BSF_0666_HCW2YDRXX_3_samples/

for bam_file in $(ls ./data/raw/run/BSF_0670_HFCK2DRXX/BSF_0670_HFCK2DRXX_1_samples);do
  new_name=$(echo $bam_file|sed 's/BSF_0670_HFCK2DRXX_1/BSF_0666_HCW2YDRXX_3/')
  ln -s ../../BSF_0670_HFCK2DRXX/BSF_0670_HFCK2DRXX_1_samples/$bam_file ./data/raw/run/BSF_0666_HCW2YDRXX/BSF_0666_HCW2YDRXX_3_samples/$new_name
done

if [ -d ./data/raw/run/BSF_0666_HCW2YDRXX/BSF_0666_HCW2YDRXX_4_samples/ ];then
  rm -r ./data/raw/run/BSF_0666_HCW2YDRXX/BSF_0666_HCW2YDRXX_4_samples/
fi
mkdir -p ./data/raw/run/BSF_0666_HCW2YDRXX/BSF_0666_HCW2YDRXX_4_samples/

for bam_file in $(ls ./data/raw/run/BSF_0670_HFCK2DRXX/BSF_0670_HFCK2DRXX_2_samples);do
  new_name=$(echo $bam_file|sed 's/BSF_0670_HFCK2DRXX_2/BSF_0666_HCW2YDRXX_4/')
  ln -s ../../BSF_0670_HFCK2DRXX/BSF_0670_HFCK2DRXX_2_samples/$bam_file ./data/raw/run/BSF_0666_HCW2YDRXX/BSF_0666_HCW2YDRXX_4_samples/$new_name
done
chmod u-w ./data/raw/run/BSF_0666_HCW2YDRXX

chmod u+w ./data/raw/run/BSF_0674_HL7L5DMXX

if [ -d ./data/raw/run/BSF_0674_HL7L5DMXX/BSF_0674_HL7L5DMXX_3_samples/ ];then
  rm -r ./data/raw/run/BSF_0674_HL7L5DMXX/BSF_0674_HL7L5DMXX_3_samples/
fi
mkdir -p ./data/raw/run/BSF_0674_HL7L5DMXX/BSF_0674_HL7L5DMXX_3_samples/
for bam_file in $(ls ./data/raw/run/BSF_0678_HM7GMDMXX/BSF_0678_HM7GMDMXX_1_samples);do
  new_name=$(echo $bam_file|sed 's/BSF_0678_HM7GMDMXX_1/BSF_0674_HL7L5DMXX_3/')
  ln -s ../../BSF_0678_HM7GMDMXX/BSF_0678_HM7GMDMXX_1_samples/$bam_file ./data/raw/run/BSF_0674_HL7L5DMXX/BSF_0674_HL7L5DMXX_3_samples/$new_name
done

if [ -d ./data/raw/run/BSF_0674_HL7L5DMXX/BSF_0674_HL7L5DMXX_4_samples/ ];then
  rm -r ./data/raw/run/BSF_0674_HL7L5DMXX/BSF_0674_HL7L5DMXX_4_samples/
fi
mkdir -p ./data/raw/run/BSF_0674_HL7L5DMXX/BSF_0674_HL7L5DMXX_4_samples/
for bam_file in $(ls ./data/raw/run/BSF_0678_HM7GMDMXX/BSF_0678_HM7GMDMXX_2_samples);do
  new_name=$(echo $bam_file|sed 's/BSF_0678_HM7GMDMXX_2/BSF_0674_HL7L5DMXX_4/')
  ln -s ../../BSF_0678_HM7GMDMXX/BSF_0678_HM7GMDMXX_2_samples/$bam_file ./data/raw/run/BSF_0674_HL7L5DMXX/BSF_0674_HL7L5DMXX_4_samples/$new_name
done
chmod u-w ./data/raw/run/BSF_0674_HL7L5DMXX

chmod u+w ./data/raw/run/BSF_0684_HGNMYDRXX
if [ -d ./data/raw/run/BSF_0684_HGNMYDRXX/BSF_0684_HGNMYDRXX_3_samples/ ];then
  rm -r ./data/raw/run/BSF_0684_HGNMYDRXX/BSF_0684_HGNMYDRXX_3_samples/
fi
mkdir -p ./data/raw/run/BSF_0684_HGNMYDRXX/BSF_0684_HGNMYDRXX_3_samples/

for bam_file in $(ls ./data/raw/run/BSF_0688_HMCGTDMXX/BSF_0688_HMCGTDMXX_1_samples);do
  new_name=$(echo $bam_file|sed 's/BSF_0688_HMCGTDMXX_1/BSF_0684_HGNMYDRXX_3/')
  ln -s ../../BSF_0688_HMCGTDMXX/BSF_0688_HMCGTDMXX_1_samples/$bam_file ./data/raw/run/BSF_0684_HGNMYDRXX/BSF_0684_HGNMYDRXX_3_samples/$new_name
done

if [ -d ./data/raw/run/BSF_0684_HGNMYDRXX/BSF_0684_HGNMYDRXX_4_samples/ ];then
  rm -r ./data/raw/run/BSF_0684_HGNMYDRXX/BSF_0684_HGNMYDRXX_4_samples/
fi
mkdir -p ./data/raw/run/BSF_0684_HGNMYDRXX/BSF_0684_HGNMYDRXX_4_samples/

for bam_file in $(ls ./data/raw/run/BSF_0688_HMCGTDMXX/BSF_0688_HMCGTDMXX_2_samples);do
  new_name=$(echo $bam_file|sed 's/BSF_0688_HMCGTDMXX_2/BSF_0684_HGNMYDRXX_4/')
  ln -s ../../BSF_0688_HMCGTDMXX/BSF_0688_HMCGTDMXX_2_samples/$bam_file ./data/raw/run/BSF_0684_HGNMYDRXX/BSF_0684_HGNMYDRXX_4_samples/$new_name
done
chmod u-w ./data/raw/run/BSF_0684_HGNMYDRXX
