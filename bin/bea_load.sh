#!/usr/bin/env bash

#./../bea/bin/cea_load.sh http://0:3030 1613 /home/vagrant/cxgn /home/vagrant/cxgn/cassbase/bin true true localhost cea postgres postgres cass_index_CASS_6Genotypes_Sampling_2015
#./../bea/bin/cea_load.sh https://cassbase.org 1611 /home/vagrant/cxgn /home/vagrant/cxgn/cassbase/bin false true localhost cea postgres postgres cass_index_T100

host_address=$1 #http://0:3000
trial_ids=$2 #[1613,1612]
accession_ids=$3
trait_ids=$4
root_dir=$5 #/home/vagrant/cxgn
correlation_index_dir=$6 #/home/vagrant/cxgn/indexed_files/correlation_indexes
expression_index_dir=$7 #/home/vagrant/cxgn/indexed_files/expression_indexes
description_index_dir=$8 #/home/vagrant/cxgn/indexed_files/loci_and_description_index
processing_file_dir=$9 #/home/vagrant/cxgn/cassbase/bin
clear_old_database=${10} #true or false
clear_old_index_files=${11} #true or false
database_host=${12}
database_name=${13}
database_user=${14}
database_password=${15}
directory_name=${16} #cass_index_<Project_name>
export_type=${17} #1 or 2
project_name=${18} #<Project_name>
sgn_session_id=${19}
user_name=${20}
search_params_description=${21}

wget --no-check-certificate --output-document=${processing_file_dir}/phenotype_download.csv "${host_address}/breeders/trials/phenotype/download?trial_list=${trial_ids}&format=csv&timestamp=0&trait_list=${trait_ids}&accession_list=${accession_ids}&dataLevel=plot&search_type=complete&sgn_session_id=${sgn_session_id}&has_header=0"

perl ${root_dir}/bea/bin/transform_phenotypes_to_expression_atlas_lucy_index.pl -i ${processing_file_dir}/phenotype_download.csv -o ${processing_file_dir}/lucy.tsv -c ${processing_file_dir}/pre_corr.tsv -f ${processing_file_dir}/corr.tsv -p ${processing_file_dir}/project.txt -d ${processing_file_dir}/desc.tsv -v ${export_type} -n ${project_name} -t ${processing_file_dir} -u ${user_name} -s ${search_params_description}

if ${clear_old_database}; then
PGPASSWORD=${database_password} dropdb -U ${database_user} -h ${database_host} ${database_name}
PGPASSWORD=${database_password} createdb -U ${database_user} -h ${database_host} ${database_name}
#DROP SCHEMA public CASCADE;
#CREATE SCHEMA public;
PGPASSWORD=${database_password} psql -U ${database_user} -d ${database_name} -h ${database_host} -a -f ${root_dir}/Tea/import_project/create_tea_schema.sql
fi

mkdir ${correlation_index_dir}/${directory_name}
mkdir ${expression_index_dir}/${directory_name}
mkdir ${description_index_dir}/${directory_name}

if ${clear_old_index_files}; then
rm -R ${correlation_index_dir}/${directory_name}/*
rm -R ${expression_index_dir}/${directory_name}/*
rm -R ${description_index_dir}/${directory_name}/*
fi

perl ${root_dir}/Tea/import_project/TEA_import_project_metadata.pl -d ${database_name} -H ${database_host} -u ${database_user} -p ${database_password} -t ${processing_file_dir}/project.txt -n

perl ${root_dir}/Tea/import_project/index_correlation_file.pl ${processing_file_dir}/corr.tsv ${correlation_index_dir}/${directory_name}/
perl ${root_dir}/Tea/import_project/index_expression_file.pl ${processing_file_dir}/lucy.tsv ${expression_index_dir}/${directory_name}/
perl ${root_dir}/Tea/import_project/index_description_file.pl ${processing_file_dir}/desc.tsv ${description_index_dir}/${directory_name}/

