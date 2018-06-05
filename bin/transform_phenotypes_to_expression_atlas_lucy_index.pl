#!/usr/bin/perl

=head1

transform_phenotypes_to_expression_atlas_lucy_index.pl 

=head1 SYNOPSIS

    transform_phenotypes_to_expression_atlas_lucy_index.pl  -i [infile]
    
    perl bin/transform_phenotypes_to_expression_atlas_lucy_index.pl -i /home/vagrant/Downloads/cass_phenotype.csv -o /home/vagrant/cxgn/bea/bin/lucy.tsv -c /home/vagrant/cxgn/bea/bin/pre_corr.tsv -f /home/vagrant/cxgn/bea/bin/corr.tsv -p /home/vagrant/cxgn/bea/bin/project.txt -d /home/vagrant/cxgn/beambase/bin/desc.tsv -v 1 -n project_name

=head1 COMMAND-LINE OPTIONS
  ARGUMENTS
 -i phenotype csv file downloaded directly from website
 -o output lucy.tsv file that can be indexed using TEA script
 -c output pre-correlation file. This is the file fed into the R script below.
 -f output correlation.tsv file that can be indexed using TEA script
 -p output project.txt file that can be loaded into database using TEA script
 -d output description file that can be indexed by TEA script
 -v version number for how to group variables (e.g. grouping traits with tissues). currently 1 to 5
 -n project name

=head1 DESCRIPTION


=head1 AUTHOR

 Nicolas Morales (nm529@cornell.edu)

=cut

use strict;

use Getopt::Std;
use Data::Dumper;
use Carp qw /croak/ ;
use Pod::Usage;
use Spreadsheet::ParseExcel;
use Statistics::Basic qw(:all);
use Statistics::R;
use Text::CSV;
use JSON;
use LWP::Simple;

our ($opt_i, $opt_p, $opt_o, $opt_c, $opt_f, $opt_d, $opt_v, $opt_n, $opt_t, $opt_u, $opt_s);

getopts('i:p:o:c:f:d:v:n:t:u:s:');

if (!$opt_i || !$opt_p || !$opt_o || !$opt_c || !$opt_f || !$opt_d || !$opt_v || !$opt_n || !$opt_t || !$opt_u || !$opt_s) {
    die "Must provide options -i (input file) -p (project file out) -o (lucy out file) -c (corr pre-3col out file) -f (corr out file) -d (metabolite description oufile -v (script_version) -n (project name) -t (temp dir) -u (user_name) -s (search param description)\n";
}

my $csv = Text::CSV->new({ sep_char => ',' });

open(my $fh, '<', $opt_i)
    or die "Could not open file '$opt_i' $!";

#my $brapi_json_response = <$fh>;
#my $brapi_response = decode_json $brapi_json_response;
#print STDERR Dumper $brapi_response;
#my $remote_file_path = $brapi_response->{metadata}->{datafiles}->[0];
#print STDERR $remote_file_path."\n";
#my $phenotype_download = $opt_t."/phenotype_download.csv";
#my $status_code = mirror($remote_file_path, $phenotype_download);
#print STDERR $status_code."\n";


#open(my $fh, '<', $phenotype_download)
#    or die "Could not open file '$phenotype_download' $!";

my $trait_row = <$fh>;
my @columns;
if ($csv->parse($trait_row)) {
    @columns = $csv->fields();
} else {
    die "Could Not Parse Line: $trait_row\n";
}
my $col_max = scalar(@columns)-1;

#my $parser   = Spreadsheet::ParseExcel->new();
#my $excel_obj = $parser->parse($opt_i);

#my $worksheet = ( $excel_obj->worksheets() )[0]; #support only one worksheet
#my ( $row_min, $row_max ) = $worksheet->row_range();
#my ( $col_min, $col_max ) = $worksheet->col_range();

my @data_out;
my @traits;
for my $col ( 30 .. $col_max) {
    push @traits, $columns[$col];
}

#print STDERR Dumper \@traits;
#print STDERR scalar(@traits)."\n";

my %intermed;
my %corr_steps;
my %project_info;
my %accession_info_hash;
my $project_name = $opt_n;

my %unique_designs;
my %unique_locations;
my %unique_years;
my %unique_trial_names;

while ( my $row = <$fh> ){
    my @columns;
    if ($csv->parse($row)) {
        @columns = $csv->fields();
    } else {
        die "Could not parse line $row\n";
    }

    my $accession_name = $columns[18];
    my $trial_name = $columns[5];
    my $project_design = $columns[7];
    my $project_location = $columns[16];
    my $project_year = $columns[0];

    if (!exists($unique_designs{$project_design})){
        $unique_designs{$project_design} = 1;
        push @{$project_info{$project_name}->{designs}}, $project_design;
    }
    if (!exists($unique_locations{$project_location})){
        $unique_locations{$project_location} = 1;
        push @{$project_info{$project_name}->{locations}}, $project_location;
    }
    if (!exists($unique_years{$project_year})){
        $unique_years{$project_year} = 1;
        push @{$project_info{$project_name}->{years}}, $project_year;
    }
    if (!exists($unique_trial_names{$trial_name})){
        $unique_trial_names{$trial_name} = 1;
        push @{$project_info{$project_name}->{trial_names}}, $trial_name;
    }

    #print STDERR $accession_name."\n";
    my $plot_name = $columns[22];
    my $rep_number = $columns[23];
    my $block_number = $columns[24];
    my $plot_number = $columns[25];
    my $row_number = "Row_".$columns[26];
    my $col_number = "Col_".$columns[27];

    for( my $i=0; $i<scalar(@traits); $i++) {
        my $trait_col = $i + 30;
        #print STDERR "$row $trait_col\n";
        my $value = '';
        if ($columns[$trait_col]) {
            $value = $columns[$trait_col];
        }
        my $trait_term = $traits[$i];
        $trait_term =~ s/\s/_/g;
        $trait_term =~ s/ /_/g;
        $trait_term =~ s/\(//g;
        $trait_term =~ s/\)//g;
        $trait_term =~ s/\|/_/g;
        $trait_term =~ s/\://g;
        $trait_term =~ s/\///g;
        $trait_term =~ s/\-//g;

        my $stage;
        my $temp_key;
        my $step2;
        my $corr_step;
        if ($opt_v == 6){
            if (!$row_number) {
                print STDERR "No Row Number, skipping\n";
                next;
            }
            if (!$col_number) {
                print STDERR "No Col Number, skipping\n";
                next;
            }

            $temp_key = "$trait_term, $plot_name";
            $corr_step = "$plot_name, $row_number, $col_number";
            $accession_info_hash{$project_name}->{$col_number}->{$col_number}->{$row_number} = 1;

            if (exists($intermed{$temp_key})) {
                my $values = $intermed{$temp_key}->[3];
                push @$values, $value;
                $intermed{$temp_key}->[3] = $values;
            } else {
                $intermed{$temp_key} = [$trait_term, $row_number, $col_number, [$value], $corr_step];
            }
            $corr_steps{$corr_step} = 1;
        }
    }

}
#print STDERR Dumper \%intermed;
#print STDERR Dumper keys %intermed;
#print STDERR Dumper \%project_info;

foreach my $project_name (keys %accession_info_hash) {
    $project_info{$project_name}->{accessions} = $accession_info_hash{$project_name};
}

my @corr_steps_sorted;
foreach (sort keys %corr_steps) {
    push @corr_steps_sorted, $_;
}

my %corr_out;
my %unique_traits;
foreach (sort keys %intermed) {
    my $trait_term = $intermed{$_}->[0];
    my $row_number = $intermed{$_}->[1];
    my $col_number = $intermed{$_}->[2];
    my $values = $intermed{$_}->[3];
    my $corr_step = $intermed{$_}->[4];

    my @non_empty_values = grep($_ ne "", @$values);
    #print STDERR Dumper \@non_empty_values;
    my $average = mean(\@non_empty_values);
    my $stddev = stddev(\@non_empty_values);

    my $display_average = sprintf("%.2f", $average->query);
    my $display_stddev = sprintf("%.2f", $stddev->query);
    my @non_empty_values_formatted;
    foreach (@non_empty_values) {
        push @non_empty_values_formatted, sprintf("%.2f", $_);
    }

    push @data_out, [$trait_term, $row_number, $col_number, $display_average, $display_stddev, \@non_empty_values_formatted];

    $corr_out{$trait_term}->{$corr_step} = $display_average;
    $unique_traits{$trait_term}++;
}

#print STDERR Dumper \@data_out;
#print STDERR Dumper \%corr_out;

open(my $fh, ">", $opt_d);
    print STDERR $opt_d."\n";
    foreach (keys %unique_traits) {
        print $fh "1\t$_\t$_\n";
    }
close $fh;

open(my $fh, ">", $opt_o);
    print STDERR $opt_o."\n";
    foreach (@data_out) {
        my $values = $_->[5];
        my $trait_term = $_->[0];
        print $fh "$trait_term\t$_->[2]\t$_->[1]\t$_->[3]\t$_->[4]\t".join(',', @$values),"\n";
    }
close $fh;

open(my $fh, ">", $opt_c) || die("\nERROR:\n");
    print STDERR $opt_c."\n";
    print $fh "Metabolites\t", join("\t", @corr_steps_sorted), "\n";
    foreach my $trait_term (sort keys %corr_out) {
        print $fh "$trait_term\t";
        my $vals = $corr_out{$trait_term};
        my $step = 1;
        foreach my $corr_step (@corr_steps_sorted) {
            my $c = $vals->{$corr_step};
            if (!$c) { $c = 0;}
            print $fh "$c";
            if ($step < scalar(@corr_steps_sorted)) {
                print $fh "\t";
            }
            $step++;
        }
        print $fh "\n";
    }
close $fh;

my $R = Statistics::R->new();
my $out1 = $R->run(
    qq`data<-read.delim(exp <- "$opt_c", header=TRUE)`,
    qq`rownames(data)<-data[,1]`,
    qq`data<-data[,-1]`,
    qq`minexp <- 0`,
    qq`CorrMat<-as.data.frame(cor(t(data[rowSums(data)>minexp,]), method="spearman"))`,
    qq`write.table(CorrMat, file="$opt_f", sep="\t")`,
);

open (my $file_fh, "<", "$opt_f") || die ("\nERROR: the file $opt_f could not be found\n");
    my $header = <$file_fh>;
    chomp($header);
    my @gene_header = split("\t",$header);
    unshift(@gene_header,"");

    my %pairs;
    my @final_out;
    my $c = 1;
    while (my $line = <$file_fh>) {
        chomp($line);
        my @line = split("\t",$line);
        print STDERR Dumper \@line;

        for (my $n = 1; $n < $c; $n++) {
            $line[0] =~ s/\"//g;
            $gene_header[$n] =~ s/\"//g;

            my $hash_key = $line[0]."_".$gene_header[$n];

            if ($line[$n] && $line[0] ne $gene_header[$n]) {
                my $corr = $line[$n];
                if ($corr < 0) { $corr = $corr * -1; }
                push @final_out, "$line[0]\t$gene_header[$n]\t".sprintf("%.2f",$corr);
            }
        }
        $c++;
        #print STDERR $c."\n";
    }
close $file_fh;

open (my $file_fh, ">", "$opt_f") || die ("\nERROR:\n");
    print STDERR $opt_f."\n";
    foreach (@final_out) {
        print $file_fh $_;
        print $file_fh "\n";
    }
close $file_fh;

my $stage_ordinal = 1;
my $tissue_ordinal = 1;

open (my $file_fh, ">", "$opt_p") || die ("\nERROR:\n");
    print STDERR $opt_p."\n";
    print $file_fh "#organism\norganism_species: Manihot esculenta\norganism_variety: \norganism_description: Manihot esculenta\n# organism - end\n\n";
    foreach my $project_name (keys %project_info) {
        my $project_name_index_dir = $project_name;
        $project_name_index_dir =~ s/ //g;
        $project_name_index_dir =~ s/\s//g;
        my $project_years = $project_info{$project_name}->{years};
		my $project_years_text = join ',', @$project_years;
        my $project_designs = $project_info{$project_name}->{designs};
		my $project_designs_text = join ',', @$project_designs;
        my $project_locations = $project_info{$project_name}->{locations};
		my $project_locations_text = join ',', @$project_locations;
        my $trial_names = $project_info{$project_name}->{trial_names};
		my $trial_names_text = join ',', @$trial_names;

        print $file_fh "#project\nproject_name: $project_name\nproject_contact: $opt_u\nproject_description: $opt_s Returned Dataset Includes ( Trial(s) '$trial_names_text', Location(s): '$project_locations_text', Year(s): '$project_years_text', Breeding Program: 'CASS', Project Design(s): '$project_designs_text' )\nexpr_unit: varies\nindex_dir_name: cass_index_$project_name_index_dir\n# project - end\n\n";

        my $accession_hash = $project_info{$project_name}->{accessions};
        print STDERR Dumper $accession_hash;
        foreach my $accession (sort keys %$accession_hash) {

            print $file_fh "# figure --- All info needed for a cluster of images (usually includes a stage and all its tissues). Copy this block as many times as you need (including as many tissue layer blocks as you need).\nfigure_name: $accession\ncube_stage_name: $accession\nconditions:\n# write figure metadata\n\n";

            #if ($opt_v == 1){
                my $stage_hash = $accession_hash->{$accession};
                print STDERR Dumper $stage_hash;
                foreach my $stage (sort keys %$stage_hash) {
                    print $file_fh "#stage layer\nlayer_name: $accession\nlayer_description:\nlayer_type: stage\nbg_color:\nlayer_image: plant_background.png\nimage_width: 250\nimage_height: 500\ncube_ordinal: $stage_ordinal\nimg_ordinal: $stage_ordinal\norgan: $stage\n# layer - end\n\n";
                    $stage_ordinal++;

                    my $tissue_hash = $stage_hash->{$stage};
                    foreach my $tissue (sort keys %$tissue_hash) {
                        print $file_fh "#tissue layer\nlayer_name: $tissue\nlayer_description: $tissue\nlayer_type: tissue\nbg_color:\nlayer_image: $tissue.png\nimage_width: 250\nimage_height: 500\ncube_ordinal: $tissue_ordinal\nimg_ordinal: $tissue_ordinal\norgan: $stage\n# layer - end\n\n";
                        $tissue_ordinal++;
                    }

                }
            #}

            print $file_fh "# figure - end\n\n";
        }

    }
close $file_fh;

print STDERR "Script Complete.\n";
