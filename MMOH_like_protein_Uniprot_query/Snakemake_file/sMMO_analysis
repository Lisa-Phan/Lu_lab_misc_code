"""
Creating a sample snake file for pocket analysis


Sample script used by author for pocket param characterization
$ROSETTA_DEV/pocket_measure.linuxgccrelease -s 
${pdb_id}.pdb -central_relax_pdb_num ${rosetta_res} -pocket_num_angles 100 -
pocket_dump_pdbs -pocket_filter_by_exemplar 1 -pocket_grid_size 15 -
ignore_unrecognized_res

"""
rule get_files:
### Create a directory of post-filtered pdb files
### Custom function is a shell script with download function 
### which takes an input file of AF code and downloading them
    input: "snakemake_sMMO/Filtered_small_set/AF_id_files.txt"
    output: "snakemake_sMMO/Alphafold_structure/Filtered_AF_files"
    shell: "source /home/09069/dhp563/custom_func.sh; download {input}"


rule create_database:
    input: 
    output: 


rule foldseek_get_superimposition:
    input:
        ""
    output:


rule rosetta_pocket: 
    input:
        ""
    output:
    shell:
        /work/09069/dhp563/ls6/rosetta_bin/rosetta_bin_linux_2021.16.61629_bundle/main/source/bin/pocket_measure.linuxgccrelease / 
        -s {pdb} 
        -central_relax_pdb_num {rosetta_res} 
        -pocket_num_angles 100 
        -