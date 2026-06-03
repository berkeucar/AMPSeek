#! /usr/bin/env nextflow
nextflow.enable.dsl = 2

process DOWNLOADSEQUENCES {
    tag "Preparing for execution / Downloading requested data"
    label "cpu_process"
    def data_path = file(params.data_path)
    def is_file_target = data_path.exists() ? data_path.isFile() : params.data_path.matches('.*\\.(fasta|fa|fna)$')

    output:
    path( is_file_target ? "${data_path.name}" : "*.{fasta,fa,fna}" )

    script:
    if (params.download_from) {
        if (is_file_target)
        """
        # --- Backup Logic ---
        target_file="${params.data_path}"

        if [ -f "\$target_file" ]; then
            checksum=\$(md5sum "\$target_file" | cut -d ' ' -f1)
            base_path="\${target_file%.*}"
            ext="\${target_file##*.}"
            
            if [ "\$target_file" == "\$ext" ] || [ "\$base_path" == "\$target_file" ]; then
                backup_file="\${target_file}_backup_\${checksum}"
            else
                backup_file="\${base_path}_backup_\${checksum}.\${ext}"
            fi
            
            if [ -f "\$backup_file" ]; then
                echo "Backup already exists (\$backup_file), skipping new backup."
            else
                echo "Backing up existing data to: \$backup_file"
                mv "\$target_file" "\$backup_file"
            fi
        fi

        wget -O "\$target_file" "${params.download_from}"
        ln -s "\$target_file" ./
        """
        else 
        """
        # --- Backup Logic ---
        if [ -d "${params.data_path}" ]; then
            checksum=\$(find "${params.data_path}" -type f -exec md5sum {} + | md5sum | cut -d ' ' -f1)
            backup_dir="${params.data_path}_backup_\${checksum}"
            if [ -d "\$backup_dir" ]; then
                echo "Backup already exists (\$backup_dir), skipping new backup."
            else
                echo "Backing up existing data to: \$backup_dir"
                mv "${params.data_path}" "\$backup_dir"
            fi
        fi

        mkdir -p "${params.data_path}"
        wget -O "${params.data_path}/downloaded_input.fa" "${params.download_from}"
        ln -s "${params.data_path}/downloaded_input.fa" ./ 
        """
    }
    else {
        if (is_file_target)
        """
        target_file="${params.data_path}"
        
        if [ -f "\$target_file" ]; then
            echo "Successfully found exactly one sequence file named ${params.data_path}"
            ln -s "\$target_file" ./
        else 
            echo "Error: There is no file named ${params.data_path}"
            exit 1
        fi 
        """
        else
        """
        find "${params.data_path}" -maxdepth 1 -type f -name "*.fasta" -exec ln -s {} ./ \\;
        find "${params.data_path}" -maxdepth 1 -type f -name "*.fa" -exec ln -s {} ./ \\;
        find "${params.data_path}" -maxdepth 1 -type f -name "*.fna" -exec ln -s {} ./ \\;

        file_count=\$(ls *.fasta *.fa *.fna 2>/dev/null | wc -l)

        if [ "\$file_count" -eq 1 ]; then
            echo "Successfully found and linked exactly one sequence file from ${params.data_path}"
        elif [ "\$file_count" -gt 1 ]; then
            echo "Error: Found multiple sequence files (\${file_count}) in ${params.data_path}! Please ensure only one is present."
            exit 1
        else
            echo "Error: No download URL provided, and no existing fasta file found in ${params.data_path}!"
            exit 1
        fi
        """
    }
}

process RUNAMPLIFY {
    tag "Running AMPlify"
    label "cpu_process"
    
    input:
    path data_path

    output:
    path "*.tsv"

    script:
    """
    AMPlify -m balanced -s $data_path -sub on -att on
    """
}

process RUNCOLABFOLD{
    tag "Running colabfold"
    label workflow.profile.contains('gpu') ? "gpu_process" : "cpu_process" 
    publishDir "${params.output_path}", mode: 'copy'

    input:
    path data_path

    output:
    path "foldings"

    script:
    def jax_env = workflow.profile.contains('gpu') ? "" : "export JAX_PLATFORMS=cpu"
    """
    $jax_env
    mkdir -p ./tmp_home
    mkdir -p ./colabfold_cache
    mkdir -p ./matplotlib_config
    export HOME=\$PWD/tmp_home
    export COLABFOLD_CACHE_DIR=\$PWD/colabfold_cache
    export MPLCONFIGDIR=\$PWD/matplotlib_config

    colabfold_batch --amber --zip $data_path foldings
    cd foldings
    find . -maxdepth 1 -type f ! -name '*.zip' -delete
    rm -rf */
    """
}

process RUNTAMPER {
    tag "Running tAMPer"
    label workflow.profile.contains('gpu') ? "gpu_process" : "cpu_process"
    
    input:
    path input_data
    path structure_data

    output:
    path "results.csv"

    script:
    """
    export MPLCONFIGDIR=\$PWD/matplotlib_config
    python /opt/tAMPer/src/predict_tAMPer.py -seqs $input_data -pdbs $structure_data -chkpnt /opt/tAMPer/checkpoints/trained/chkpnt.pt -out .
    """
}

process COMPILERESULTS{
    tag "Compiling results"
    label "cpu_process"

    publishDir "${params.output_path}", mode: 'copy'

    input:
    path amplify
    path tamper
    path colabfold
    path compiler_path
    path imgs
    path templates

    output:
    path "*.html"

    script:
    output_name = params.output_file ?: "results.html"
    """
    export MPLCONFIGDIR=\$(pwd)/matplotlib_config
    python $compiler_path $amplify $tamper $colabfold $output_name $imgs $templates
    """
}

// create the output directory if necessary before the pipeline starts
def myDir = file(params.output_path)
if (!myDir.exists()) {
    myDir.mkdirs()
    println "Created output directory: ${myDir.name}"
}

// create the log directory if running on slurm
if (workflow.profile.contains('slurm') && !file(params.slurm_log_dir).exists()){
    file(params.slurm_log_dir).mkdirs()
    println "Created SLURM log directory: ${params.slurm_log_dir}"
}

workflow{
    input_data_ch = DOWNLOADSEQUENCES()
    output_data_ch = Channel.fromPath("$params.output_path")
    compiler_path = Channel.fromPath("$projectDir/src/make_report.py")
    template_path = Channel.fromPath("$projectDir/templates/report_template.html")
    img_path = Channel.fromPath("$projectDir/imgs/Logo.png")
    
    output_amplify = RUNAMPLIFY(input_data_ch)
    output_colabfold = RUNCOLABFOLD(input_data_ch)
    output_tamper = RUNTAMPER(input_data_ch, output_colabfold)
    COMPILERESULTS(output_amplify, output_tamper, output_colabfold, compiler_path, img_path, template_path)
}

workflow.onComplete {
    log.info ( workflow.success ? "\nDone! The output is in --> $params.output_path\n" : "Oops .. something went wrong" )
}
