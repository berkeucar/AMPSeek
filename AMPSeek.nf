 #! /usr/bin/env nextflow
nextflow.enable.dsl = 2

process DOWNLOADSEQUENCES {
    tag "Preparing for execution / Downloading requested data"
    publishDir "$params.output_path", mode: 'copy'
    label "cpu_process"

    output:
    path "*.{fasta,fa,fna}", emit: input_file

    script:
    if (params.download_from)
        """
        # --- Backup Logic ---
        if [ -d "$params.data_path" ]; then
            checksum=\$(find "$params.data_path" -type f -exec md5sum {} + | md5sum | cut -d ' ' -f1)
            backup_dir="${params.data_path}_backup_\${checksum}"
            if [ -d "\$backup_dir" ]; then
                echo "Backup already exists (\$backup_dir), skipping new backup."
            else
                echo "Backing up existing data to: \$backup_dir"
                mv "$params.data_path" "\$backup_dir"
            fi
        fi

        mkdir -p "$params.data_path"
        
        # --- Download Logic ---
        wget -O "downloaded_input.fa" "$params.download_from"
        
        # Copy to the external data path for storage
        cp "downloaded_input.fa" "$params.data_path/"
        """
    else
        """
        find "$params.data_path" -maxdepth 1 -type f -name "*.fasta" -exec ln -s {} ./ \\;
        find "$params.data_path" -maxdepth 1 -type f -name "*.fa" -exec ln -s {} ./ \\;
        find "$params.data_path" -maxdepth 1 -type f -name "*.fna" -exec ln -s {} ./ \\;

        # Count how many valid files we just linked into the working directory
        file_count=\$(ls *.fasta *.fa *.fna 2>/dev/null | wc -l)

        if [ "\$file_count" -eq 1 ]; then
            echo "Successfully found and linked exactly one sequence file from $params.data_path"
        elif [ "\$file_count" -gt 1 ]; then
            echo "Error: Found multiple sequence files (\${file_count}) in $params.data_path! Please ensure only one is present."
            exit 1
        else
            echo "Error: No download URL provided, and no existing fasta file found in $params.data_path!"
            exit 1
        fi
        """
}

process RUNAMPLIFY {
    tag "Running AMPlify"
    label "cpu_process"
    
    input:
    path data_path
    path output_path

    output:
    path "$output_path/*.tsv"

    script:
    """
    AMPlify -m balanced -s $data_path -od $output_path -sub on -att on
    """
}

process RUNCOLABFOLD{
    tag "Running colabfold"
    label workflow.profile.contains('gpu') ? "gpu_process" : "cpu_process" 
    
    input:
    path data_path
    path output_path

    output:
    path "$output_path/foldings"

    script:
    """
    mkdir -p ./tmp_home
    mkdir -p ./colabfold_cache
    mkdir -p ./matplotlib_config
    export HOME=\$PWD/tmp_home
    export COLABFOLD_CACHE_DIR=\$PWD/colabfold_cache
    export MPLCONFIGDIR=\$PWD/matplotlib_config

    colabfold_batch --amber --zip $data_path $output_path/foldings
    """
}

process RUNTAMPER {
    tag "Running tAMPer"
    label workflow.profile.contains('gpu') ? "gpu_process" : "cpu_process"
    
    input:
    path input_data
    path structure_data
    path output_path

    output:
    path "$output_path/results.csv"

    script:
    """
    python /opt/tAMPer/src/predict_tAMPer.py -seqs $input_data -pdbs $structure_data -chkpnt /opt/tAMPer/checkpoints/trained/chkpnt.pt -out $output_path
    find $structure_data -type f ! -name '*.zip' -delete
    """
}

process COMPILERESULTS{
    tag "Compiling results"
    label "cpu_process"

    input:
    path amplify
    path tamper
    path colabfold
    path compiler_path
    path imgs
    path templates

    script:
    if(params.output_file)
        """
        python $compiler_path $amplify $tamper $colabfold $params.output_path/$params.output_file $imgs $templates
        rm -f $params.output_path/AMPlify*.tsv
        rm -f $params.output_path/*.csv
        """
    else
        """
        python $compiler_path $amplify $tamper $colabfold $params.output_path/results.html $imgs $templates
        rm -f $params.output_path/AMPlify*.tsv
        rm -f $params.output_path/*.csv
        """
}

workflow{
    prep_out = DOWNLOADSEQUENCES()

    input_data_ch = prep_out.input_file
    output_data_ch = Channel.fromPath("$params.output_path")
    compiler_path = Channel.fromPath("$projectDir/src/make_report.py")
    template_path = Channel.fromPath("$projectDir/templates/report_template.html")
    img_path = Channel.fromPath("$projectDir/imgs/Logo.png")
    
    output_amplify = RUNAMPLIFY(input_data_ch, output_data_ch)
    output_colabfold = RUNCOLABFOLD(input_data_ch, output_data_ch)
    output_tamper = RUNTAMPER(input_data_ch, output_colabfold, output_data_ch)
    COMPILERESULTS(output_amplify, output_tamper, output_colabfold, compiler_path, img_path, template_path)
}

workflow.onComplete {
    log.info ( workflow.success ? "\nDone! The output is in --> $params.output_path\n" : "Oops .. something went wrong" )
}
