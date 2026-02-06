nextflow.enable.dsl = 2

process VERSIONS {
    publishDir "${params.results_dir}", mode: 'copy'

    output:
    path "software_version.yml"

    script:
    """
    echo "pipeline:" > software_version.yml
    echo "  name: tb_mutation_analysis_mini" >> software_version.yml
    echo "  description: Deeplex and VCF module" >> software_version.yml
    echo "  URL: https://github.com/oucru-id/tb-to-fhir-mini" >> software_version.yml
    echo "  version: ${params.version}" >> software_version.yml
    echo "  nextflow: $nextflow.version" >> software_version.yml
    
    echo "databases:" >> software_version.yml
    echo "  reference: \$(basename ${params.reference})" >> software_version.yml
    echo "  mutation_db: \$(basename ${params.mutation_db})" >> software_version.yml
    echo "  repetitive_regions: \$(basename ${params.repetitive_regions})" >> software_version.yml
    
    echo "processing_settings:" >> software_version.yml
    echo "  vcf:" >> software_version.yml
    echo "    filter_min_depth: 5" >> software_version.yml
    echo "    filter_min_gq: 20" >> software_version.yml
    
    export BASE_DIR="${baseDir}"
    python3 $baseDir/scripts/get_versions.py >> software_version.yml
    """
}