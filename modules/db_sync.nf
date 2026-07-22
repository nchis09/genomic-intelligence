process DB_SYNC {
    tag "db_update"

    input:
    val db_url
    val pathogens
    val sources
    val dry_run

    output:
    stdout

    script:
    def dry_flag = dry_run ? "--dry-run" : ""
    """
    export PYTHONPATH="${projectDir}"
    ${params.python} "${projectDir}/scripts/db_sync.py" \
        --db-url "${db_url}" \
        --pathogens "${pathogens}" \
        --sources "${sources}" \
        ${dry_flag}
    """
}
