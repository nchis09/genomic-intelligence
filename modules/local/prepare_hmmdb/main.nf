/*
 * Local module: PREPARE_HMMDB
 *
 * Download a generic HMM database once, decompress/concatenate the
 * individual *.hmm profiles, and press the result with hmmpress so it
 * can be re-used by hmmscan across pipeline runs.
 *
 * The output directory is stored in params.hmm_db_cache_dir via storeDir
 * in conf/modules.config.
 */

process PREPARE_HMMDB {
    tag 'hmmdb'
    label 'process_low'

    conda "${projectDir}/envs/pgirl_hmmer.yml"
    container null

    input:
    val hmm_db_url

    output:
    path "hmmdb", emit: hmm_db_dir

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mkdir -p hmmdb/extracted
    wget -q -O hmmdb/hmm_download ${hmm_db_url}
    if [[ "${hmm_db_url}" == *.tar.gz ]]; then
        tar -xzf hmmdb/hmm_download -C hmmdb/extracted
        find hmmdb/extracted -type f -name "*.hmm" -print0 | xargs -0 cat > hmmdb/hmm_db.hmm
    elif [[ "${hmm_db_url}" == *.hmm.gz ]]; then
        gunzip -c hmmdb/hmm_download > hmmdb/hmm_db.hmm
    else
        mv hmmdb/hmm_download hmmdb/hmm_db.hmm
    fi

    rm -rf hmmdb/extracted hmmdb/hmm_download
    hmmpress hmmdb/hmm_db.hmm
    """
}
