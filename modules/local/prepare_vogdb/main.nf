/*
 * Local module: PREPARE_VOGDB
 *
 * Download the full VOGDB HMM tarball once, extract/concatenate the
 * individual *.hmm profiles, and press the result with hmmpress so it
 * can be re-used by hmmscan across pipeline runs.
 *
 * The output directory is stored in params.vogdb_cache_dir via storeDir
 * in conf/modules.config.
 */

process PREPARE_VOGDB {
    tag 'vogdb'
    label 'process_low'

    conda "${projectDir}/envs/pgirl_hmmer.yml"
    container null

    input:
    val vogdb_url

    output:
    path "vogdb", emit: vogdb_dir

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mkdir -p vogdb/extracted
    wget -q -O vogdb/vog.hmm.tar.gz ${vogdb_url}
    tar -xzf vogdb/vog.hmm.tar.gz -C vogdb/extracted
    find vogdb/extracted -type f -name "*.hmm" -print0 | xargs -0 cat > vogdb/vog.hmm
    rm -rf vogdb/extracted vogdb/vog.hmm.tar.gz
    hmmpress vogdb/vog.hmm
    """
}
