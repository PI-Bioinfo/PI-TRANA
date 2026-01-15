include { DIAMOND_MAKEDB } from '../../modules/nf-core/diamond/makedb/main'            
include { DIAMOND_BLASTX } from '../../modules/nf-core/diamond/blastx/main'          
include { DIAMOND_FILTER } from '../../modules/local/filter_amr'

workflow AMR_GENES {

    take:
    ch_processed_optionally_sampled_reads
    fasta_database
    file_taxonmap
    file_taxonnodes
    file_taxonnames
    
    main:
    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()

    //
    // MODULE: Run DIAMOND_MAKED if needed
    //
    if (params.diamond_index == null) {
    DIAMOND_MAKEDB([channel.empty(),fasta_database],
        file_taxonmap,
        file_taxonnodes,
        file_taxonnames)  
    diamond_db = DIAMOND_MAKEDB.out.db
    } else {
        diamond_db = params.diamond_index
    }
    ch_diamond_db = Channel.fromPath(params.diamond_index)
        .map { db -> [[:], db] }

    DIAMOND_BLASTX(ch_processed_optionally_sampled_reads,
                ch_diamond_db,
                'txt',
                '')

    DIAMOND_FILTER(DIAMOND_BLASTX.out.txt)
}