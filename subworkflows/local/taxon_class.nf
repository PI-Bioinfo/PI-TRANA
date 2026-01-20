include { TRANSLATE_TAXIDS                       } from '../../modules/local/translate_taxids/main.nf'
include { KRAKEN2_KRAKEN2                        } from '../../modules/nf-core/kraken2/kraken2'
include { ASSIGNMENT_HEATMAP                     } from '../../modules/local/assignment_heatmap/main.nf'
include { EMU_ABUNDANCE                          } from '../../modules/local/emu/abundance/main.nf'
include { EMU_COMBINE_OUTPUTS                    } from '../../modules/local/emu/combine_outputs/main.nf'
include { KRONA_KTIMPORTTAXONOMY                 } from '../../modules/nf-core/krona/ktimporttaxonomy/main.nf'

workflow TAXONOMIC_CLASSIFICATION {
    take:
        ch_processed_optionally_sampled_reads

    main:
        ch_versions = Channel.empty()
        //
        // MODULE: Run Kraken2 
        if (params.run_kraken2) {
            KRAKEN2_KRAKEN2(
                ch_processed_optionally_sampled_reads,
                params.kraken2_db,
                params.kraken2_save_output_fastqs,
                params.kraken2_save_read_assignment
            )
        }
        //
        // MODULE: run EMU abundance calculation
        EMU_ABUNDANCE(ch_processed_optionally_sampled_reads)
        ch_versions = ch_versions.mix(EMU_ABUNDANCE.out.versions)
        //
        // MODULE: run emu combine-outputs
        ch_emu_combine_input_files = channel.empty()
        ch_emu_combine_input_files = EMU_ABUNDANCE.out.report
            .map { it[1] }  // Extract only the file path from the tuple (meta, path)
            .collect()
            .set { collected_files }
        EMU_COMBINE_OUTPUTS(collected_files)
        //
        // MODULE: Run krona plot
        if (params.run_krona) {
            KRONA_KTIMPORTTAXONOMY(
                EMU_ABUNDANCE.out.report,
                file(params.krona_taxonomy_tab, checkExists: true)
            )
            ch_versions = ch_versions.mix(KRONA_KTIMPORTTAXONOMY.out.versions.first())
        }

        //
        // MODULE: run translate_taxids
        if (params.make_heatmap) {
            TRANSLATE_TAXIDS(EMU_ABUNDANCE.out.assignment_report)
            ch_versions = ch_versions.mix(TRANSLATE_TAXIDS.out.versions)

            //
            // Module: run assignment_heatmap
            ASSIGNMENT_HEATMAP(TRANSLATE_TAXIDS.out.assignment_translated_report)
            ch_versions = ch_versions.mix(ASSIGNMENT_HEATMAP.out.versions)
        }
        ch_versions = ch_versions.mix(EMU_COMBINE_OUTPUTS.out.versions)

    emit:
        report                  = ch_emu_combine_input_files
        combined_report         = EMU_COMBINE_OUTPUTS.out.combined_report
        combined_counts_report  = EMU_COMBINE_OUTPUTS.out.combined_counts_report
        versions                = ch_versions
}