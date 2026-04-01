"""Quick smoke-test: digest a test PEFF file and print the resulting DataFrame."""
import polars as pl

from peff_digest import DigestConfig, digest, InternalMod

if __name__ == "__main__":


    internal_mods = [
        # InternalModification(residue="P", modification="+0.0001", mod_type=ModType.STATIC),
        InternalMod(residue="M", modification="U:Oxidation", mod_type="variable"),
        InternalMod(residue="K", modification="U:Propyl", mod_type="variable"),
        InternalMod(residue="K", modification="U:Butyryl", mod_type="variable"),
    ]
    internal_mods = []

    config = DigestConfig(
        input_file="/home/patrick-garrett/Repos/peff_uniprot_fetcher/p01_argc.peff",
        cleave_on="R",
        missed_cleavages=1,
        min_length=4,
        max_length=50,
        max_ptm_per_peptide=2,
        use_psi_mods=True,
        use_mod_names=True,
        use_unimod_output=True,
        workers=10,
        consolidate_proteins=True,
        internal_mods=internal_mods,
        drop_invalid_mass=True,
    )

    df = digest(config, show_progress=True)
    assert df is not None
    print(df)
    print(f"\n{len(df)} peptides from {df['protein_id'].n_unique()} proteins")

    # save to parquet for easier inspection in Python or R
    df.write_parquet("peptides.parquet")

    # convert mass_array to string and save to csv
    df.with_columns(
        pl.col("mass_array").list.eval(pl.element().cast(pl.String)).list.join(" ").alias("mass_array_str")
    ).drop("mass_array").write_csv("peptides.csv")
