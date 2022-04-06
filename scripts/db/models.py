from enum import Enum as PyEnum

from sqlalchemy import (
    Boolean,
    Column,
    DateTime,
    Enum,
    FetchedValue,
    Float,
    ForeignKey,
    Integer,
    LargeBinary,
    String,
)
from sqlalchemy.orm import relationship
from sqlalchemy.ext.declarative import declarative_base

Base = declarative_base()
metadata = Base.metadata


class InputFileType(PyEnum):
    UNKNOWN = "UNKNOWN"
    BAM = "BAM"
    FASTQ = "FASTQ"

    def __str__(self):
        return str(self.name)


class Workflow(PyEnum):
    UNKNOWN = "UNKNOWN"
    ILLUMINA_ARTIC = "ILLUMINA_ARTIC"
    MEDAKA_ARTIC = "MEDAKA_ARTIC"

    def __str__(self):
        return str(self.name)


class PangolinStatus(PyEnum):
    UNKNOWN = "UNKNOWN"
    FAIL = "FAIL"
    PASS = "PASS"

    def __str__(self):
        return str(self.name)


class AnalysisRun(Base):  # type: ignore
    __tablename__ = "analysis_run"
    __table_args__ = {"schema": "psga", "comment": "The analysis run of samples"}

    analysis_run_id = Column(
        Integer,
        primary_key=True,
        server_default=FetchedValue(),
        comment="Primary key, serial sequence generated in the database",
    )
    analysis_run_name = Column(String, comment="The name of this run")
    primer_scheme_name = Column(String, comment="The primer scheme name")
    primer_scheme_version = Column(String, comment="The primer scheme version")
    input_file_type = Column(
        Enum(InputFileType, name="input_file_type"),
        nullable=False,
        server_default=FetchedValue(),
        comment="The type of input files",
    )
    workflow = Column(
        Enum(Workflow, name="workflow"),
        nullable=False,
        server_default=FetchedValue(),
        comment="The name of the workflow",
    )
    pipeline_version = Column(String, comment="PSGA pipeline version")
    pangolin_version = Column(String, comment="Pangolin version")
    pangolin_data_version = Column(
        String,
        comment="A version number that represents both pangolin-data version number, which as of pangolin 4.0 "
        "corresponds to the pango-designation version used to prepare the inference files. For example:"
        "(A) PANGO-1.2 indicates an identical sequence has been previously designated this lineage, and has so "
        "gone through manual curation. The number 1.2 indicates the version of pango-designation that this "
        "assignment is based on. These hashes and pango-designation version is bundled with the pangoLEARN "
        "and UShER models. "
        "(B) PLEARN-1.2 indicates that this sequence is different from any previously designated and that "
        "the pangoLEARN model was used as an inference engine to predict the most likely lineage based on "
        "the given version of pango-designation upon which the pangoLEARN model was trained. (C) PUSHER-1.2 "
        "indicates that this sequence is different from any previously designated and that UShER was used as "
        "an inference engine with fast tree placement and parsimony-based lineage assignment, based on a guide "
        "tree (protobuf) file built from the data in a given pango-designation release version.",
    )
    constellation_version = Column(
        String, comment="The version of constellations that scorpio has used to curate the lineage assignment."
    )
    scorpio_version = Column(String, comment="Scorpio version used by Pangolin")


class Sample(Base):  # type: ignore
    __tablename__ = "sample"
    __table_args__ = {"schema": "psga", "comment": "Sample data"}

    sample_id = Column(
        Integer,
        primary_key=True,
        server_default=FetchedValue(),
        comment="Primary key, serial sequence generated in the database",
    )
    analysis_run_id = Column(
        Integer,
        ForeignKey("psga.analysis_run.analysis_run_id", deferrable=True, initially="DEFERRED"),
        comment="Foreign key to analysis_run table",
    )
    sample_name = Column(String, comment="Lab sample identifier")
    date_collected = Column(DateTime, comment="timestamp of sample collection")
    pangolin_lineage = Column(
        String,
        comment="The most likely lineage assigned to a given sequence based on the inference engine used and the "
        "SARS-CoV-2 diversity designated. This assignment may be is sensitive to missing data at key sites.",
    )
    pangolin_conflict = Column(
        Float,
        comment="In the pangoLEARN model, a given sequence gets assigned to the most likely category based on known "
        "diversity. If a sequence can fit into more than one category, the conflict score will be greater "
        "than 0 and reflect the number of categories the sequence could fit into. If the conflict score is 0, "
        "this means that within the current decision tree there is only one category that the sequence could "
        "be assigned to.",
    )
    pangolin_ambiguity_score = Column(
        Float,
        comment="This score is a function of the quantity of missing data in a sequence. It represents the proportion "
        "of relevant sites in a sequnece which were imputed to the reference values. A score of 1 indicates "
        "that no sites were imputed, while a score of 0 indicates that more sites were imputed than were not "
        "imputed. This score only includes sites which are used by the decision tree to classify a sequence.",
    )
    pangolin_status = Column(
        Enum(PangolinStatus, name="pangolin_status"),
        nullable=False,
        server_default=FetchedValue(),
        comment="Indicates whether the sequence passed the QC thresholds for minimum length and maximum N content.",
    )
    scorpio_call = Column(
        String,
        comment="If a query is assigned a constellation by scorpio this call is output in this column. "
        "The full set of constellations searched by default can be found at the constellations repository.",
    )
    scorpio_support = Column(
        Float,
        comment=(
            "The support score is the proportion of defining variants which have the alternative "
            "allele in the sequence.",
        ),
    )
    scorpio_conflict = Column(
        Float,
        comment=(
            "The conflict score is the proportion of defining variants which have the reference allele in the "
            "sequence. Ambiguous/other non-ref/alt bases at each of the variant positions contribute only to the "
            "denominators of these scores.",
        ),
    )
    scorpio_notes = Column(
        String,
        comment=(
            "If any conflicts from the decision tree, this field will output the alternative assignments. "
            "If the sequence failed QC this field will describe why. If the sequence met the SNP thresholds "
            "for scorpio to call a constellation, it’ll describe the exact SNP counts of Alt, Ref and Amb "
            "(Alternative, reference and ambiguous) alleles for that call.",
        ),
    )
    is_designated = Column(
        Boolean,
        comment=(
            "A boolean (True/False) column indicating whether that particular sequence has been offically "
            "designated a lineage."
        ),
    )
    qc_notes = Column(
        String,
        comment=("Notes specific to the QC checks run on the sequences."),
    )
    note = Column(
        String,
        comment=(
            "If any conflicts from the decision tree, this field will output the alternative assignments. "
            "If the sequence failed QC this field will describe why. If the sequence met the SNP thresholds "
            "for scorpio to call a constellation, it’ll describe the exact SNP counts of Alt, Ref and Amb "
            "(Alternative, reference and ambiguous) alleles for that call."
        ),
    )
    genbank_submit_id = Column(
        String, comment="Unique identifier of GenBank submission " "session, which was used to submit a sample."
    )
    gisaid_id = Column(
        String,
        comment=(
            "GISAID identifier (GISAID is a global science initiative and primary source that provides open-access "
            "to genomic data of the novel coronavirus responsible for COVID-19)",
        ),
    )
    genome_length = Column(Integer, comment="Number of base pair in the virus genome")
    metadata_loaded = Column(Boolean, comment="Metadata was loaded")

    sample_qc = relationship("SampleQC", uselist=False)


class SampleQC(Base):  # type: ignore
    __tablename__ = "sample_qc"
    __table_args__ = {
        "schema": "psga",
        "comment": "Sample Quality Control, QC data for a sample",
    }

    sample_id = Column(
        Integer,
        ForeignKey("psga.sample.sample_id", deferrable=True, initially="DEFERRED"),
        primary_key=True,
        comment="Primary key, and foreign key to sample table",
    )
    pct_n_bases = Column(Float(53), comment="percentage of N bases in sequenced genome")
    pct_covered_bases = Column(Float(53), comment="percentage of covered bases in sequenced genome")
    longest_no_n_run = Column(Integer, comment="longest sequence run without N base")
    num_aligned_reads = Column(Float(53), comment="Number of reads succesfully aligned and mapped")
    qc_pass = Column(Boolean, comment="Has sample passed QC")
    qc_plot = Column(LargeBinary, comment="QC plot image")
