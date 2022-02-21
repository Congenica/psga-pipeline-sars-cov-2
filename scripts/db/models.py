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
    unknown = "unknown"
    bam = "bam"
    fastq = "fastq"

    def __str__(self):
        return str(self.name)


class Workflow(PyEnum):
    unknown = "unknown"
    illumina_artic = "illumina_artic"
    medaka_artic = "medaka_artic"

    def __str__(self):
        return str(self.name)


class PangolinStatus(PyEnum):
    unknown = "unknown"
    fail = "fail"
    passed_qc = "passed_qc"

    def __str__(self):
        return str(self.name)


class AnalysisRun(Base):  # type: ignore
    __tablename__ = "analysis_run"
    __table_args__ = {"schema": "sars_cov_2", "comment": "The analysis run of samples"}

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
    pipeline_version = Column(String, comment="mapping pipeline version")
    pangolin_version = Column(String, comment="pangolin version")
    pangolearn_version = Column(String, comment="pangoLEARN version used by Pangolin")
    pango_version = Column(String, comment="pango version used by Pangolin")


class Sample(Base):  # type: ignore
    __tablename__ = "sample"
    __table_args__ = {"schema": "sars_cov_2", "comment": "Sample data"}

    sample_id = Column(
        Integer,
        primary_key=True,
        server_default=FetchedValue(),
        comment="Primary key, serial sequence generated in the database",
    )
    analysis_run_id = Column(
        Integer,
        ForeignKey("sars_cov_2.analysis_run.analysis_run_id", deferrable=True, initially="DEFERRED"),
        comment="Foreign key to analysis_run table",
    )
    sample_name = Column(String, comment="Lab sample identifier")
    date_collected = Column(DateTime, comment="timestamp of sample collection")
    pangolin_lineage = Column(String, comment="Viral phylogenetic lineage")
    pangolin_conflict = Column(
        Float, comment="Positive if a sequence can fit into more than 1 category based on known diversity"
    )
    pangolin_ambiguity_score = Column(Float, comment="A function of the quantity of missing data in a sequence")
    pangolin_status = Column(
        Enum(PangolinStatus, name="pangolin_status"),
        nullable=False,
        server_default=FetchedValue(),
        comment="Reported pangolin lineage status",
    )
    scorpio_call = Column(
        String, comment="If a query is assigned a constellation by scorpio this call is output in this column"
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
            "The conflict score is the proportion of defining variants which have the reference allele "
            "in the sequence.",
        ),
    )
    note = Column(
        String,
        comment=(
            "If any conflicts from the decision tree, this field will output the alternative assignments. "
            "If the sequence failed QC this field will describe why. If the sequence met the SNP thresholds "
            "for scorpio to call a constellation, it’ll describe the exact SNP counts of Alt, Ref and Amb "
            "(Alternative, reference and ambiguous) alleles for that call.",
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
    metadata_loaded = Column(Boolean, comment="Metadata was loaded from Bahrain medical records system I-SEHA")

    sample_qc = relationship("SampleQC", uselist=False)


class SampleQC(Base):  # type: ignore
    __tablename__ = "sample_qc"
    __table_args__ = {
        "schema": "sars_cov_2",
        "comment": "Sample Quality Control, QC data for a sample",
    }

    sample_id = Column(
        Integer,
        ForeignKey("sars_cov_2.sample.sample_id", deferrable=True, initially="DEFERRED"),
        primary_key=True,
        comment="Primary key, and foreign key to sample table",
    )
    pct_n_bases = Column(Float(53), comment="percentage of N bases in sequenced genome")
    pct_covered_bases = Column(Float(53), comment="percentage of covered bases in sequenced genome")
    longest_no_n_run = Column(Integer, comment="longest sequence run without N base")
    num_aligned_reads = Column(Float(53), comment="Number of reads succesfully aligned and mapped")
    qc_pass = Column(Boolean, comment="Has sample passed QC")
    qc_plot = Column(LargeBinary, comment="QC plot image")
