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
    Table,
    Text,
)
from sqlalchemy.orm import relationship
from sqlalchemy.ext.declarative import declarative_base

Base = declarative_base()
metadata = Base.metadata


class GovernorateEnum(PyEnum):
    Capital = "Capital"
    Northern = "Northern"
    Southern = "Southern"
    Muharraq = "Muharraq"

    def __str__(self):
        return str(self.name)


class GenderEnum(PyEnum):
    M = "M"
    F = "F"
    U = "U"

    def __str__(self):
        return str(self.name)


class HospitalAdmittanceEnum(PyEnum):
    No = "No"
    Standard = "Standard"
    ICU = "ICU"

    def __str__(self):
        return str(self.name)


class PangolinStatus(PyEnum):
    unknown = "unknown"
    fail = "fail"
    passed_qc = "passed_qc"

    def __str__(self):
        return str(self.name)


class Area(Base):  # type: ignore
    __tablename__ = "area"
    __table_args__ = {"schema": "sars_cov_2", "comment": 'Geographic "Area" location'}

    name = Column(
        String,
        primary_key=True,
        comment="Primary key, the area name, capitol level geolocation from where the sample was collected",
    )
    latitude = Column(
        Float(53),
        comment="The angular distance of the location north or south of the equator",
    )
    longitude = Column(
        Float(53),
        comment="the angular distance of the location east or west of the Greenwich meridian",
    )
    colour = Column(LargeBinary, comment="The colour of this block in reports")


class Block(Base):  # type: ignore
    __tablename__ = "block"
    __table_args__ = {"schema": "sars_cov_2", "comment": 'Geographic "Block" location'}

    number = Column(
        String,
        primary_key=True,
        comment="Primary key, the block number, low level geolocation from where the sample was collected",
    )
    latitude = Column(
        Float(53),
        comment="The angular distance of the location north or south of the equator",
    )
    longitude = Column(
        Float(53),
        comment="the angular distance of the location east or west of the Greenwich meridian",
    )
    colour = Column(LargeBinary, comment="The colour of this block in reports")


class Comorbidity(Base):  # type: ignore
    __tablename__ = "comorbidity"
    __table_args__ = {
        "schema": "sars_cov_2",
        "comment": "Coded diseases presented in patient (eg Diabetes)",
    }

    comorbidity_id = Column(Integer, primary_key=True, comment="Primary key")
    description = Column(Text, comment="Description of the comorbidity")

    samples = relationship("Sample", secondary="sars_cov_2.sample_comorbidity")


class Governorate(Base):  # type: ignore
    __tablename__ = "governorate"
    __table_args__ = {
        "schema": "sars_cov_2",
        "comment": 'Geographic "Governorate" location',
    }

    name = Column(
        Enum(GovernorateEnum, name="governorate_name"),
        primary_key=True,
    )
    latitude = Column(Float(53))
    longitude = Column(Float(53))
    colour = Column(LargeBinary)


class Sample(Base):  # type: ignore
    __tablename__ = "sample"
    __table_args__ = {"schema": "sars_cov_2", "comment": "Sample data"}

    sample_id = Column(
        Integer,
        primary_key=True,
        server_default=FetchedValue(),
        comment="Primary key, serial sequence generated in the database",
    )
    lab_id = Column(String, comment="Lab sample identifier")
    date_collected = Column(DateTime, comment="timestamp of sample collection")
    data_sequenced = Column(DateTime, comment="timestamp of sample sequencing")
    sequencing_run = Column(String, comment="identifier of sequencing run")
    gender = Column(
        Enum(GenderEnum, name="gender"),
        nullable=False,
        server_default=FetchedValue(),
        comment="Gender of host person",
    )
    age = Column(Integer, comment="Age of host person in years")
    nationality = Column(String, comment="Nationality of host person")
    governorate_name = Column(
        ForeignKey("sars_cov_2.governorate.name", deferrable=True, initially="DEFERRED"),
        comment="high level geolocation from where sample was collected",
    )
    area_name = Column(
        ForeignKey("sars_cov_2.area.name", deferrable=True, initially="DEFERRED"),
        comment="foreign key to area table, mid level geolocation from where sample was collected",
    )
    block_number = Column(String, comment="low level geolocation from where sample was collected")
    sample_number = Column(Integer, comment="Lab sample identifier (may be redundant with lab_id)")
    ct_value = Column(Float(53), comment="RNA quality measure")
    symptoms = Column(Text, comment="Symptoms observed")
    travel_exposure = Column(String, comment="Country host visited recently")
    hospital_admittance = Column(
        Enum(HospitalAdmittanceEnum, name="hospital_admittance"),
        comment="Hospital type host admited to",
    )
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
    mrn = Column(String, comment="Individual National Identity number")
    metadata_loaded = Column(Boolean, comment="Metadata was loaded from Bahrain medical records system I-SEHA")

    area = relationship("Area")
    governorate = relationship("Governorate")
    comorbidities = relationship("Comorbidity", secondary="sars_cov_2.sample_comorbidity")
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
    pipeline_version = Column(String, comment="mapping pipeline version")
    pangolearn_version = Column(String, comment="pangoLEARN version used by Pangolin")


t_sample_comorbidity = Table(
    "sample_comorbidity",
    metadata,
    Column(
        "sample_id",
        ForeignKey("sars_cov_2.sample.sample_id", deferrable=True, initially="DEFERRED"),
        primary_key=True,
        nullable=False,
        comment="Primary key with comorbidity_id, and foreign key to sample table",
    ),
    Column(
        "comorbidity_id",
        ForeignKey(
            "sars_cov_2.comorbidity.comorbidity_id",
            deferrable=True,
            initially="DEFERRED",
        ),
        primary_key=True,
        nullable=False,
        comment="Primary key with sample_id, and foreign key to comorbidity table",
    ),
    schema="sars_cov_2",
    comment="Sample comorbitiy, a linking table between the sample an comorbidity record",
)
